//! C ABI used by longcallD.
//!
//! longcallD stores bases in the usual 4-bit encoding (`0=A, 1=C, 2=G, 3=T,
//! 4=N`).  This module is a thin, panic-safe shim that decodes those into the
//! ASCII sequences the POA engine expects, runs a single- or multi-allele
//! consensus, and encodes the result back.  The C side owns the output buffer
//! and must release it with [`poa_consensus_output_free`].
//!
//! The ABI deliberately exposes only plain arrays and `usize`/integers so it
//! can be called from C/C++ without bindgen.

#![allow(unsafe_op_in_unsafe_fn)]

use std::ffi::{CString, c_char};
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::ptr;
use std::slice;

use crate::{ConsensusMode, PoaConfig, consensus, consensus_multi};

/// Result of one [`poa_consensus_run`] call.
///
/// All pointers are owned by the result and released by
/// [`poa_consensus_output_free`].  On error `n_cons` is 0 and `error` holds a
/// NUL-terminated message (or is null if allocation of the message failed).
#[repr(C)]
pub struct PoaConsensusOutput {
    /// Number of consensus sequences (0 on error).
    pub n_cons: usize,
    /// `n_cons` sequence lengths.
    pub cons_lens: *mut usize,
    /// `n_cons` encoded consensus sequences (`0..4`, as in longcallD).
    pub cons_seqs: *mut *mut u8,
    /// `n_cons` cluster sizes.
    pub cluster_n_reads: *mut usize,
    /// `n_cons` arrays of input read indices, one per cluster.
    pub cluster_read_indices: *mut *mut usize,
    /// Error message on failure (NUL-terminated), otherwise null.
    pub error: *mut c_char,
}

impl PoaConsensusOutput {
    fn empty() -> Self {
        Self {
            n_cons: 0,
            cons_lens: ptr::null_mut(),
            cons_seqs: ptr::null_mut(),
            cluster_n_reads: ptr::null_mut(),
            cluster_read_indices: ptr::null_mut(),
            error: ptr::null_mut(),
        }
    }
}

fn decode_read(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        })
        .collect()
}

fn encode_cons(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4,
        })
        .collect()
}

struct InnerResult {
    sequences: Vec<Vec<u8>>,
    clusters: Vec<Vec<usize>>,
}

/// Build `n_cons` consensus sequences and the input-read indices belonging to
/// each one.
///
/// `multi == 0` builds one consensus from all reads.  `multi != 0` runs the
/// multi-allele splitter; at most `max_n_cons` alleles are returned, and any
/// read left outside them is attached to the returned consensus with the
/// closest length (longcallD only models a diploid locus).
unsafe fn run_inner(
    n_reads: usize,
    read_seqs: *const *const u8,
    read_lens: *const usize,
    multi: i32,
    min_reads: usize,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    band_width: usize,
    min_allele_freq: f64,
    max_n_cons: usize,
) -> Result<InnerResult, String> {
    if n_reads == 0 {
        return Err("poa-consensus: no reads provided".to_string());
    }
    if read_seqs.is_null() || read_lens.is_null() {
        return Err("poa-consensus: null reads pointer".to_string());
    }

    // Decode longcallD's 4-bit bases into ASCII for the POA engine.
    let owned: Vec<Vec<u8>> = (0..n_reads)
        .map(|i| {
            let len = *read_lens.add(i);
            if len == 0 {
                Vec::new()
            } else {
                let ptr = *read_seqs.add(i);
                if ptr.is_null() {
                    Vec::new()
                } else {
                    decode_read(slice::from_raw_parts(ptr, len))
                }
            }
        })
        .collect();
    let refs: Vec<&[u8]> = owned.iter().map(|r| r.as_slice()).collect();

    let mut cfg = PoaConfig {
        // abPOA was configured with ABPOA_MF in longcallD, so keep the same
        // consensus extraction rule here.
        consensus_mode: ConsensusMode::MajorityFrequency,
        warn_on_long_unbanded: false,
        ..PoaConfig::default()
    };
    cfg.min_reads = min_reads.max(1);
    cfg.match_score = match_score;
    cfg.mismatch_score = mismatch_score;
    cfg.gap_open = gap_open;
    cfg.gap_extend = gap_extend;
    cfg.band_width = band_width;
    cfg.min_allele_freq = min_allele_freq.clamp(0.0, 1.0);
    if max_n_cons > 0 {
        // no direct field: applied below by truncation
    }

    if multi == 0 {
        let c = consensus(&refs, 0, &cfg).map_err(|e| e.to_string())?;
        return Ok(InnerResult {
            sequences: vec![encode_cons(&c.sequence)],
            clusters: vec![(0..n_reads).collect()],
        });
    }

    let mut alleles = consensus_multi(&refs, 0, &cfg).map_err(|e| e.to_string())?;
    if alleles.is_empty() {
        return Err("poa-consensus: no consensus returned".to_string());
    }
    // Keep the best-supported alleles.  The caller (longcallD) only has room
    // for two haplotypes.
    alleles.sort_by(|a, b| b.read_indices.len().cmp(&a.read_indices.len()));
    alleles.truncate(max_n_cons.max(1));

    let sequences: Vec<Vec<u8>> = alleles.iter().map(|c| encode_cons(&c.sequence)).collect();
    let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); sequences.len()];
    let mut assigned = vec![usize::MAX; n_reads];
    for (ci, allele) in alleles.iter().enumerate() {
        for &ri in &allele.read_indices {
            if ri < n_reads && assigned[ri] == usize::MAX {
                assigned[ri] = ci;
                clusters[ci].push(ri);
            }
        }
    }
    // Reads that the splitter did not place (e.g. dropped third alleles or
    // uninformative reads) still need an alignment target.  Attach each to the
    // allele whose length it is closest to; that is the right call for the
    // length-variable repeats longcallD refines.
    for ri in 0..n_reads {
        if assigned[ri] != usize::MAX {
            continue;
        }
        let mut best = 0usize;
        let mut best_d = usize::MAX;
        for (ci, seq) in sequences.iter().enumerate() {
            let d = owned[ri].len().abs_diff(seq.len());
            if d < best_d {
                best_d = d;
                best = ci;
            }
        }
        clusters[best].push(ri);
    }
    for c in &mut clusters {
        c.sort_unstable();
    }
    Ok(InnerResult {
        sequences,
        clusters,
    })
}

/// Run POA consensus.
///
/// * `read_seqs` — `n_reads` pointers to encoded bases (`0=A, 1=C, 2=G, 3=T,
///   4=N`).
/// * `read_lens` — `n_reads` sequence lengths.
/// * `multi` — 0: single consensus; non-zero: multi-allele splitter.
/// * `match_score` — positive; `mismatch_score`, `gap_open`, `gap_extend` are
///   negative.
///
/// Returns 0 on success and non-zero on error.  `out` is always initialised;
/// release it with [`poa_consensus_output_free`].
#[unsafe(no_mangle)]
pub unsafe extern "C" fn poa_consensus_run(
    n_reads: usize,
    read_seqs: *const *const u8,
    read_lens: *const usize,
    multi: i32,
    min_reads: usize,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    band_width: usize,
    min_allele_freq: f64,
    max_n_cons: usize,
    out: *mut PoaConsensusOutput,
) -> i32 {
    if out.is_null() {
        return -1;
    }
    *out = PoaConsensusOutput::empty();

    let outcome = catch_unwind(AssertUnwindSafe(|| {
        run_inner(
            n_reads,
            read_seqs,
            read_lens,
            multi,
            min_reads,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            band_width,
            min_allele_freq,
            max_n_cons,
        )
    }));

    match outcome {
        Ok(Ok(result)) => {
            let n_cons = result.sequences.len();
            let mut cons_lens: Vec<usize> = Vec::with_capacity(n_cons);
            let mut cons_ptrs: Vec<*mut u8> = Vec::with_capacity(n_cons);
            for seq in result.sequences {
                cons_lens.push(seq.len());
                let boxed = seq.into_boxed_slice();
                cons_ptrs.push(Box::into_raw(boxed) as *mut u8);
            }
            let mut cluster_n: Vec<usize> = Vec::with_capacity(n_cons);
            let mut cluster_ptrs: Vec<*mut usize> = Vec::with_capacity(n_cons);
            for cluster in result.clusters {
                cluster_n.push(cluster.len());
                let boxed = cluster.into_boxed_slice();
                cluster_ptrs.push(Box::into_raw(boxed) as *mut usize);
            }

            (*out).n_cons = n_cons;
            (*out).cons_lens = Box::into_raw(cons_lens.into_boxed_slice()) as *mut usize;
            (*out).cons_seqs = Box::into_raw(cons_ptrs.into_boxed_slice()) as *mut *mut u8;
            (*out).cluster_n_reads = Box::into_raw(cluster_n.into_boxed_slice()) as *mut usize;
            (*out).cluster_read_indices =
                Box::into_raw(cluster_ptrs.into_boxed_slice()) as *mut *mut usize;
            0
        }
        Ok(Err(msg)) => {
            if let Ok(cmsg) = CString::new(msg) {
                (*out).error = cmsg.into_raw();
            }
            -1
        }
        Err(_) => {
            if let Ok(cmsg) = CString::new("poa-consensus: internal panic") {
                (*out).error = cmsg.into_raw();
            }
            -2
        }
    }
}

/// Release all memory owned by a result produced by [`poa_consensus_run`].
/// Safe to call more than once on an already-freed (zeroed) result.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn poa_consensus_output_free(out: *mut PoaConsensusOutput) {
    if out.is_null() {
        return;
    }
    let o = &mut *out;
    let n = o.n_cons;

    if !o.cons_seqs.is_null() {
        let ptrs = Box::from_raw(slice::from_raw_parts_mut(o.cons_seqs, n));
        let lens = if o.cons_lens.is_null() {
            None
        } else {
            Some(Box::from_raw(slice::from_raw_parts_mut(o.cons_lens, n)))
        };
        for (i, p) in ptrs.into_vec().into_iter().enumerate() {
            if !p.is_null() {
                let len = lens.as_ref().map_or(0, |l| l[i]);
                drop(Vec::from_raw_parts(p, len, len));
            }
        }
    } else if !o.cons_lens.is_null() {
        drop(Box::from_raw(slice::from_raw_parts_mut(o.cons_lens, n)));
    }

    if !o.cluster_read_indices.is_null() {
        let ptrs = Box::from_raw(slice::from_raw_parts_mut(o.cluster_read_indices, n));
        let counts = if o.cluster_n_reads.is_null() {
            None
        } else {
            Some(Box::from_raw(slice::from_raw_parts_mut(
                o.cluster_n_reads,
                n,
            )))
        };
        for (i, p) in ptrs.into_vec().into_iter().enumerate() {
            if !p.is_null() {
                let c = counts.as_ref().map_or(0, |l| l[i]);
                drop(Vec::from_raw_parts(p, c, c));
            }
        }
    } else if !o.cluster_n_reads.is_null() {
        drop(Box::from_raw(slice::from_raw_parts_mut(
            o.cluster_n_reads,
            n,
        )));
    }

    if !o.error.is_null() {
        drop(CString::from_raw(o.error));
    }
    *out = PoaConsensusOutput::empty();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(
        reads: &[&[u8]],
        multi: i32,
        max_n_cons: usize,
    ) -> (usize, Vec<Vec<u8>>, Vec<Vec<usize>>) {
        let owned: Vec<Vec<u8>> = reads
            .iter()
            .map(|r| {
                r.iter()
                    .map(|&b| match b {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        b'T' => 3,
                        _ => 4,
                    })
                    .collect()
            })
            .collect();
        let ptrs: Vec<*const u8> = owned.iter().map(|r| r.as_ptr()).collect();
        let lens: Vec<usize> = owned.iter().map(|r| r.len()).collect();
        let mut out = PoaConsensusOutput::empty();
        let ret = unsafe {
            poa_consensus_run(
                owned.len(),
                ptrs.as_ptr(),
                lens.as_ptr(),
                multi,
                1,
                2,
                -6,
                -6,
                -2,
                50,
                0.2,
                max_n_cons,
                &mut out,
            )
        };
        assert_eq!(ret, 0);
        assert!(out.n_cons >= 1);
        let n_cons = out.n_cons;
        let mut seqs = Vec::new();
        let mut clusters = Vec::new();
        unsafe {
            for i in 0..n_cons {
                let len = *out.cons_lens.add(i);
                seqs.push(slice::from_raw_parts(*out.cons_seqs.add(i), len).to_vec());
                let n = *out.cluster_n_reads.add(i);
                clusters.push(slice::from_raw_parts(*out.cluster_read_indices.add(i), n).to_vec());
            }
            poa_consensus_output_free(&mut out);
        }
        (n_cons, seqs, clusters)
    }

    #[test]
    fn single_consensus_roundtrip() {
        let reads: &[&[u8]] = &[
            b"ACGTTGCAACGTTGCAAGG",
            b"ACGTTGCAACGTTGCAAGG",
            b"ACGTTGCAACGTTGCAAGG",
            b"ACGTTGCAACGTTGCAAGG",
        ];
        let (n, seqs, clusters) = run(reads, 0, 1);
        assert_eq!(n, 1);
        let expected: Vec<u8> = b"ACGTTGCAACGTTGCAAGG"
            .iter()
            .map(|&b| match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect();
        assert_eq!(seqs[0], expected);
        assert_eq!(clusters[0], vec![0, 1, 2, 3]);
    }

    #[test]
    fn multi_assigns_every_read() {
        let reads: &[&[u8]] = &[
            b"ACGTTGCAACGTTGCAAGG",
            b"ACGTTGCAACGTTGCAAGG",
            b"ACGTTGCAACGTTGCAAGG",
            b"ACGTTGCAACGTTGCAAGG",
            b"ACGTTGCAAGGTTTGACCA",
            b"ACGTTGCAAGGTTTGACCA",
            b"ACGTTGCAAGGTTTGACCA",
            b"ACGTTGCAAGGTTTGACCA",
        ];
        let (n, _seqs, clusters) = run(reads, 1, 2);
        assert_eq!(n, 2);
        let mut all: Vec<usize> = clusters.iter().flatten().copied().collect();
        all.sort_unstable();
        assert_eq!(all, (0..reads.len()).collect::<Vec<_>>());
    }
}
