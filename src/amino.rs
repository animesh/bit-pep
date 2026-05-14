/// amino.rs — 5-bit amino acid encoding for FM-index
///
/// Replaces the 2-bit DNA encoding in lib.rs.
/// Standard 20 AAs → codes 1–20 (0 reserved for sentinel '$').
/// Unknown / ambiguous residues (B, X, Z, U, O) → code 21.
/// Code 0 is the BWT sentinel; libsais requires it to be the lexicographically
/// smallest character, so all AA codes must be ≥ 1.
///
/// ALPHA_SIZE = 22 (sentinel + 20 AAs + unknown bucket).

pub const ALPHA_SIZE: usize = 22;

/// Ordered standard amino-acid alphabet (matches code 1–20).
pub const AA_ORDER: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// Encode a single residue to its 5-bit code.
/// Returns 0 only for the BWT sentinel '$'; never returned for normal input.
#[inline]
pub fn encode_aa(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' =>  1,
        b'C' =>  2,
        b'D' =>  3,
        b'E' =>  4,
        b'F' =>  5,
        b'G' =>  6,
        b'H' =>  7,
        b'I' =>  8,
        b'K' =>  9,
        b'L' => 10,
        b'M' => 11,
        b'N' => 12,
        b'P' => 13,
        b'Q' => 14,
        b'R' => 15,
        b'S' => 16,
        b'T' => 17,
        b'V' => 18,
        b'W' => 19,
        b'Y' => 20,
        // Ambiguous / non-standard: B(D|N), Z(E|Q), X(any), U(Sec), O(Pyl)
        _ =>    21,
    }
}

/// Decode a code back to single-letter AA (for display).
#[inline]
pub fn decode_aa(code: u8) -> u8 {
    match code {
        1..=20 => AA_ORDER[(code - 1) as usize],
        _ => b'X',
    }
}

/// Encode an entire peptide/protein sequence into a Vec<u8> of AA codes.
/// Stops at '*' (stop codon) if present — treats it as end of sequence.
pub fn encode_sequence(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &c in seq {
        if c == b'*' { break; }
        out.push(encode_aa(c));
    }
    out
}

/// Pack 6 encoded AAs into one u32 (5 bits × 6 = 30 bits; 2 bits unused).
/// Used for XOR-based approximate matching in align.rs.
/// `encoded` must already be encode_sequence() output.
#[inline]
pub fn pack_u32(encoded: &[u8]) -> u32 {
    debug_assert!(encoded.len() <= 6);
    let mut word: u32 = 0;
    for (i, &code) in encoded.iter().enumerate() {
        word |= (code as u32) << (i * 5);
    }
    word
}

/// Count mismatching AA positions in two packed u32 words.
/// Each 5-bit field is non-zero iff the residues differ.
#[inline]
pub fn count_mismatches_u32(a: u32, b: u32) -> u32 {
    const MASK5: u32 = 0b11111;
    let xor = a ^ b;
    let mut mismatches = 0u32;
    for i in 0..6u32 {
        if (xor >> (i * 5)) & MASK5 != 0 {
            mismatches += 1;
        }
    }
    mismatches
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_decode_roundtrip() {
        for &aa in AA_ORDER.iter() {
            let code = encode_aa(aa);
            assert!(code >= 1 && code <= 20, "code out of range for {}", aa as char);
            assert_eq!(decode_aa(code), aa);
        }
    }

    #[test]
    fn ambiguous_residues_get_code_21() {
        assert_eq!(encode_aa(b'X'), 21);
        assert_eq!(encode_aa(b'B'), 21);
        assert_eq!(encode_aa(b'Z'), 21);
    }

    #[test]
    fn pack_xor_identical_is_zero_mismatches() {
        let seq = b"ACDEFG";
        let enc: Vec<u8> = seq.iter().map(|&c| encode_aa(c)).collect();
        let w = pack_u32(&enc);
        assert_eq!(count_mismatches_u32(w, w), 0);
    }

    #[test]
    fn pack_xor_counts_one_mismatch() {
        let a: Vec<u8> = b"ACDEFG".iter().map(|&c| encode_aa(c)).collect();
        let b_seq: Vec<u8> = b"ACDEFH".iter().map(|&c| encode_aa(c)).collect(); // last AA differs
        assert_eq!(count_mismatches_u32(pack_u32(&a), pack_u32(&b_seq)), 1);
    }
}
