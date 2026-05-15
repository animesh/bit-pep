/// amino.rs — 5-bit amino acid encoding for FM-index
///
/// Supports ALL 26 IUPAC single-letter amino acid codes (A–Z).
/// https://en.wikipedia.org/wiki/Amino_acid
///
/// Code assignment:
///   0        → BWT sentinel '$' (reserved; never appears in input)
///   1–20     → Standard 20 proteinogenic AAs
///   21       → U  Selenocysteine  (21st proteinogenic AA)
///   22       → O  Pyrrolysine     (22nd proteinogenic AA)
///   23       → B  Asx  (Asp or Asn — ambiguous)
///   24       → J  Xle  (Leu or Ile — common in MS/MS proteomics; LC-MS can't distinguish)
///   25       → Z  Glx  (Glu or Gln — ambiguous)
///   26       → X  Xaa  (any / unknown)
///
/// 5 bits per residue (range 0–31) → pack 6 residues per u32.
/// ALPHA_SIZE = 27  (codes 0–26).

pub const ALPHA_SIZE: usize = 27;

/// Ordered standard amino-acid alphabet (matches codes 1–20).
pub const AA_ORDER: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// Encode a single IUPAC residue to its 5-bit code.
/// Code 0 is never returned for normal input (reserved for BWT sentinel).
#[inline]
pub fn encode_aa(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        // Standard 20
        b'A' =>  1,  // Alanine
        b'C' =>  2,  // Cysteine
        b'D' =>  3,  // Aspartic acid
        b'E' =>  4,  // Glutamic acid
        b'F' =>  5,  // Phenylalanine
        b'G' =>  6,  // Glycine
        b'H' =>  7,  // Histidine
        b'I' =>  8,  // Isoleucine
        b'K' =>  9,  // Lysine
        b'L' => 10,  // Leucine
        b'M' => 11,  // Methionine
        b'N' => 12,  // Asparagine
        b'P' => 13,  // Proline
        b'Q' => 14,  // Glutamine
        b'R' => 15,  // Arginine
        b'S' => 16,  // Serine
        b'T' => 17,  // Threonine
        b'V' => 18,  // Valine
        b'W' => 19,  // Tryptophan
        b'Y' => 20,  // Tyrosine
        // Special proteinogenic
        b'U' => 21,  // Selenocysteine (Sec)
        b'O' => 22,  // Pyrrolysine    (Pyl)
        // IUPAC ambiguity codes
        b'B' => 23,  // Asx  = Asp or Asn
        b'J' => 24,  // Xle  = Leu or Ile  (important for MS/MS)
        b'Z' => 25,  // Glx  = Glu or Gln
        b'X' => 26,  // Xaa  = any amino acid
        _    => 26,  // treat anything else as X
    }
}

/// Decode a code back to its single-letter IUPAC symbol.
#[inline]
pub fn decode_aa(code: u8) -> u8 {
    match code {
        1..=20 => AA_ORDER[(code - 1) as usize],
        21 => b'U',
        22 => b'O',
        23 => b'B',
        24 => b'J',
        25 => b'Z',
        _  => b'X',
    }
}

/// Encode an entire peptide/protein sequence into a Vec<u8> of AA codes.
/// Stops at '*' (stop codon) if present.
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
    fn encode_decode_roundtrip_standard_20() {
        for &aa in AA_ORDER.iter() {
            let code = encode_aa(aa);
            assert!(code >= 1 && code <= 20, "code out of range for {}", aa as char);
            assert_eq!(decode_aa(code), aa);
        }
    }

    #[test]
    fn all_26_iupac_codes_are_nonzero() {
        for c in b'A'..=b'Z' {
            let code = encode_aa(c);
            assert!(code >= 1, "code 0 (sentinel) leaked for '{}'", c as char);
            assert!(code <= 26, "code > 26 for '{}'", c as char);
        }
    }

    #[test]
    fn j_encodes_as_xle() {
        let code = encode_aa(b'J');
        assert_eq!(code, 24, "J (Xle) should be code 24");
        assert_eq!(decode_aa(24), b'J');
    }

    #[test]
    fn special_codes_roundtrip() {
        assert_eq!(encode_aa(b'U'), 21); assert_eq!(decode_aa(21), b'U');
        assert_eq!(encode_aa(b'O'), 22); assert_eq!(decode_aa(22), b'O');
        assert_eq!(encode_aa(b'B'), 23); assert_eq!(decode_aa(23), b'B');
        assert_eq!(encode_aa(b'J'), 24); assert_eq!(decode_aa(24), b'J');
        assert_eq!(encode_aa(b'Z'), 25); assert_eq!(decode_aa(25), b'Z');
        assert_eq!(encode_aa(b'X'), 26); assert_eq!(decode_aa(26), b'X');
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
        let a: Vec<u8>   = b"ACDEFG".iter().map(|&c| encode_aa(c)).collect();
        let b_seq: Vec<u8> = b"ACDEFH".iter().map(|&c| encode_aa(c)).collect();
        assert_eq!(count_mismatches_u32(pack_u32(&a), pack_u32(&b_seq)), 1);
    }
}
