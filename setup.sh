#!/usr/bin/env bash
# setup_prot.sh (corrected) — patches Cargo.toml + lib.rs, compiles, tests
# Run from root of the bit-pop repo:
#   chmod +x setup_prot_v2.sh && ./setup_prot_v2.sh

set -euo pipefail
RESET='\033[0m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'
ok()   { echo -e "${GREEN}✓ $*${RESET}"; }
warn() { echo -e "${YELLOW}⚠ $*${RESET}"; }
die()  { echo -e "${RED}✗ $*${RESET}"; exit 1; }

[[ -f Cargo.toml ]] || die "Run this from the bit-pop repo root"
[[ -f src/lib.rs ]] || die "src/lib.rs not found"
for f in amino peptide_input uniprot pep_output; do
    [[ -f "src/${f}.rs" ]] || die "src/${f}.rs not found"
done
ok "New source files present"

echo ""
echo "=== Cargo.toml ==="

add_dep() {
    local name="$1" line="$2"
    if grep -qE "^${name}\s*=" Cargo.toml; then
        ok "${name} already present"
    else
        sed -i "/^\[dependencies\]/a ${line}" Cargo.toml
        ok "Added: ${line}"
    fi
}

add_dep "reqwest"     'reqwest     = { version = "0.12", features = ["json", "stream"] }'
add_dep "urlencoding" 'urlencoding  = "2.1"'
add_dep "tokio"       'tokio       = { version = "1", features = ["full"] }'
add_dep "serde"       'serde       = { version = "1", features = ["derive"] }'
add_dep "serde_json"  'serde_json  = "1"'
add_dep "anyhow"      'anyhow      = "1"'

echo ""
echo "=== lib.rs module block ==="

# Insert after "pub mod serialize;" in REVERSE order so final order is
# amino / peptide_input / uniprot / pep_output
for mod in pep_output uniprot peptide_input amino; do
    if grep -q "^pub mod ${mod};" src/lib.rs; then
        ok "pub mod ${mod} already declared"
    else
        sed -i '/^pub mod serialize;/a pub mod '"${mod}"';' src/lib.rs
        ok "Added: pub mod ${mod};"
    fi
done

echo ""
echo "Module block now:"
grep "^pub mod" src/lib.rs

echo ""
echo "=== cargo check ==="
cargo check 2>&1 | tail -20

echo ""
echo "=== Unit tests (no network) ==="
for mod in amino peptide_input pep_output; do
    echo "--- ${mod} ---"
    cargo test --lib ${mod} -- --nocapture 2>&1 | grep -E "(test |running|FAILED|error\[)" || true
done

ok "Done. For network tests: cargo test --test uniprot_integration -- --nocapture"
