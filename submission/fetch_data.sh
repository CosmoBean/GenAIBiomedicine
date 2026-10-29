#!/usr/bin/env bash
# Download every input file used in the submission into ../data/ and verify checksums.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
D="$ROOT/data"; mkdir -p "$D/visium" "$D/protein"; cd "$D"
B=https://github.com/GenAIBioMed/GenAIBioMedAssignment3/raw/refs/heads/main/datasets
get() { [ -s "$2" ] || curl -sSL -A "Mozilla/5.0" -o "$2" "$1"; }

# Task 2 / Task 3 (course repository copies of 10x Genomics data)
get $B/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz pbmc3k_filtered_gene_bc_matrices.tar.gz
get $B/visium_mouse_brain/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5 visium/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5
get $B/visium_mouse_brain/V1_Adult_Mouse_Brain_spatial.tar.gz V1_Adult_Mouse_Brain_spatial.tar.gz
sha256sum -c <<'EOF'
847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5  pbmc3k_filtered_gene_bc_matrices.tar.gz
eb78379e02dcf48036abf05b67233e73ecb0d880787feb82f76ff16f6ce01eb3  visium/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5
46d6b05ba740f232d6bf4b27b9a8846815851e000985fb878f1364bab04e5bd4  V1_Adult_Mouse_Brain_spatial.tar.gz
EOF
[ -d filtered_gene_bc_matrices ] || tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
[ -d visium/spatial ] || tar -xzf V1_Adult_Mouse_Brain_spatial.tar.gz -C visium

# Task 1
cd protein
get $B/protein/4hjo.pdb 4hjo.pdb
echo "95b122b5e152e705fcbcb4314be07ce9378c5a64ab3226e06046ca24ca674948  4hjo.pdb" | sha256sum -c
for id in 1M17 2ITX 2ITN 2ITV 2JIT 4I22; do get https://files.rcsb.org/download/$id.pdb "$(echo $id | tr A-Z a-z).pdb"; done
get https://rest.uniprot.org/uniprotkb/P00533.json P00533.json
get https://rest.uniprot.org/uniprotkb/P00533.fasta P00533.fasta
for p in P04626 P21860 Q15303; do curl -sSL https://rest.uniprot.org/uniprotkb/$p.fasta; done > erbb_paralogs.fasta
AF=https://alphafold.ebi.ac.uk/files
get $AF/AF-P00533-F1-model_v6.pdb AF-P00533-F1-model_v6.pdb
get $AF/AF-P00533-F1-confidence_v6.json AF-P00533-F1-confidence_v6.json
get $AF/AF-P00533-F1-aa-substitutions.csv AF-P00533-F1-aa-substitutions.csv
echo "data ready in $D"
