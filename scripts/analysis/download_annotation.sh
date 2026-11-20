#!/bin/bash
# Public annotation files used by scripts/analysis/selection_overlap_annot.py (Result 4), saved in annotation/.
# The published results used the versions downloaded on 18 November 2026. GO, HGNC, Reactome and the GWAS Catalog
# are updated over time, so a later download can change the enrichment results slightly; the known-AD-gene list used
# is in results/selection_overlap/annot/gwas_ad_genes.csv.
set -e
cd "$(dirname "$0")/../.."
mkdir -p annotation/gwas_ad && cd annotation
get() { [ -s "$2" ] || curl -fsSL --retry 5 -o "$2" "$1"; }
get https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.basic.annotation.gtf.gz gencode.v44lift37.basic.annotation.gtf.gz
get https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt hgnc_complete_set.txt
get http://current.geneontology.org/annotations/goa_human.gaf.gz goa_human.gaf.gz
get http://purl.obolibrary.org/obo/go/go-basic.obo go-basic.obo
get https://reactome.org/download/current/ReactomePathways.gmt.zip ReactomePathways.gmt.zip
# GWAS Catalog REST API v2: all associations with Alzheimer disease (MONDO_0004975), 100 per page
p=0
while :; do
  f=gwas_ad/page_$p.json
  get "https://www.ebi.ac.uk/gwas/rest/api/v2/associations?efo_id=MONDO_0004975&size=100&page=$p" "$f"
  python3 -c "import json,sys; d=json.load(open('$f')); pg=d.get('page',{}); sys.exit(0 if pg.get('number',0)+1 < pg.get('totalPages',0) else 1)" || break
  p=$((p+1))
done
echo "annotation files in $(pwd)"
