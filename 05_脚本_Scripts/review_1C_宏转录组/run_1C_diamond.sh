#!/bin/bash
set -e
cd ./metatx
DIAMOND=./downloads/diamond
# build pooled metatranscriptome protein DB from downloaded samples
echo "building pooled DB from $(ls cds/*.faa.gz | wc -l) samples..."
zcat cds/*.faa.gz > pooled_metatx.faa
echo "pooled ORFs: $(grep -c '>' pooled_metatx.faa)"
$DIAMOND makedb --in pooled_metatx.faa -d metatx_db --threads 24 2>/dev/null
# blastp candidates vs metatranscriptome ORFs
$DIAMOND blastp -q query_passers.faa -d metatx_db -o hits.tsv \
  --very-sensitive --max-target-seqs 25 --evalue 1e-5 --id 60 \
  --query-cover 0 --threads 24 \
  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore 2>/dev/null
echo "hits: $(wc -l < hits.tsv)"
rm -f pooled_metatx.faa  # free space
