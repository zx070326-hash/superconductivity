#!/bin/bash
# Supplemental extraction: the per-material `lambda` file (README: "Eliashberg lambda
# calculated for different smearings") was NOT in fetch_all.sh's tar filter.
# Runs over whatever tarballs are already md5-verified in tar/; safe to re-run.
# Does NOT touch fetch_all.sh, which is still running (bash reads scripts incrementally,
# so editing the running file would corrupt it).
#   nohup /zx/cerqueira6912/extract_lambda.sh > /zx/cerqueira6912/log/lambda.out 2>&1 &
cd /zx/cerqueira6912 || exit 1
for f in tar/batch-*.tar.bz2; do
  b=$(basename "$f" .tar.bz2)          # batch-x
  # only batches fetch_all.sh already finished (marker: the extracted dir has McMillan.dat)
  n=$(find "extracted/$b" -name McMillan.dat 2>/dev/null | wc -l)
  [ "$n" -eq 0 ] && { echo "[$(date -Is)] SKIP $b (not yet extracted)"; continue; }
  have=$(find "extracted/$b" -name lambda 2>/dev/null | wc -l)
  if [ "$have" -eq "$n" ]; then echo "[$(date -Is)] SKIP $b (lambda already $have/$n)"; continue; fi
  echo "[$(date -Is)] START $b  (McMillan=$n, lambda=$have)"
  tar -xjf "$f" -C "extracted/$b" --wildcards --no-anchored "lambda" 2>>log/tarlam-$b.log
  after=$(find "extracted/$b" -name lambda 2>/dev/null | wc -l)
  echo "[$(date -Is)] DONE  $b  lambda=$after/$n"
done
echo "[$(date -Is)] LAMBDA EXTRACTION FINISHED"
