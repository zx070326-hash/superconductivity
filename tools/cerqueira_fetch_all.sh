#!/bin/bash
# Cerqueira 6912 (Materials Cloud record 3kbt5-r3n56, CC-BY-4.0): fetch the 16 remaining
# batches, verify byte size AND md5 against the record metadata (api/records JSON,
# read 2026-09-21), extract only the per-material summary files
# (McMillan.dat / Eliashberg.dat / Tc.dat / a2F.dos*), keep the tarballs.
# Run on petlove-train:  nohup /zx/cerqueira6912/fetch_all.sh > /zx/cerqueira6912/log/nohup.out 2>&1 &
cd /zx/cerqueira6912 || exit 1
mkdir -p tar extracted log
declare -A SZ=( [a]=3227483011 [b]=3438986099 [c]=3424761873 [d]=3494234906 [e]=3378142502
  [f]=3242369863 [g]=3121983436 [h]=3221363829 [i]=3028598664 [j]=3097509131 [k]=3026215208
  [l]=3023483754 [m]=3025352999 [n]=3120868222 [o]=2992368798 [p]=2623380238 )
declare -A MD=( [a]=b7f75b330e0b6e0573d15f172c322707 [b]=63ee110453a8d2588c5271d71134e5a3
  [c]=77a3927603fa4654f003742b9f37e740 [d]=00d933a8bf6041997c0242ff7dd53f5a
  [e]=5e4724b49246796a9b8125f8ac0b3467 [f]=96801dd4489b6c3351edee326b58f7bb
  [g]=ff32a36bd5e7917b0167a803db073756 [h]=c9739508448664a4628116acc2c86912
  [i]=0d10ed458c403bd9c31705cc823135d3 [j]=b1ab88ba69dd299b99550700e628d46d
  [k]=d4f879b6d33864777b08380c61cdb8ab [l]=dda7c023c250a343f7de7753de6ed67c
  [m]=0bbecf3bf187552fbd4a76e92f926ce3 [n]=94adedd49761a0a02d603f02f6b64869
  [o]=7a5fe557d8f0aae014a5fad7e1050359 [p]=a2b8ca9b8116879ac41a984f738a1980 )
# smallest first so usable data appears early
for b in p o n m l k j i h g f e d c b a; do
  f=tar/batch-$b.tar.bz2
  url="https://archive.materialscloud.org/api/records/3kbt5-r3n56/files/batch-$b.tar.bz2/content"
  # a previous wrong-URL run may have left a 404 HTML page here
  if [ -f "$f" ] && [ "$(stat -c %s "$f")" -lt 1000000 ]; then rm -f "$f"; fi
  echo "[$(date -Is)] START batch-$b" >> log/progress.log
  for attempt in 1 2 3 4 5 6 7 8; do
    curl -sL -C - --retry 20 --retry-delay 2 --retry-all-errors --speed-limit 1 --speed-time 30 \
         -o "$f" "$url" 2>>log/curl-$b.log
    have=$(stat -c %s "$f" 2>/dev/null || echo 0)
    [ "$have" -eq "${SZ[$b]}" ] && break
    echo "[$(date -Is)] batch-$b attempt $attempt size $have / ${SZ[$b]}, resuming" >> log/progress.log
    sleep 10
  done
  have=$(stat -c %s "$f" 2>/dev/null || echo 0)
  if [ "$have" -ne "${SZ[$b]}" ]; then
    echo "[$(date -Is)] FAIL-SIZE batch-$b size $have" >> log/progress.log
    continue
  fi
  sum=$(md5sum "$f" | cut -d" " -f1)
  if [ "$sum" != "${MD[$b]}" ]; then
    echo "[$(date -Is)] FAIL-MD5 batch-$b $sum != ${MD[$b]}" >> log/progress.log
    continue
  fi
  echo "[$(date -Is)] DONE-DL batch-$b md5 OK" >> log/progress.log
  mkdir -p extracted/batch-$b
  tar -xjf "$f" -C extracted/batch-$b --wildcards --no-anchored \
      "McMillan.dat" "Eliashberg.dat" "Tc.dat" "a2F.dos*" 2>>log/tar-$b.log
  n=$(find extracted/batch-$b -name McMillan.dat | wc -l)
  echo "[$(date -Is)] DONE-EXTRACT batch-$b materials=$n" >> log/progress.log
done
echo "[$(date -Is)] ALL FINISHED" >> log/progress.log
