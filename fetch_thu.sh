#!/bin/bash
# 并行分片下载 THU-ATOM/dekois_airdd.zip（hf-mirror 单连接被限速 ~160KB/s，分片绕过）
set -e
cd "$(dirname "$0")"
URL="https://hf-mirror.com/datasets/THU-ATOM/dekois_airdd/resolve/main/dekois_airdd.zip"
SIZE=711510511
N=12
CHUNK=$(( (SIZE + N - 1) / N ))
mkdir -p data/raw/thu_parts
echo "start $(date)"
for i in $(seq 0 $((N-1))); do
  start=$((i*CHUNK)); end=$(( (i+1)*CHUNK - 1 ))
  [ $end -ge $SIZE ] && end=$((SIZE-1))
  ( for try in 1 2 3; do
      curl -sL -r $start-$end -o data/raw/thu_parts/part_$(printf %02d $i) "$URL" && break
      sleep 3
    done ) &
done
wait
echo "parts done $(date)"
cat data/raw/thu_parts/part_* > data/raw/dekois_airdd.zip
ls -la data/raw/dekois_airdd.zip
unzip -t data/raw/dekois_airdd.zip > /tmp/thu_ziptest.log 2>&1 && echo "ZIP_OK" || echo "ZIP_BROKEN"
echo "end $(date)"
