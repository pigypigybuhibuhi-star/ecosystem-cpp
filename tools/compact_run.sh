#!/bin/bash
# 終わった実行のフォルダを小さくする(作者が 2026-11-17 に追加)
#  - 一番新しい autosave*.txt だけを final_save.txt.gz に圧縮して残し、他の autosave は消す
#    (元の名前は final_save.from に書く。再開や試験器に使うときは gunzip -k final_save.txt.gz で戻す)
#  - 1MB を超える brains_*.txt と *.log を gzip で圧縮する
#  - 30分以内に書かれたセーブがあり、phase3 が動いていれば、まだ実行中かもしれないので触らない
# 使い方: tools/compact_run.sh [--dry] フォルダ...   (--dry は消さずに、何をするかだけを表示)
dry=0; [ "$1" = "--dry" ] && { dry=1; shift; }
for D in "$@"; do
  [ -d "$D" ] || { echo "無い: $D"; continue; }
  if [ -n "$(find "$D" -maxdepth 1 -name 'autosave*.txt' -mmin -30)" ] && pgrep -x phase3 >/dev/null; then
    echo "触らない(実行中かもしれない): $D"; continue
  fi
  before=$(du -sk "$D" | cut -f1)
  newest=$(ls -t "$D"/autosave*.txt 2>/dev/null | head -1)
  big=()
  for f in "$D"/brains_*.txt "$D"/*.log; do [ -f "$f" ] && [ "$(wc -c < "$f")" -gt 1048576 ] && big+=("$f"); done
  if [ $dry = 1 ]; then
    old=$(ls "$D"/autosave*.txt 2>/dev/null | grep -v -x -F "$newest" | xargs -n1 basename 2>/dev/null | tr '\n' ' ')
    echo "[見積もり] $D ($((before/1024))MB): 残して圧縮=$(basename "${newest:-無し}")  消す=${old:-無し}  圧縮=${#big[@]}個"
    continue
  fi
  if [ -n "$newest" ]; then
    gzip -6 -c "$newest" > "$D/final_save.txt.gz.part" && mv "$D/final_save.txt.gz.part" "$D/final_save.txt.gz" \
      && basename "$newest" > "$D/final_save.from" && rm -f "$D"/autosave*.txt
  fi
  for f in "${big[@]}"; do gzip -6 "$f"; done
  after=$(du -sk "$D" | cut -f1)
  echo "$D: $((before/1024))MB -> $((after/1024))MB"
done
