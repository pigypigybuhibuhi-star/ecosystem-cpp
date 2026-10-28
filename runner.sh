#!/bin/bash
# 週末の待ち行列を一つずつ実行する。計画は Claude が runs/queue.txt に書く。この台本は Claude を使わない
export PATH=/opt/homebrew/bin:/usr/local/bin:$PATH
R=~/research/evosim; Q=$R/runs/queue.txt; DONE=$R/runs/queue_done.txt
pgrep -f "caffeinate -i" >/dev/null && exit 0          # 何か動いていれば待つ(夜の部の実行も含む)
[ -f "$Q" ] || exit 0
line=$(grep -v -e '^[[:space:]]*#' -e '^[[:space:]]*$' "$Q" | head -1)
[ -z "$line" ] && exit 0
IFS='|' read -r name commit hours args <<< "$line"
name=$(echo "$name" | xargs); commit=$(echo "$commit" | xargs); hours=$(echo "$hours" | xargs)
grep -v -F -x "$line" "$Q" > "$Q.tmp"; mv "$Q.tmp" "$Q"
echo "$(date '+%F %T') start | $line" >> "$DONE"
D=$R/runs/queue/$name; mkdir -p "$D"; cd "$D" || exit 1
git -C "$R" fetch -q
git -C "$R" archive "$commit" | tar -x -C "$D" || { echo "$(date '+%F %T') fail(source) | $name" >> "$DONE"; exit 1; }
git -C "$R" rev-parse "$commit" > commit.txt
clang++ -std=c++17 -O2 phase3.cpp -o phase3 $(pkg-config --cflags --libs sfml-graphics) -framework IOKit -framework CoreFoundation > build.log 2>&1 \
  || { echo "$(date '+%F %T') fail(build) | $name" >> "$DONE"; exit 1; }
read -ra A <<< "$args"
nohup sh -c 'caffeinate -i ./phase3 "$@" > run.log 2>&1; echo "exit=$? $(date)" > done.txt' sh "${A[@]}" < /dev/null > /dev/null 2>&1 &
