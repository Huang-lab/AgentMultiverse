#!/usr/bin/env bash
# Pack a finished study into the blinded archive as <study>.zip, verify it, then delete the unpacked copy.
# Usage (from the project root): tools/archive_study.sh STUDY
#   STUDY is a plain folder name, found in results/ or, if already moved, in the archive.
# This is the sanctioned way to put studies into the archive: it prints only sizes and file counts,
# never archive contents, so agents in this project stay blind to earlier results.
set -euo pipefail
ARCHIVE="${AMX_ARCHIVE:-$HOME/Projects/AgentMultiverse_archive}"   # override only for testing
study=${1:?usage: tools/archive_study.sh STUDY}
[[ "$study" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "STUDY must be a plain folder name" >&2; exit 1; }
mkdir -p "$ARCHIVE"

if [ -d "results/$study" ]; then parent="$PWD/results"
elif [ -d "$ARCHIVE/$study" ]; then parent="$ARCHIVE"
else echo "no folder named $study in results/ or in the archive" >&2; exit 1; fi
zipf="$ARCHIVE/$study.zip"
[ -e "$zipf" ] && { echo "$study.zip already exists in the archive" >&2; exit 1; }

before=$(du -sh "$parent/$study" | cut -f1)

# Byte-identical copies of shared-prep files (listed in prep/*/MD5SUMS) are replaced by a small
# <file>.prepref pointer holding the prep path and md5, so archives do not store LD matrices again.
n_ref=0
if compgen -G "prep/*/MD5SUMS" > /dev/null; then
  while IFS= read -r -d '' f; do
    h=$(md5 -q "$f" 2>/dev/null || md5sum "$f" | cut -d' ' -f1)
    hit=$(awk -v h="$h" '$1 == h { d = FILENAME; sub(/MD5SUMS$/, "", d); print d $2; exit }' prep/*/MD5SUMS)
    if [ -n "$hit" ]; then
      printf '%s\t%s\n' "$hit" "$h" > "$f.prepref" && rm -f "$f"
      n_ref=$((n_ref + 1))
    fi
  done < <(find "$parent/$study" -type f -size +10M -print0)
fi
[ "$n_ref" -gt 0 ] && echo "replaced $n_ref copies of shared-prep files with .prepref pointers"

n_files=$(find "$parent/$study" \( -type f -o -type l \) | wc -l | tr -d ' ')
rm -f "$zipf.part"
(cd "$parent" && zip -r -q -y "$zipf.part" "$study")
unzip -tqq "$zipf.part" > /dev/null || { echo "zip test failed; unpacked copy kept" >&2; exit 1; }
n_zip=$(zipinfo -1 "$zipf.part" | grep -vc '/$')
[ "$n_zip" = "$n_files" ] || { echo "file count mismatch ($n_zip in zip, $n_files on disk); unpacked copy kept" >&2; exit 1; }
mv "$zipf.part" "$zipf"
chmod -R u+w "$parent/$study" && rm -rf "${parent:?}/$study"
echo "archived $study: $n_files files, $before unpacked -> $(du -h "$zipf" | cut -f1) zip (verified)"
