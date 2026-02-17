#!/usr/bin/env bash
set -euo pipefail

base_dir="${1:-}"
if [ -z "$base_dir" ]; then
    echo "Usage: $0 <directory>" >&2
    exit 1
fi

outputs_dir="$base_dir/outputs"
if [ ! -d "$outputs_dir" ]; then
    echo "Missing outputs directory: $outputs_dir" >&2
    exit 1
fi

# Collect files safely using null-delimited output into an array (portable)
all_files=()
while IFS= read -r -d '' f; do
    all_files+=("$f")
done < <(find "$outputs_dir" -type f -name '*.root' -print0)

numFiles=${#all_files[@]}
echo "Merging $numFiles files from $outputs_dir into ${base_dir}/merged_output.root"

if (( numFiles == 0 )); then
    echo "No .root files found in $outputs_dir. Nothing to merge."
    exit 0
fi

if (( numFiles > 100 )); then
    tmp_merge_dir="$(mktemp -d -p "$base_dir" submerge.XXXXXX)"
    # Ensure cleanup happens even if the script exits early
    cleanup() {
        rm -rf -- "$tmp_merge_dir" 2>/dev/null || true
    }
    trap cleanup EXIT

    tmp_merged_files=()

    # Merge each immediate subdirectory into a temporary merged file
    while IFS= read -r -d '' dir; do
        subfiles=()
        while IFS= read -r -d '' file; do
            subfiles+=("$file")
        done < <(find "$dir" -type f -name '*.root' -print0)
        if (( ${#subfiles[@]} == 0 )); then
            continue
        fi
        sub_out="$tmp_merge_dir/$(basename "$dir")_merged.root"
        hadd -f "$sub_out" "${subfiles[@]}"
        tmp_merged_files+=("$sub_out")
    done < <(find "$outputs_dir" -mindepth 1 -maxdepth 1 -type d -print0)

    # Also merge any root files directly under outputs_dir
    topfiles=()
    while IFS= read -r -d '' file; do
        topfiles+=("$file")
    done < <(find "$outputs_dir" -maxdepth 1 -type f -name '*.root' -print0)

    if (( ${#topfiles[@]} > 0 )); then
        top_out="$tmp_merge_dir/_rootlevel_merged.root"
        hadd -f "$top_out" "${topfiles[@]}"
        tmp_merged_files+=("$top_out")
    fi

    # Final merge of all temporary merged files
    hadd -f "${base_dir}/merged_output.root" "${tmp_merged_files[@]}"
else
    hadd -f "${base_dir}/merged_output.root" "${all_files[@]}"
fi
