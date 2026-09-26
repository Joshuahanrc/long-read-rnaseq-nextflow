#!/usr/bin/env bash

set -euo pipefail

qc_dir="${1:-task5_results/qc}"
output_file="${2:-qc_summary.tsv}"

if [[ ! -d "$qc_dir" ]]; then
    echo "Error: QC directory not found: $qc_dir" >&2
    exit 1
fi

printf "sample\traw_sequences\tprimary_mapped\tmapped_pct\tsecondary\tsupplementary\tavg_length\tmax_length\tavg_quality\terror_rate\n" \
    > "$output_file"

shopt -s nullglob
stats_files=("$qc_dir"/*.stats.txt)

if [[ ${#stats_files[@]} -eq 0 ]]; then
    echo "Error: no .stats.txt files found in $qc_dir" >&2
    exit 1
fi

for stats_file in "${stats_files[@]}"; do
    sample="$(basename "$stats_file" .stats.txt)"
    flagstat_file="$qc_dir/${sample}.flagstat.txt"

    if [[ ! -f "$flagstat_file" ]]; then
        echo "Warning: flagstat file missing for $sample" >&2
        continue
    fi

    raw_sequences=$(awk -F '\t' \
        '$1 == "SN" && $2 == "raw total sequences:" {print $3; exit}' \
        "$stats_file")

    avg_length=$(awk -F '\t' \
        '$1 == "SN" && $2 == "average length:" {print $3; exit}' \
        "$stats_file")

    max_length=$(awk -F '\t' \
        '$1 == "SN" && $2 == "maximum length:" {print $3; exit}' \
        "$stats_file")

    avg_quality=$(awk -F '\t' \
        '$1 == "SN" && $2 == "average quality:" {print $3; exit}' \
        "$stats_file")

    error_rate=$(awk -F '\t' \
        '$1 == "SN" && $2 == "error rate:" {
            printf "%.2f%%", $3 * 100
            exit
        }' "$stats_file")

    primary_mapped=$(awk \
        '$4 == "primary" && $5 == "mapped" {print $1; exit}' \
        "$flagstat_file")

    mapped_pct=$(awk \
        '$4 == "primary" && $5 == "mapped" {
            gsub(/[()]/, "", $6)
            print $6
            exit
        }' "$flagstat_file")

    secondary=$(awk \
        '$4 == "secondary" {print $1; exit}' \
        "$flagstat_file")

    supplementary=$(awk \
        '$4 == "supplementary" {print $1; exit}' \
        "$flagstat_file")

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$raw_sequences" \
        "$primary_mapped" \
        "$mapped_pct" \
        "$secondary" \
        "$supplementary" \
        "$avg_length" \
        "$max_length" \
        "$avg_quality" \
        "$error_rate" \
        >> "$output_file"
done

echo "QC summary written to: $output_file"
