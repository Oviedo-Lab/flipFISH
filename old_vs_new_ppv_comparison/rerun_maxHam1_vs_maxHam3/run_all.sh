#!/bin/bash
set -e
cd /Users/michaelbarkasi/projects_Oviedo_lab/R_packages/flipFISH

MICE=(m1 m2 m3 m4 m8)
OUTDIR=old_vs_new_ppv_comparison/rerun_maxHam1_vs_maxHam3/qc_outputs

for maxHam in 3 1; do
  for mouse in "${MICE[@]}"; do
    out="${OUTDIR}/qc_real_${mouse}_maxHam${maxHam}.rds"
    if [ -f "$out" ]; then
      echo "=== SKIP (exists): mouse=$mouse maxHam=$maxHam ==="
      continue
    fi
    echo "=== START: mouse=$mouse maxHam=$maxHam $(date) ==="
    Rscript old_vs_new_ppv_comparison/rerun_maxHam1_vs_maxHam3/run_one.R "$mouse" "$maxHam" "$out"
    echo "=== DONE: mouse=$mouse maxHam=$maxHam $(date) ==="
  done
done

echo "ALL RUNS COMPLETE"
