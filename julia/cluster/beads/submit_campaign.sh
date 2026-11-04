#!/bin/bash
# Master script of the bead campaign: submits every stage NOW, each one waiting for the
# previous one through Slurm dependencies. Nothing runs on the login node except the sbatch
# calls; after this script ends you can log out. Slurm starts each stage when its
# predecessor ends.
#
#   cd /oscar/data/dharri15/eaguerov/Github/km-disk-vibrating
#   bash julia/cluster/beads/submit_campaign.sh --dry-run     # print what would be submitted
#   bash julia/cluster/beads/submit_campaign.sh               # submit
#
# Chain (A -> B means B starts only if A ended with exit code 0; "any" = B starts when A has
# ended in any way, so a failure of A does not block B):
#
#   00a self-test of the DTN generator (small domain, ~1 min)
#     -> 00 regenerate the D25Q200 cache (24 workers)
#     -> 02 pilot (3 real cases; fails if any case errors)
#     -> 03 acceleration test, Gamma = 0.1, 0.25, 0.5, 1.0, PTFE at R, 34 frequencies
#     -> 03b Gamma = 1.5, 2.0 (flagged unphysical; a failure here blocks nothing)
#     -> 06 numerics check at Gamma = MAIN_GAMMA (Euler; 30 steps/period)   [needs 03; 03b: any]
#     -> 04a gate: does the curve depend on Gamma?                          [needs 03; 06: any]
#          independent -> continue.  dependent or not enough data -> the gate FAILS, 04/05 are
#          cancelled, you get an email, and the data stay for analysis.
#     -> 04 sweep at radius R, four beads, Gamma = MAIN_GAMMA              [needs gate]
#     -> 05 sweep at radius rc, four beads                                  [needs gate; 04: any]
#     -> notify (final): runs after everything ended, mails a status summary.
#
# Email (to MAILTO): Slurm mails END, FAIL and TIME_LIMIT for every job; the final notify job
# also sends a summary with the content (job states, peak memory, gate report, quality counts).
#
# Environment overrides: MAIN_GAMMA (default 0.5; must be one of the acceleration-test values),
# MAILTO, SCHEME, DOMAIN ("25,200"), TRES (60).

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../../.."

DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1

MAILTO="${MAILTO:-elvis_vera@brown.edu}"
MAIN_GAMMA="${MAIN_GAMMA:-0.5}"
ACC_GAMMAS="0.1,0.25,0.5,1.0"
HI_GAMMAS="1.5,2.0"
B=julia/cluster/beads
MAIL=(--mail-type=END,FAIL,TIME_LIMIT --mail-user="$MAILTO")

[[ ",$ACC_GAMMAS," == *",$MAIN_GAMMA,"* ]] || { echo "MAIN_GAMMA=$MAIN_GAMMA is not in the acceleration-test list ($ACC_GAMMAS)." >&2; exit 1; }
[[ -f "$B/run_bead_sweep.jl" ]] || { echo "Run from the repository (files not found)." >&2; exit 1; }

CAMP="output/beads/campaign_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$CAMP"
: > "$CAMP/job_ids.txt"
DRYN=0

# sub <name> <dependency or ""> <script> [VAR=value ...]  -> prints the job id
sub() {
    local name="$1" dep="$2" script="$3"; shift 3
    local args=(--parsable --job-name="$name" "${MAIL[@]}")
    [[ -n "$dep" ]] && args+=(--dependency="$dep" --kill-on-invalid-dep=yes)
    local id
    if [[ "$DRY" == 1 ]]; then
        DRYN=$(( $(wc -l < "$CAMP/job_ids.txt") + 1 )); id="DRY$DRYN"
        echo "DRY  env CAMP=$CAMP $* | sbatch ${args[*]} $script" >&2
    else
        id=$(env CAMP="$CAMP" MAILTO="$MAILTO" MAIN_GAMMA="$MAIN_GAMMA" "$@" sbatch "${args[@]}" "$script")
    fi
    echo "$name $id" >> "$CAMP/job_ids.txt"
    echo "$id"
}

echo "Campaign directory: $CAMP   (dry run: $DRY)"
J_start=$(sub bd-notify-start "" $B/notify.sh MODE=start)

if [[ -f julia/data/dtn_cache/D25_Quant200_nr2500.jld2 ]]; then
    echo "The D25Q200 cache file already exists: skipping stages 00a and 00." >&2
    DEP_PILOT=""
else
    J_self=$(sub bd-dtn-selftest "" $B/00a_dtn_selftest.sh)
    J_dtn=$(sub bd-dtn-D25Q200 "afterok:$J_self" $B/00_dtn_D25Q200.sh)
    DEP_PILOT="afterok:$J_dtn"
fi

J_pilot=$(sub bd-pilot "$DEP_PILOT" $B/02_pilot.sh RUN_DIR="$CAMP/02_pilot")
J_acc=$(sub bd-accel "afterok:$J_pilot" $B/03_accel.sh RUN_DIR="$CAMP/03_accel" GAMMAS="$ACC_GAMMAS")
J_hi=$(sub bd-accel-gt1 "afterok:$J_acc" $B/03_accel.sh RUN_DIR="$CAMP/03b_accel_gt1" GAMMAS="$HI_GAMMAS")
J_num=$(sub bd-numerics "afterok:$J_acc,afterany:$J_hi" $B/06_numerics_check.sh RUN_DIR="$CAMP/06_numerics" GAMMAS="$MAIN_GAMMA")
J_gate=$(sub bd-gate "afterok:$J_acc,afterany:$J_num" $B/04a_gate.sh)
J_R=$(sub bd-sweep-R "afterok:$J_gate" $B/04_sweep_R.sh RUN_DIR="$CAMP/04_sweep_R" GAMMAS="$MAIN_GAMMA")
J_rc=$(sub bd-sweep-rc "afterok:$J_gate,afterany:$J_R" $B/05_sweep_rc.sh RUN_DIR="$CAMP/05_sweep_rc" GAMMAS="$MAIN_GAMMA")

ALL=$(awk '$1!="bd-notify-start"{printf "%s%s", (n++?":":""), $2}' "$CAMP/job_ids.txt")
J_final=$(sub bd-notify-final "afterany:$ALL" $B/notify.sh MODE=final)

echo
echo "Submitted. Job ids are in $CAMP/job_ids.txt:"
cat "$CAMP/job_ids.txt"
echo
echo "Final summary (also mailed to $MAILTO): $CAMP/CAMPAIGN_STATUS_final.txt"
echo "Gate report:                            $CAMP/gate/gate_report.md"
