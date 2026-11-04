#!/bin/bash
#SBATCH --job-name=bd-notify
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M
#SBATCH --time=00:05:00
#SBATCH --output=output/beads/99_report/%x-%j.out
#SBATCH --error=output/beads/99_report/%x-%j.err
#
# Campaign notifications, in two modes (environment MODE):
#   MODE=start : tiny "campaign started" message (tests that mail from a compute node works).
#   MODE=final : runs after every other job has ended (afterany), writes
#                <campaign>/CAMPAIGN_STATUS_final.txt (state, exit code, time and peak memory of
#                every job, the gate report, quality counts of the sweeps) and mails it.
# It never fails the campaign: every command is allowed to fail. The mail is sent with
# /usr/sbin/sendmail; if that does not work on the node, Slurm's own END mail still arrives
# (the job is submitted with --mail-type), it just lacks the content.

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?run with sbatch}"
: "${CAMP:?set CAMP}"
MODE="${MODE:-final}"
MAILTO="${MAILTO:-elvis_vera@brown.edu}"
OUT="$CAMP/CAMPAIGN_STATUS_${MODE}.txt"

{
    echo "Bead campaign: $MODE   ($(date))"
    echo "Campaign directory: $PWD/$CAMP"
    if [[ "$MODE" == "final" && -f "$CAMP/job_ids.txt" ]]; then
        echo
        echo "== Jobs =="
        printf '%-16s %-10s %-14s %-9s %-10s %s\n' stage jobid state exit elapsed peak_mem
        n_bad=0
        while read -r name id; do
            [[ "$id" == "$SLURM_JOB_ID" ]] && continue
            line=$(sacct -X -n -P -j "$id" -o State,ExitCode,Elapsed 2>/dev/null | head -1)
            state=$(cut -d'|' -f1 <<<"$line"); ec=$(cut -d'|' -f2 <<<"$line"); el=$(cut -d'|' -f3 <<<"$line")
            mem=$(sacct -n -P -j "$id" -o MaxRSS 2>/dev/null | grep -E '[0-9]' | sort -h | tail -1)
            printf '%-16s %-10s %-14s %-9s %-10s %s\n' "$name" "$id" "${state:-?}" "${ec:-?}" "${el:-?}" "${mem:--}"
            [[ "${state%% *}" == "COMPLETED" ]] || { [[ "$name" == bd-notify-start ]] || n_bad=$((n_bad+1)); }
        done < "$CAMP/job_ids.txt"
        echo
        [[ -f "$CAMP/gate/gate_report.md" ]] && { echo "== Gate report =="; cat "$CAMP/gate/gate_report.md"; echo; }
        for d in "$CAMP"/0*/ ; do
            f="$d/gain_phase_all.csv"
            [[ -f "$f" ]] || continue
            echo "== Quality counts: $(basename "$d") =="
            awk -F, 'NR==1{for(i=1;i<=NF;i++) if($i=="quality") q=i; next} {c[$q]++} END{for(k in c) printf "  %-16s %d\n", k, c[k]}' "$f" | sort
        done
    fi
} > "$OUT" 2>&1

if [[ "$MODE" == "final" ]]; then
    if [[ "${n_bad:-1}" -eq 0 ]]; then SUBJECT="[bead campaign] ALL STAGES COMPLETED"; else SUBJECT="[bead campaign] ATTENTION: ${n_bad:-?} stage(s) did not complete"; fi
else
    SUBJECT="[bead campaign] started"
fi
{ echo "To: $MAILTO"; echo "Subject: $SUBJECT"; echo; cat "$OUT"; } | /usr/sbin/sendmail -t \
    && echo "sendmail: message handed over" || echo "sendmail failed (non-fatal)"
exit 0
