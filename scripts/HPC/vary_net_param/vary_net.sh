var=$1
inc=$2
jobs=$3

SUMMARY_FILE="summary.csv"

# Writes the header
if [[ -f "$SUMMARY_FILE" ]]; then
    rm -f "$SUMMARY_FILE"
fi

echo "shape_name,pi_pi_energy" > "$SUMMARY_FILE"
# write summary file
for i in $(seq 1 $jobs); do
    value=$(echo "$var + $inc * ($i-1)" | bc)
    # Prepare line to append
    LINE="$i,$value"
    echo "$LINE" >> "$SUMMARY_FILE"
done

# submit job
sbatch --array="1-$jobs" vary_net_param.sh $var $inc

