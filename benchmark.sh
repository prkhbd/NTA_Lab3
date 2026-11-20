
#!/bin/bash

set -u

RUNS=5
THREADS=4

CSV="benchmark_results.csv"
LOG_DIR="benchmark_logs"

if [ ! -x "./index_calculus" ] || [ ! -x "./index_calculus_parallel" ]; then
    echo "Compile both programs first."
    exit 1
fi

mkdir -p "$LOG_DIR"

echo "p,alpha,beta,B,threads,run,version,total_seconds" > "$CSV"

tests=(
    "101 2 10 10"
    "1009 11 10 15"
    "10007 5 10 20"
    "100003 2 10 25"
    "1000003 2 10 30"
    "100000007 5 10 40"
        "999999937 11 10 50"
)

for test in "${tests[@]}"; do
    read -r p alpha beta B <<< "$test"

    echo
    echo "=========================================="
    echo "Testing p=$p, alpha=$alpha, beta=$beta, B=$B"
    echo "=========================================="

    for run in $(seq 1 "$RUNS"); do
        echo
        echo "Run $run/$RUNS: sequential"

        input="$p
$alpha
$beta
$B"

        if ./index_calculus <<< "$input" > "$LOG_DIR/seq_${p}_${run}.txt"; then
            total=$(awk '/Total:/ {print $2}' "$LOG_DIR/seq_${p}_${run}.txt")

            if [ -n "$total" ]; then
                echo "$p,$alpha,$beta,$B,0,$run,sequential,$total" >> "$CSV"
                echo "Sequential total: $total s"
            else
                echo "Could not extract sequential time."
            fi
        else
            echo "Sequential run failed."
        fi

        echo
        echo "Run $run/$RUNS: parallel"

        parallel_input="$p
$alpha
$beta
$B
$THREADS"

        if ./index_calculus_parallel <<< "$parallel_input" > "$LOG_DIR/par_${p}_${run}.txt"; then
            total=$(awk '/Total:/ {print $2}' "$LOG_DIR/par_${p}_${run}.txt")

            if [ -n "$total" ]; then
                echo "$p,$alpha,$beta,$B,$THREADS,$run,parallel,$total" >> "$CSV"
                echo "Parallel total: $total s"
            else
                echo "Could not extract parallel time."
            fi
        else
            echo "Parallel run failed."
        fi
    done
done

echo
echo "=========================================="
echo "Benchmark completed."
echo "Results: $CSV"
echo "Logs: $LOG_DIR"
echo "=========================================="
