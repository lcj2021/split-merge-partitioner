#!/bin/bash

command="$1"
output_file="$2"

# temp file: capture all output from Shell
temp_file=$(mktemp)

exec > >(tee "$temp_file") 2>&1

echo "Output redirected to $output_file"

$command >> "$output_file" 2>&1 &
pid=$!
echo "Command PID: $pid"

timestamp=$(date +"%Y%m%d%H%M%S")
log_file="memory_$timestamp.log"
echo "Timestamp    Memory (kB)" > "$log_file"

while true; do
    memory=$(awk '/VmHWM/ {print $2}' "/proc/$pid/status")
    timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "$timestamp    $memory" >> "$log_file"
    sleep 1
    
    if ! ps -p $pid > /dev/null; then
        break
    fi
done

max_memory=$(awk 'NR>1{print $3}' "$log_file" | sort -nr | head -n1)

max_memory=$(awk 'NR>1{print $3}' "$log_file" | sort -nr | head -n1)
max_memory_gb=$(awk "BEGIN { printf \"%.2f\", $max_memory / 1024 / 1024 }")

echo "Max Memory: $max_memory_gb GB" >> "$output_file"

cat "$temp_file" >> "$output_file"

rm "$temp_file" "$log_file"