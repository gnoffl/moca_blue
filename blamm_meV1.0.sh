#!/bin/bash

# Define variables
out_dir="random_sequences"
motifs_file="out/rdf5_MSR_ARTHM1X0.75dC7K25f_cwm-motifs.jaspar"
sequences="ref_seq/sequences_random.mf"
pt_score=0.0001


# Parse input arguments
if [ $# -gt 5 ]; then
  echo "Usage: $0 [out_dir] [motifs_file] [sequences] [pt_score]"
  exit 1
fi
if [ $# -ge 1 ]; then
  out_dir="$1"
fi
if [ $# -ge 2 ]; then
  motifs_file="$2"
fi
if [ $# -ge 3 ]; then
  sequences="$3"
fi
if [ $# -ge 4 ]; then
  pt_score="$4"
fi

#
# Start recording runtime and resource information
start_time=$(date +%s)
start_resources=$(ps -o pid,%cpu,%mem,vsz,rss,tty,stat,start_time --no-headers $$)


# Run the commands with variables
blamm dict "$sequences"
blamm hist -e "$motifs_file" "$sequences" #-e generated empirical PWM scores
blamm scan -rc -pt "$pt_score" "$motifs_file" "$sequences"

#
echo "CLEANING UP"

mkdir -p $out_dir/histograms
mv hist_* $out_dir/histograms/
mv occurrences.txt $out_dir/
mv PWMthresholds.txt $out_dir/

#
# Stop recording runtime and resource information
end_time=$(date +%s)
end_resources=$(ps -o pid,%cpu,%mem,vsz,rss,tty,stat,start_time --no-headers $$)

# Calculate runtime
runtime=$((end_time - start_time))
# Print runtime and resource information
echo "Script runtime: $runtime seconds"
echo "Resource usage:"
echo "$start_resources" | awk '{print "Start:", $0}'
echo "$end_resources" | awk '{print "End:", $0}'
