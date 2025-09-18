#!/usr/bin/env bash
# Usage: ./make_sqs_w5ta4nb.sh  # requires ATAT binaries (corrdump, mcsqs) on PATH

set -euo pipefail

WORKDIR=$(pwd)/sqs_W5Ta4Nb
mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "Working in $WORKDIR"

# 1) Create rndstr.in describing a conventional BCC cubic cell (2 atoms per conventional cell)
#    Occupancies on each lattice site: W=0.5, Ta=0.4, Nb=0.1
cat > rndstr.in <<'EOF'
1.0 1.0 1.0 90 90 90
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
0.0 0.0 0.0 W=0.5,Ta=0.4,Nb=0.1
0.5 0.5 0.5 W=0.5,Ta=0.4,Nb=0.1
EOF

# 2) Create sqscell.out to make a 5x1x1 supercell => 2 * 5 = 10 atoms (composition 5:4:1)
#    Adjust this matrix to change size (e.g., 5 2 1 -> 20*2=40 atoms, etc.)
cat > sqscell.out <<'EOF'
1

5 0 0
0 1 0
0 0 1
EOF

# 3) Run corrdump to generate cluster information (pairs/triplets/quadruplets up to a given distance)
#    adjust -2 -3 -4 distances as desired for more/less cluster range
echo "Running corrdump..."
corrdump -l=rndstr.in -ro -noe -nop -clus -2=1.1 -3=1.1 -4=1.1 > corrdump.log 2>&1 || {
  echo "corrdump failed — check corrdump.log"
  exit 1
}
echo "corrdump completed. log -> corrdump.log"

# 4) Run mcsqs to search for SQS on the supercell(s) listed in sqscell.out
#    -ip lets you run parallel instances (set to number of CPU cores you want), here 4 as example.
#    If you have many cores, increase -ip. You can also drop the trailing '&' lines and run sequentially.
echo "Launching mcsqs instances (parallel search). Change -ip value if you want different concurrency."
# Launch 4 instances (adjust if you prefer)
mcsqs -rc -ip=1 > mcsqs_ip1.log 2>&1 &
mcsqs -rc -ip=2 > mcsqs_ip2.log 2>&1 &
mcsqs -rc -ip=3 > mcsqs_ip3.log 2>&1 &
mcsqs -rc -ip=4 > mcsqs_ip4.log 2>&1 &

# Wait for background mcsqs jobs to finish
echo "Waiting for mcsqs background jobs..."
wait

echo "All mcsqs instances finished. Logs: mcsqs_ip*.log"

# 5) Choose the best structure among generated outputs
echo "Selecting best SQS (mcsqs -best)..."
mcsqs -best > mcsqs_best.log 2>&1 || {
  echo "mcsqs -best failed — check mcsqs_best.log"
  exit 1
}
echo "Best SQS selected. Logs -> mcsqs_best.log"

# 6) The bestsqs.out should now exist. Show a short preview:
if [ -f bestsqs.out ]; then
  echo "Preview of bestsqs.out (first 60 lines):"
  sed -n '1,60p' bestsqs.out
else
  echo "bestsqs.out not found — check logs."
  exit 1
fi

echo "Done. If you want POSCAR/VASP format, run the optional python converter provided separately."
