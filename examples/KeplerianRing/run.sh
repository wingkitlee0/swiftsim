#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e initial_conditions.hdf5 ]
then
    echo "Generating initial conditions for the keplerian ring example..."
    python3 generate_ics.py -m=667.428 -pm=1
fi

rm -rf keplerian_ring_*.hdf5
../swift -g -s -t 1 keplerian_ring.yml 2>&1 | tee output.log
