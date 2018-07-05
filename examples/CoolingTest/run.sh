#!/bin/bash

# Get the cooling tables if necessary
if [ ! -e coolingtables ]
then
    echo "Fetching EAGLE cooling tables."
    ./getTables.sh
fi

# Run the test at z=1. with half-solar abundance
./cooling_test 1. 0.5

# Plot the results
python plot_rates.py
