#!/bin/sh
# import MSP 
../../JUMP/bin/jump -msp2h5 -v MoNA-export-ReSpect.msp.bz2 MoNA-export-ReSpect.h5
# run speectral search and store results in "results.csv"
../../JUMP/bin/jump -sl jump_sl.params sample.dtas --output=- --output-format=csv | tee results.csv