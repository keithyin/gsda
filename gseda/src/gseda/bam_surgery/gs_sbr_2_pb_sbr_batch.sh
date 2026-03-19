#!/bin/bash

for filep in /data1/ccs_data/speed-test/benchmark-data/*.subreads.bam; do
    output_name="${filep%.bam}.pb.bam"
    echo "$filep --> $output_name"
    python gs_sbr_2_pb_sbr.py $filep $output_name
   
done