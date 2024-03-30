#!/bin/bash
########################################################
# job_builder.sh
# Job Arrays without the Resource Manager # Version 1.0.0 ########################################################
# James Joseph Balamuta
# balamut2@illinois.edu
########################################################
# ## Example
#
# # Allow the builder script to work on the file system # chmod +x job_builder.sh # # # Run the job builder # ./job_builder.sh ########################################################


### Builds the job index
# Create a sequential range
array_values=`seq 50`

# Launch the job and then remove the temporarily created qsub file.
for i in $array_values
do 
# This submits the single job to the resource manager
sbatch BCPlaidSimulations.sub $i

done
