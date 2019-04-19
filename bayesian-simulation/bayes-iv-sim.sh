#!/bin/sh
#
#

# request bash shell as shell for job
#$ -S /bin/bash
#
# set working directory
cd ~/Thomas/bayes-iv-model

# set local variables
outfile='R-bayes-iv-sim.out'.$$

# execute command(s)
R CMD BATCH bayes-iv-sim-transparent_0918.R $outfile

# notify when script has finished
echo "Outfile is `pwd`/$outfile" | mail -v -s "bayes-iv-sim.R.sh finished" jiaqiyin@uw.edu
