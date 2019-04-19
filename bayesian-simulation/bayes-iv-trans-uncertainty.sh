#!/bin/sh
#
#

# request bash shell as shell for job
#$ -S /bin/bash
#
# set working directory
cd ~/Thomas/bayes-iv-model

# set local variables
outfile='bayes-iv-trans-uncertainty.out'.$$
        
# execute command(s)
R CMD BATCH bayes-iv-trans-uncertainty.R $outfile

# notify when script has finished
echo "Outfile is `pwd`/$outfile" | mail -v -s "bayes-iv-trans-uncertainty.sh finished" jiaqiyin@uw.edu
