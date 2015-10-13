#!/bin/bash
#SBATCH -J MatchTest
#SBATCH -A nesi00186
#SBATCH --time=01:00:00
#SBATCH --mem=7G
#SBATCH --nodes=1
#SBATCH --output=matchout_%a.txt
#SBATCH --error=matcherr_%a.txt
#SBATCH --cpus-per-task=16
#SBATCH --array=0-12

module load OpenCV/2.4.9-intel-2015a
module load ObsPy/0.10.3rc1-intel-2015a-Python-2.7.9
module load joblib/0.8.4-intel-2015a-Python-2.7.9

srun python2.7 LFEsearch.py --splits 13 --instance $SLURM_ARRAY_TASK_ID
#srun python2.7 LFEsearch.py
