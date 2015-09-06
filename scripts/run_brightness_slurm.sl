#!/bin/bash
#SBATCH -J PythonTest
#SBATCH -A nesi00186
#SBATCH --time=00:01:00
#SBATCH --mem-per-cpu=1024
#SBATCH --output=stdout.txt
#SBATCH --error=stderr.txt

module load OpenCV/2.4.9-intel-2015a
module load ObsPy/0.10.3rc1-intel-2015a-Python-2.7.9
module load joblib/0.8.4-intel-2015a-Python-2.7.9

python2.7 Test_install.py 2>&1
