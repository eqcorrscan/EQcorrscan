#!/bin/env bash
#===============================================================================
#
#          FILE: run_all_time_matches.sh
#
#         USAGE: ./run_all_time_matches.sh
#
#   DESCRIPTION:
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (),
#  ORGANIZATION:
#       CREATED: 20/07/15 02:15:10 UTC
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

# python LFEsearch.py > all_detections.log 2>&1 &
splits=10 # Number of instances
i=0 # Beginning iterator
while ((i<$splits)); do
    python LFEsearch.py --splits $splits --instance $i > instance_$i.log 2>&1 &
    let i++
done

# python LFEsearch.py --splits 2 --instance 0 > instance_4.log 2>&1 &
# python LFEsearch.py --splits 2 --instance 1 > instance_11.log 2>&1 &
# python LFEsearch.py --splits 2 --instance 1 > instance_2.log 2>&1
# python LFEsearch.py --splits 5 --instance 2 > instance_3.log 2>&1 &
# python LFEsearch.py --splits 5 --instance 3 > instance_4.log 2>&1 &
# python LFEsearch.py --splits 5 --instance 4 > instance_5.log 2>&1
# python LFEsearch.py --splits 6 --instance 5 > instance_6.log 2>&1
