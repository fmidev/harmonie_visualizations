#!/bin/sh

# Run python script to plot Harmonie weather charts
# Sonja Murto / Feb 20180221
# 
# Modified by MKa
# - structure changes / 20180223
# - added file name suffix, chnage of timestep parameter numbers / 20180601
# - Python 3 system version 20220119

appPath=${1-"noneGiven"}          # script directory
Cycle=${2-"noneGiven"}            # forecast cycle : HH
fnSuffix=${3-"noneGiven"}         # forecast file suffix : member number (e.g. '_mbr000' or '' if none specified)
N0=${4-"noneGiven"}               # timesteps to plot : start step
N1=${5-"noneGiven"}               # timesteps to plot : stop step

echo $Cycle
echo $fnSuffix
echo $N0
echo $N1

echo ' '
echo '== Starting'
echo ' '

for K in `seq $N0 $N1`
do
    echo " "
    echo " "
    printf -v ZK '%03d' "$K"      # pad with leading zeros to HHH format
    echo "Plotting timestep : $ZK"
    echo $Cycle
    echo $ZK
    echo $fnDuffix
    echo "Scandinavia:"
    python3 $appPath/WMP.py $Cycle $ZK $fnSuffix
    # echo "Southern Finland:"
    # python3 $appPath/WMP_southFIN.py $Cycle $ZK $fnSuffix
    # echo "Northern Finland:"
    # python3 $appPath/WMP_northFIN.py $Cycle $ZK $fnSuffix
done

echo ' '
echo '== Done plotting'
echo ' '
