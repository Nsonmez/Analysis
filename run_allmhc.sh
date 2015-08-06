#!/bin/bash


for channel in s t2 t3 tw
do

    echo "now is making analysis for $channel"
#./skimming_eventsv5 list_files/input_signal_${channel}_mhc80.par
    ./skimming_eventsv5 list_files/input_signal_${channel}_mhc85.par
    ./skimming_eventsv5 list_files/input_signal_${channel}_mhc90.par
    ./skimming_eventsv5 list_files/input_signal_${channel}_mhc100.par
    ./skimming_eventsv5 list_files/input_signal_${channel}_mhc110.par
    ./skimming_eventsv5 list_files/input_signal_${channel}_mhc120.par
    ./skimming_eventsv5 list_files/input_signal_${channel}_mhc130.par

done


lumi=20ifb

hadd signal_mhc85_${lumi}.root  singletop_*_mhc85.root
hadd signal_mhc90_${lumi}.root  singletop_*_mhc90.root
hadd signal_mhc100_${lumi}.root singletop_*_mhc100.root
hadd signal_mhc110_${lumi}.root singletop_*_mhc110.root
hadd signal_mhc120_${lumi}.root singletop_*_mhc120.root
hadd signal_mhc130_${lumi}.root singletop_*_mhc130.root


mkdir filtered_events/signal_for_${lumi}
mv singletop_*_mhc*.root filtered_events/signal_for_${lumi}
mv signal_mhc*_${lumi}.root filtered_events/

echo "\n the analysis is done  \n"
#exit 0

