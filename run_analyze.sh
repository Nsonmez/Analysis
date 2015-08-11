#!/bin/bash

make analyze_events

mhc=10000

./analyze_events analyze_list_files/input_signal_s_mhc100.par
./analyze_events analyze_list_files/input_signal_t2_mhc100.par
./analyze_events analyze_list_files/input_signal_t3_mhc100.par
./analyze_events analyze_list_files/input_signal_tw_mhc100.par


gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_back_signal_s.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_back_signal_t2.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_back_signal_t3.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_back_signal_tw.par

./analyze_events analyze_list_files/input_back_signal_s.par
./analyze_events analyze_list_files/input_back_signal_t2.par
./analyze_events analyze_list_files/input_back_signal_t3.par
./analyze_events analyze_list_files/input_back_signal_tw.par

gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_ttjets_lept.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_ttjets_semilept.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_ttjets_hadr.par

./analyze_events analyze_list_files/input_ttjets_lept.par
./analyze_events analyze_list_files/input_ttjets_semilept.par
./analyze_events analyze_list_files/input_ttjets_hadr.par

gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_wc0jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_wc1jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_wc2jets.par

./analyze_events analyze_list_files/input_wc0jets.par
./analyze_events analyze_list_files/input_wc1jets.par
./analyze_events analyze_list_files/input_wc2jets.par

gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_wbb0jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_wbb1jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_wbb2jets.par

./analyze_events analyze_list_files/input_wbb0jets.par
./analyze_events analyze_list_files/input_wbb1jets.par
./analyze_events analyze_list_files/input_wbb2jets.par

gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_wmp0jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_wmp1jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" analyze_list_files/input_wmp2jets.par

./analyze_events analyze_list_files/input_wmp0jets.par
./analyze_events analyze_list_files/input_wmp1jets.par
./analyze_events analyze_list_files/input_wmp2jets.par

exit

#./skimming_events analyze_list_files/input_wmp3jets.par
#./skimming_events analyze_list_files/input_wc3jets.par
#./skimming_events analyze_list_files/input_wbb3jets.par


exit

./skimming_events input_qcd_250To500.par
./skimming_events input_qcd_1000ToInf.par
./skimming_events input_qcd_500To1000.par
./skimming_events input_qcd_100To250.par 
