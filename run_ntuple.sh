#!/bin/bash

mhc=10000

./create_ntuple ntuple_list_files/input_signal_s_mhc100.par
./create_ntuple ntuple_list_files/input_signal_t2_mhc100.par
./create_ntuple ntuple_list_files/input_signal_t3_mhc100.par
./create_ntuple ntuple_list_files/input_signal_tw_mhc100.par


gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_back_signal_s.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_back_signal_t2.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_back_signal_t3.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_back_signal_tw.par

./create_ntuple ntuple_list_files/input_back_signal_s.par
./create_ntuple ntuple_list_files/input_back_signal_t2.par
./create_ntuple ntuple_list_files/input_back_signal_t3.par
./create_ntuple ntuple_list_files/input_back_signal_tw.par

gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_ttjets_lept.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_ttjets_semilept.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_ttjets_hadr.par

./create_ntuple ntuple_list_files/input_ttjets_lept.par
./create_ntuple ntuple_list_files/input_ttjets_semilept.par
./create_ntuple ntuple_list_files/input_ttjets_hadr.par

gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_wc0jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_wc1jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_wc2jets.par

./create_ntuple ntuple_list_files/input_wc0jets.par
./create_ntuple ntuple_list_files/input_wc1jets.par
./create_ntuple ntuple_list_files/input_wc2jets.par

gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_wbb0jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_wbb1jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_wbb2jets.par

./create_ntuple ntuple_list_files/input_wbb0jets.par
./create_ntuple ntuple_list_files/input_wbb1jets.par
./create_ntuple ntuple_list_files/input_wbb2jets.par

gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_wmp0jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_wmp1jets.par
gsed -i "s/MH2 = [0-9]*/MH2 = ${mhc}/g" ntuple_list_files/input_wmp2jets.par

./create_ntuple ntuple_list_files/input_wmp0jets.par
./create_ntuple ntuple_list_files/input_wmp1jets.par
./create_ntuple ntuple_list_files/input_wmp2jets.par

exit

#./skimming_events ntuple_list_files/input_wmp3jets.par
#./skimming_events ntuple_list_files/input_wc3jets.par
#./skimming_events ntuple_list_files/input_wbb3jets.par


exit

./skimming_events input_qcd_250To500.par
./skimming_events input_qcd_1000ToInf.par
./skimming_events input_qcd_500To1000.par
./skimming_events input_qcd_100To250.par 
