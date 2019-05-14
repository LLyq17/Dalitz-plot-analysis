#!/bin/bash

nevts=5000000
#nevts=200000

#for i in `seq 21 100`; do
#for i in `seq 1 100`; do
for i in `seq 0 99`; do
    echo "Job #$i ..............................."
./DalitzFit_Ds3pi_pwa  toyFiles/FastMC_DsPiPiPi_PWA_11K_30pt_${i}.root $nevts | tee toylogs/toysigfit_n${nevts}_toy${i}_defstart.log
#./dalitz_Ds3pi_PWA_SB  FastMC_DsPiPiPi_PWA_11K_30pt_0.root -n $nevts -r $i  | tee toylogs/toysigsbfit_n${nevts}_startv${i}.log
done
