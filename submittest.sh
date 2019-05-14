#!/bin/bash

nevts=5000000
#nevts=200000

#for i in `seq 21 100`; do
#for i in `seq 1 100`; do
for i in `seq 1 100`; do
    #    ./dalitz_D3K_sig /home/sunla/ToyMC/FastMC_DKKK_SimpleIsobar_1M_0.root > logfiles/mctest_ex1_oldnorm_300_${i}.log
    #    ./dalitz_D3K_sig /home/sunla/ToyMC/FastMC_DKKK_SimpleIsobar_1M_0.root > logfiles/mctest_ex1_newnorm_1e6_${i}.log
    echo "Job #$i ..............................."
#    ./dalitz_D3K_sig FastMC_DKKK_SimpleIsobar_ex2_1M_0.root > logfiles/mctest_ex2_newnorm_1e5_${i}.log
#    ./dalitz_D3K_sig FastMC_DKKK_SimpleIsobar_ex2_1M_0.root > logfiles/mctest_ex2_newnorm_5e5_${i}.log
#    ./dalitz_D3K_sig FastMC_DKKK_SwaveOnly_1M_0.root > logfiles/mctest_swave_newnorm_1e5_${i}.log
#   ./dalitz_D3K_sig FastMC_sigtoy_0.root > logfiles/test100K_b50_nb1000_t_a10_r${i}.log
#./dalitz_Ds3pi_PWA_Sig  FastMC_DsPiPiPi_PWA_11K_30pt_0.root $nevts i | tee toylogs/toysigfit_n${nevts}_startv${i}.log
#./dalitz_Ds3pi_PWA_SB  FastMC_DsPiPiPi_PWA_11K_30pt_0.root -n $nevts -r $i  | tee toylogs/toysigsbfit_n${nevts}_startv${i}.log
#./dalitz_Ds3pi_PWA Data_Ds2pipipi_sw_sigreg.root -e -r $i  | tee toylogs/datafit_paraEff_histBkg_c_n${nevts}_startv_${i}_altspline.log
./dalitz_Ds3pi_PWA Data_Ds2pipipi_sw_sigreg.root -e -r $i -n $nevts  | tee toylogs/datafit_paraEff_histBkg_c_n${nevts}_startv_${i}.log
#./dalitz_Ds3pi_PWA Data_Ds2pipipi_sw_sigreg.root $nevts 1 1 | tee toylogs/datafit_paraEff_paraBkg_n${nevts}_r${i}.log
done
