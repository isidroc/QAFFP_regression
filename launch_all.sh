#!/bin/bash


for dataset in A2780 A2a ABL1 Acetylcholinesterase Androgen Aurora-A B-raf Cannabinoid Carbonic Caspase CCRF-CEM Coagulation COX-1 COX-2 Dihydrofolate Dopamine DU-145 Ephrin erbB1 Estrogen Glucocorticoid Glycogen HCT-15 HERG HL-60 JAK2 K562 KB L1210 LCK LoVo MDA-MB-231 Monoamine opioid PC-3 SK-OV-3 Vanilloid  HCT-116   HeLa  HepG2  HT-29    MDA-MB-435    NCI-H460

do 

for repe in {1..50}
do

for type_fp in Morgan physchem QAFFP  QAFFP_Morgan QAFFP_physchem 
do

# 440 QAFFP
for fp_length in 440
do
for type_QAFFP in rv_no_cp_   cp_binary_
do

if [ ! -f results/test_${type_fp}_${dataset}_${repe}_0.3 ];
then

echo results/test_${type_fp}_${dataset}_${repe}_0.3

bsub -M 10G -o o.log -e e.log -J job_name "python run_QAFFP.py ${dataset} ${type_fp} ${repe}  ${type_QAFFP}${fp_length}  ${fp_length}"
fi

done
done

# 1360 QAFFP
for fp_length in 1360
do
for type_QAFFP in rv_no_cp_
do

if [ ! -f results/test_${type_fp}_${dataset}_${repe}_0.3 ];
then

echo results/test_${type_fp}_${dataset}_${repe}_0.3

bsub -M 10G -o o.log -e e.log -J job_name "python run_QAFFP.py ${dataset} ${type_fp} ${repe}  ${type_QAFFP}${fp_length}  ${fp_length}"
fi

done
done

done
done
done
