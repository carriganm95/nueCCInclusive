#!/usr/bin/bash

dataDir="/pnfs/sbn/data/sbn_fd/poms_production/2025A_icarus_NuMI_MC/FHC_NuMI/mc/reconstructed/icaruscode_v09_89_01_02p02/flatcaf/" #FHC MC Sample
#dataDir="/pnfs/sbn/data/sbn_fd/poms_production/2025A_ICARUS_NuMI_RHC_MC/NuMI_RHC/mc/reconstructed/icaruscode_v09_89_02_00p01/flatcaf/" #RHC MC Sample
#dataDir="/pnfs/sbn/data/sbn_fd/poms_production/data/Run2reprocess/reconstructed/icaruscode_v09_89_01_02p02/numimajority/flatcaf_prescaled/" #data 10% sample
#dataDir="/pnfs/sbn/data/sbn_fd/poms_production/data/Run2reprocess/reconstructed/icaruscode_v09_89_01_02p02/numimajority/flatcaf_unblind/" # run 2 unblided NUMI data
#dataDir="/exp/icarus/data/users/rtriozzi/mc/numi_FRFIX/"
#dataDir="/pnfs/sbn/data/sbn_fd/poms_production/data/Run2reprocess/reconstructed/icaruscode_v09_89_01_02p02/offbeamnumimajority/flatcaf_prescaled/" #offbeam data
#dataDir="/pnfs/icarus/scratch/users/micarrig/showerEnergyCalCaf/outputs/"
#dataDir="/pnfs/sbn/data/sbn_fd/poms_production/mc/2025A_ICARUS_NuGraph2/NuMI_MC_FullTrainingSample/v10_06_00_04p02/flatcaf/*/" #NuGraph2 MC Sample (testing shower norm)
#dataDir="/pnfs/icarus/scratch/users/fwieler/evaluation_data_caf_NuGraph2_NuMI_MC/out/" #CVN MC dataset
#dataDir="/pnfs/icarus/scratch/users/micarrig/showerEnergyCalCafMC/outputs/" #shower energy cal MC dataset
#dataDir="/pnfs/icarus/scratch/users/micarrig/showerEnergyCalCafNormMC/outputs/" #pion selection MC normalized dataset
#dataDir="/pnfs/icarus/scratch/users/micarrig/showerEnergyCalCafV2/outputs/" #pion selection data dataset v2
#dataDir='/pnfs/icarus/scratch/users/micarrig/showerEnergyCalCafNorm/outputs/' # pion selection data normalized dataset
fileStr="*.flat.caf*.root"
#fileStr="/*Unblind.DONOTLOOK.flat.caf.root"
#fileStr="concat_NuMI_MC_FRFIX_"
outputStr="nueOutputs/test/"
#exe=pionSelection.C
exe="nuCCInclusive_MC.C"

centralProd=true

if [ ! -d "$outputStr" ]; then
    mkdir -p "$outputStr"
fi

#for i in {0..9} {a..f}; do
for i in {0..0}; do
    export BEARER_TOKEN_FILE=/tmp/bt_u$(id -u) && htgettoken -a htvaultprod.fnal.gov -i icarus
    #thisJob="${dataDir}0${i}/${fileStr}"
    #thisJob="${dataDir}${fileStr}${i}.root"
    if [ "$centralProd" = true ] ; then
        hex=$(printf "%02x" "$i")
        thisJob="${dataDir}/$hex/*/${fileStr}"
        outputFile="${outputStr}/out_${hex}.root"
        echo "Starting to run... $hex $thisJob with output file $outputFile"
    else
        if (( i < 10 )); then
            ii="0$i"
        else
            ii="$i"
        fi
        thisJob="${dataDir}${ii}/${fileStr}"
        outputFile="${outputStr}/out_${ii}.root"
        echo "Starting to run... $ii $thisJob with output file $outputFile"
    fi
    cmd="cafe -bq $exe false \"$thisJob\" $outputFile"
    echo $cmd
    eval $cmd
done
