#!/usr/bin/bash

dataDir="/pnfs/sbn/data/sbn_fd/poms_production/2025A_icarus_NuMI_MC/FHC_NuMI/mc/reconstructed/icaruscode_v09_89_01_02p02/flatcaf/" #FHC MC Sample
#dataDir="/pnfs/sbn/data/sbn_fd/poms_production/2025A_ICARUS_NuMI_RHC_MC/NuMI_RHC/mc/reconstructed/icaruscode_v09_89_02_00p01/flatcaf/" #RHC MC Sample
#dataDir="/pnfs/sbn/data/sbn_fd/poms_production/data/Run2reprocess/reconstructed/icaruscode_v09_89_01_02p02/numimajority/flatcaf_prescaled/" #data 10% sample
#dataDir="/exp/icarus/data/users/rtriozzi/mc/numi_FRFIX/"
#dataDir="/pnfs/sbn/data/sbn_fd/poms_production/data/Run2reprocess/reconstructed/icaruscode_v09_89_01_02p02/offbeamnumimajority/flatcaf_prescaled/" #offbeam data
#dataDir="/pnfs/icarus/scratch/users/micarrig/showerEnergyCalCaf/outputs/"
#dataDir="/pnfs/sbn/data/sbn_fd/poms_production/mc/2025A_ICARUS_NuGraph2/NuMI_MC_FullTrainingSample/v10_06_00_04p02/flatcaf/*/" #NuGraph2 MC Sample (testing shower norm)
#dataDir="/pnfs/icarus/scratch/users/fwieler/evaluation_data_caf_NuGraph2_NuMI_MC/out/" #CVN MC dataset
#fileStr="output_*.flat.caf.root"
fileStr="/*/*.flat.caf*.root"
#fileStr="concat_NuMI_MC_FRFIX_"
outputStr="nueOutputs/mcV5/"
#exe=pionSelection.C
exe="nuCCInclusive_MC.C"

if [ ! -d "$outputStr" ]; then
    mkdir -p "$outputStr"
fi

for i in {0..9} {a..f}; do
#for i in {1..9}; do
    export BEARER_TOKEN_FILE=/tmp/bt_u$(id -u) && htgettoken -a htvaultprod.fnal.gov -i icarus
    #thisJob="${dataDir}0${i}/${fileStr}"
    #thisJob="${dataDir}${fileStr}${i}.root"
    thisJob="${dataDir}${i}?/${fileStr}"
    outputFile="${outputStr}_${i}.root"
    echo "Starting to run... $i $thisJob with output file $outputFile"
    cafe -bq $exe false "$thisJob" "$outputFile"
done