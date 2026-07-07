#!/bin/bash
timestamp=$(date +%Y_%m_%d)

SCRIPTS_DIR="${HOME}/data/mf_2020_hap/data_analysis/2025/pathogendetection_2025_12/scripts"
RAXML_WORKING_DIR="${HOME}/data/mf_2020_hap/data_analysis/2025/pathogendetection_2025_12/analysis/case_Pub_data_resp_P07_Ehormaechei-c1_240825/${timestamp}_P07_Ehormaechei-c1_240825_gubbinsdefault_ambigALL/raxml_bootstrapped_tree1000/"
RAXML_ALIGNMENT_FILE="${HOME}/data/mf_2020_hap/data_analysis/2025/pathogendetection_2025_12/analysis/case_Pub_data_resp_P07_Ehormaechei-c1_240825/${timestamp}_P07_Ehormaechei-c1_240825_gubbinsdefault_ambigALL/snp_msa_all.fa"
RAXML_TREE_PREFIX="2026-05-30_17-45-12_ML_GTRG_snp_msa_all"
RAXML_OUTGROUP="Equasihormaechei_GCF_004331385_1"

## activate conda environment 
activate_conda_path=$(echo $CONDA_EXE | sed 's#bin/conda#bin/activate#g')
source ${activate_conda_path} anapy_312

cd ${SCRIPTS_DIR}

echo "Starting analysis py"
date

## run python unbuffered to get immediate stdout prints
python -u hap2020_analysispy_public_data_logan_v1_2_gubbins.py

date

echo "Starting RaxML"

source ${activate_conda_path} raxml

mkdir -p ${RAXML_WORKING_DIR}
cd ${RAXML_WORKING_DIR}

echo "full tree, 1000 bootstraps (Felsenstein bootstrap and transfer bootstrap expectation)"
raxml-ng --all --threads 64 --seed 123 --msa ${RAXML_ALIGNMENT_FILE} --outgroup ${RAXML_OUTGROUP} --model GTR+G --prefix ${RAXML_TREE_PREFIX}_1000bs_ML_all_s123 --bs-trees 1000 --bs-metric fbp,tbe
