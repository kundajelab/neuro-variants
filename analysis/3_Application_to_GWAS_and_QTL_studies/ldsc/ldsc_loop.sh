# List of traits
traitslist=("Anorexia_nervosa_Watson_2019" "ALS_Rheenen_2021" "ASD_Grove_2019" "PTSD_Nievergelt_2019" "ADHD_Demontis_2019" "Bipolar_Mullins_2021" "MajorDepression_Meng_2024" "Tourette_Yu_2019" "Alzheimers_Bellenguez_2022" "BMI_Yengo_2018" "Schizophrenia_PGCWave3_2022" "CRC" "BMI_NealePanUKB_2021" "CAD_Tcheandjieu_2022" "CAD_Aragam_2022")

# List of annotations
annotationslist=("CRC" "ameen_2022" "corces_2020" "domcke_2020" "encode_2024" "roussos_2024" "trevino_2021" "IGVFCoronaryArteriesMultiome_2024" "turner_2022")
annotationslist=("encode_2024")

cd /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/scripts/3_Application_to_GWAS_and_QTL_studies/ldsc
for annotation in ${annotationslist[@]}; do
./build_ldcts.sh $annotation
for trait in ${traitslist[@]}; do
./run_ldsc.sh $trait $annotation
echo $annotation $trait
done
done

