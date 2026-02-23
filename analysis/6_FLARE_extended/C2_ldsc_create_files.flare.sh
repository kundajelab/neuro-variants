#!/bin/bash

# --- constants ---
path_to_ldsc="/oak/stanford/groups/smontgom/amarder/bin/ldsc"
ldsc_py="$path_to_ldsc/ldsc.py"
hapmap_dir="$path_to_ldsc/hapmap/hapmap3_snps"
geno_prefix="/oak/stanford/groups/smontgom/amarder/bin/ldsc/1kg/1000G_EUR_Phase3_plink/1000G.EUR.QC"
annot_dir="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/ldsc/annot/hg19/extended"
annot_dir="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/ldsc/annot/hg19/extended2"

# activate the LDSC conda env (Python 2 or 3 depending on your install)
source ~/.bashrc
conda activate ldsc

# --- main loop ---
for subdir in extended3; do
annot_dir="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/ldsc/annot/hg19/$subdir"
for annot_name in FLARE_fb.thres95 FLARE_fb.thres99 FLARE_fb.thres999 FLARE_heart.thres95 FLARE_heart.thres99 FLARE_heart.thres999; do
    for chrnum in {1..22}; do
        echo
        echo "==== $annot_name | chr $chrnum ===="

        hapmap="$hapmap_dir/hm.${chrnum}.snp"
        bfile="${geno_prefix}.${chrnum}"
        annot_path="${annot_dir}/${annot_name}.chr${chrnum}.annot.gz"
        out_prefix="${annot_dir}/${annot_name}.chr${chrnum}"

        # Quick checks
        [[ -f "$annot_path" ]] || { echo "Missing annot: $annot_path"; exit 1; }
        [[ -f "${bfile}.bed" ]] || { echo "Missing plink bed: ${bfile}.bed"; exit 1; }
        [[ -f "${bfile}.bim" ]] || { echo "Missing plink bim: ${bfile}.bim"; exit 1; }
        [[ -f "${bfile}.fam" ]] || { echo "Missing plink fam: ${bfile}.fam"; exit 1; }
        [[ -f "$hapmap" ]] || { echo "Missing HapMap3 list: $hapmap"; exit 1; }

        # Run LDSC
        python "$ldsc_py" \
            --l2 \
            --bfile "$bfile" \
            --ld-wind-cm 1 \
            --annot "$annot_path" \
            --thin-annot \
            --out "$out_prefix" \
            --print-snps "$hapmap"

        # Optional: quick confirm of outputs
        for suffix in ".l2.ldscore.gz" ".l2.M.gz" ".l2.log"; do
            if [[ ! -f "${out_prefix}${suffix}" ]]; then
                echo "Warning: Missing ${out_prefix}${suffix}"
            fi
        done

        echo "Done: $out_prefix"
    done
done
done
