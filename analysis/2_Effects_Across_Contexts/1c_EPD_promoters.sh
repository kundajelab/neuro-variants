conda activate bedtools

# awk 'BEGIN {OFS="\t"} { $2=$2-100; $3=$3+100; print $0 }' /oak/stanford/groups/smontgom/amarder/data/epd/EPDnew_v6_hg38.bed > /oak/stanford/groups/smontgom/amarder/data/epd/EPDnew_v6_hg38_modified.bed

epd=/oak/stanford/groups/smontgom/amarder/data/epd/EPDnew_v6_hg38_modified.bed
snp=/oak/stanford/groups/smontgom/amarder/tmp/strongest_rv.bed

bedtools intersect -a $snp -b $epd > /oak/stanford/groups/smontgom/amarder/tmp/variant_input.hg38.epd.txt




