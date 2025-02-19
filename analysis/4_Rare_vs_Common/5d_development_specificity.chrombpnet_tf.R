library(data.table)

df1 = fread("/oak/stanford/groups/akundaje/projects/neuro-variants/tables/rare_variants.motif_annotations.tsv",data.table=F,stringsAsFactors=F)
f="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/excitatory_neurons_dev.dataset.txt"
df2 = fread(f,data.table = F,stringsAsFactors = F)

df1 = df1[,c("variant_id","corces_2020.Cluster1.all_motifs","domcke_2020.fetal_brain.Excitatory_neurons.all_motifs")]
df.mg = merge(df2,df1,by.x="snp_id",by.y="variant_id")

df.mg$EGR1_a = grepl("EGR1",df.mg$corces_2020.Cluster1.all_motifs)
df.mg$EGR1_f = grepl("EGR1",df.mg$domcke_2020.fetal_brain.Excitatory_neurons.all_motifs)

table(df.mg$EGR1_a,df.mg$adult_specific)
table(df.mg$EGR1_f,df.mg$adult_specific)

df.mg$CTCF_a = grepl("CTCF",df.mg$corces_2020.Cluster1.all_motifs)
df.mg$CTCF_f = grepl("CTCF",df.mg$domcke_2020.fetal_brain.Excitatory_neurons.all_motifs)

table(df.mg$CTCF_a,df.mg$adult_specific)
table(df.mg$CTCF_f,df.mg$adult_specific)

df.mg$ATOH_a = grepl("ATOH",df.mg$corces_2020.Cluster1.all_motifs)
df.mg$ATOH_f = grepl("ATOH",df.mg$domcke_2020.fetal_brain.Excitatory_neurons.all_motifs)

table(df.mg$ATOH_a,df.mg$adult_specific)
table(df.mg$ATOH_f,df.mg$adult_specific)