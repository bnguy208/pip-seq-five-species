#!/bin/bash

#SBATCH --job-name=demux_five
#SBATCH --output=demux_five_species_%j.out
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=campus-new
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bnguyen3@fredhutch.org

# Load modules
ml CellRanger
ml STAR
ml pipseeker

# Paths to FASTA files
fasta1="./genomes/Cat/Felis_catus.Felis_catus_9.0.dna.toplevel.fa"
fasta2="./genomes/Chicken/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly_combined.fa"
fasta3="./genomes/Green_Monkey/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa"
fasta4="./genomes/Rabbit/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa"
fasta5="./genomes/Rat/Rattus_norvegicus_mRatBN7.2.dna.primary_assembly_combined.fa"

# Paths to GTF files
gtf1="./genomes/Cat/Felis_catus.Felis_catus_9.0.111.gtf"
gtf2="./genomes/Chicken/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.111.gtf"
gtf3="./genomes/Green_Monkey/Chlorocebus_sabaeus.ChlSab1.1.111.gtf"
gtf4="./genomes/Rabbit/Oryctolagus_cuniculus.OryCun2.0.111.gtf"
gtf5="./genomes/Rat/Rattus_norvegicus.mRatBN7.2.111.gtf"

# Check if files exist
for file in $fasta1 $fasta2 $fasta3 $fasta4 $fasta5 $gtf1 $gtf2 $gtf3 $gtf4 $gtf5; do
    if [[ ! -f $file ]]; then
        echo "Error: $file not found!" >&2
        exit 1
    fi
done

# # Step 1: Concatenate FASTA files for different species
# # Add prefixes for each of the species
# sed 's/^>/>fc_/' $fasta1 > ./genomes/Cat/Felis_catus_prefixed.fa
# sed 's/^>/>gg_/' $fasta2 > ./genomes/Chicken/Gallus_gallus_prefixed.fa
# sed 's/^>/>cs_/' $fasta3 > ./genomes/Green_Monkey/Chlorocebus_sabaeus_prefixed.fa
# sed 's/^>/>oc_/' $fasta4 > ./genomes/Rabbit/Oryctolagus_cuniculus_prefixed.fa
# sed 's/^>/>rn_/' $fasta5 > ./genomes/Rat/Rattus_norvegicus_prefixed.fa

# # Paths to prefixed FASTA files
# pf_fasta1="./genomes/Cat/Felis_catus_prefixed.fa"
# pf_fasta2="./genomes/Chicken/Gallus_gallus_prefixed.fa"
# pf_fasta3="./genomes/Green_Monkey/Chlorocebus_sabaeus_prefixed.fa"
# pf_fasta4="./genomes/Rabbit/Oryctolagus_cuniculus_prefixed.fa"
# pf_fasta5="./genomes/Rat/Rattus_norvegicus_prefixed.fa"

# cat $pf_fasta1 $pf_fasta2 $pf_fasta3 $pf_fasta4 $pf_fasta5 > five_species_combined_genome.fa

# # Step 2: Filter the GTF files
# # Felis catus
# cellranger \
# 	mkgtf $gtf1 ./genomes/Cat/Felis_catus.Felis_catus_9.0.111.filtered.gtf \
# 		--attribute=gene_biotype:protein_coding \
# 		--attribute=gene_biotype:lncRNA \
# 		--attribute=gene_biotype:IG_C_gene \
# 		--attribute=gene_biotype:IG_D_gene \
# 		--attribute=gene_biotype:IG_J_gene \
# 		--attribute=gene_biotype:IG_LV_gene \
# 		--attribute=gene_biotype:IG_V_gene \
# 		--attribute=gene_biotype:TR_C_gene \
# 		--attribute=gene_biotype:TR_D_gene \
# 		--attribute=gene_biotype:TR_J_gene \
# 		--attribute=gene_biotype:TR_V_gene \
# 		--attribute=gene_biotype:IG_V_pseudogene \
# 		--attribute=gene_biotype:IG_J_pseudogene \
# 		--attribute=gene_biotype:IG_C_pseudogene \
# 		--attribute=gene_biotype:TR_V_pseudogene \
# 		--attribute=gene_biotype:TR_J_pseudogene

# # Gallus gallus
# cellranger \
# 	mkgtf $gtf2 ./genomes/Chicken/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.111.filtered.gtf \
# 		--attribute=gene_biotype:protein_coding \
# 		--attribute=gene_biotype:lncRNA \
# 		--attribute=gene_biotype:IG_C_gene \
# 		--attribute=gene_biotype:IG_D_gene \
# 		--attribute=gene_biotype:IG_J_gene \
# 		--attribute=gene_biotype:IG_LV_gene \
# 		--attribute=gene_biotype:IG_V_gene \
# 		--attribute=gene_biotype:TR_C_gene \
# 		--attribute=gene_biotype:TR_D_gene \
# 		--attribute=gene_biotype:TR_J_gene \
# 		--attribute=gene_biotype:TR_V_gene \
# 		--attribute=gene_biotype:IG_V_pseudogene \
# 		--attribute=gene_biotype:IG_J_pseudogene \
# 		--attribute=gene_biotype:IG_C_pseudogene \
# 		--attribute=gene_biotype:TR_V_pseudogene \
# 		--attribute=gene_biotype:TR_J_pseudogene

# # Chlorocebus sabaeus
# cellranger \
# 	mkgtf $gtf3 ./genomes/Green_Monkey/Chlorocebus_sabaeus.ChlSab1.1.111.filtered.gtf \
# 		--attribute=gene_biotype:protein_coding \
# 		--attribute=gene_biotype:lncRNA \
# 		--attribute=gene_biotype:IG_C_gene \
# 		--attribute=gene_biotype:IG_D_gene \
# 		--attribute=gene_biotype:IG_J_gene \
# 		--attribute=gene_biotype:IG_LV_gene \
# 		--attribute=gene_biotype:IG_V_gene \
# 		--attribute=gene_biotype:TR_C_gene \
# 		--attribute=gene_biotype:TR_D_gene \
# 		--attribute=gene_biotype:TR_J_gene \
# 		--attribute=gene_biotype:TR_V_gene \
# 		--attribute=gene_biotype:IG_V_pseudogene \
# 		--attribute=gene_biotype:IG_J_pseudogene \
# 		--attribute=gene_biotype:IG_C_pseudogene \
# 		--attribute=gene_biotype:TR_V_pseudogene \
# 		--attribute=gene_biotype:TR_J_pseudogene

# # Oryctolagus cuniculus
# cellranger \
# 	mkgtf $gtf4 ./genomes/Rabbit/Oryctolagus_cuniculus.OryCun2.0.111.filtered.gtf \
# 		--attribute=gene_biotype:protein_coding \
# 		--attribute=gene_biotype:lncRNA \
# 		--attribute=gene_biotype:IG_C_gene \
# 		--attribute=gene_biotype:IG_D_gene \
# 		--attribute=gene_biotype:IG_J_gene \
# 		--attribute=gene_biotype:IG_LV_gene \
# 		--attribute=gene_biotype:IG_V_gene \
# 		--attribute=gene_biotype:TR_C_gene \
# 		--attribute=gene_biotype:TR_D_gene \
# 		--attribute=gene_biotype:TR_J_gene \
# 		--attribute=gene_biotype:TR_V_gene \
# 		--attribute=gene_biotype:IG_V_pseudogene \
# 		--attribute=gene_biotype:IG_J_pseudogene \
# 		--attribute=gene_biotype:IG_C_pseudogene \
# 		--attribute=gene_biotype:TR_V_pseudogene \
# 		--attribute=gene_biotype:TR_J_pseudogene

# # Rattus norvegicus
# cellranger \
# 	mkgtf $gtf5 ./genomes/Rat/Rattus_norvegicus.mRatBN7.2.111.filtered.gtf \
# 		--attribute=gene_biotype:protein_coding \
# 		--attribute=gene_biotype:lncRNA \
# 		--attribute=gene_biotype:IG_C_gene \
# 		--attribute=gene_biotype:IG_D_gene \
# 		--attribute=gene_biotype:IG_J_gene \
# 		--attribute=gene_biotype:IG_LV_gene \
# 		--attribute=gene_biotype:IG_V_gene \
# 		--attribute=gene_biotype:TR_C_gene \
# 		--attribute=gene_biotype:TR_D_gene \
# 		--attribute=gene_biotype:TR_J_gene \
# 		--attribute=gene_biotype:TR_V_gene \
# 		--attribute=gene_biotype:IG_V_pseudogene \
# 		--attribute=gene_biotype:IG_J_pseudogene \
# 		--attribute=gene_biotype:IG_C_pseudogene \
# 		--attribute=gene_biotype:TR_V_pseudogene \
# 		--attribute=gene_biotype:TR_J_pseudogene 

# # Step 3: Concatenate filtered GTF files for different species
# # Paths to filtered GTF files
# filtered_gtf1="./genomes/Cat/Felis_catus.Felis_catus_9.0.111.filtered.gtf"
# filtered_gtf2="./genomes/Chicken/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.111.filtered.gtf"
# filtered_gtf3="./genomes/Green_Monkey/Chlorocebus_sabaeus.ChlSab1.1.111.filtered.gtf"
# filtered_gtf4="./genomes/Rabbit/Oryctolagus_cuniculus.OryCun2.0.111.filtered.gtf"
# filtered_gtf5="./genomes/Rat/Rattus_norvegicus.mRatBN7.2.111.filtered.gtf"

# # Add prefixes for each of the species
# sed '/^#!/!s/^/fc_/; s/gene_id "\([^"]*\)"/gene_id "fc_\1"/; s/gene_name "\([^"]*\)"/gene_name "fc_\1"/' $filtered_gtf1 > ./genomes/Cat/Felis_catus_prefixed.gtf
# sed '/^#!/!s/^/gg_/; s/gene_id "\([^"]*\)"/gene_id "gg_\1"/; s/gene_name "\([^"]*\)"/gene_name "gg_\1"/' $filtered_gtf2 > ./genomes/Chicken/Gallus_gallus_prefixed.gtf
# sed '/^#!/!s/^/cs_/; s/gene_id "\([^"]*\)"/gene_id "cs_\1"/; s/gene_name "\([^"]*\)"/gene_name "cs_\1"/' $filtered_gtf3 > ./genomes/Green_Monkey/Chlorocebus_sabaeus_prefixed.gtf
# sed '/^#!/!s/^/oc_/; s/gene_id "\([^"]*\)"/gene_id "oc_\1"/; s/gene_name "\([^"]*\)"/gene_name "oc_\1"/' $filtered_gtf4 > ./genomes/Rabbit/Oryctolagus_cuniculus_prefixed.gtf
# sed '/^#!/!s/^/rn_/; s/gene_id "\([^"]*\)"/gene_id "rn_\1"/; s/gene_name "\([^"]*\)"/gene_name "rn_\1"/' $filtered_gtf5 > ./genomes/Rat/Rattus_norvegicus_prefixed.gtf

# # # From 10/21/24:
# # sed '/^#!/!s/^/fc_/; s/gene_name "\([^"]*\)"/gene_name "fc_\1"/' $filtered_gtf1 > ./genomes/Cat/Felis_catus_prefixed.gtf
# # sed '/^#!/!s/^/gg_/; s/gene_name "\([^"]*\)"/gene_name "gg_\1"/' $filtered_gtf2 > ./genomes/Chicken/Gallus_gallus_prefixed.gtf
# # sed '/^#!/!s/^/cs_/; s/gene_name "\([^"]*\)"/gene_name "cs_\1"/' $filtered_gtf3 > ./genomes/Green_Monkey/Chlorocebus_sabaeus_prefixed.gtf
# # sed '/^#!/!s/^/oc_/; s/gene_name "\([^"]*\)"/gene_name "oc_\1"/' $filtered_gtf4 > ./genomes/Rabbit/Oryctolagus_cuniculus_prefixed.gtf
# # sed '/^#!/!s/^/rn_/; s/gene_name "\([^"]*\)"/gene_name "rn_\1"/' $filtered_gtf5 > ./genomes/Rat/Rattus_norvegicus_prefixed.gtf

# # Remove metadata header from subsequent files
# sed '/^#!/d' ./genomes/Chicken/Gallus_gallus_prefixed.gtf > ./genomes/Chicken/Gallus_gallus_prefixed_noheader.gtf
# sed '/^#!/d' ./genomes/Green_Monkey/Chlorocebus_sabaeus_prefixed.gtf > ./genomes/Green_Monkey/Chlorocebus_sabaeus_prefixed_noheader.gtf
# sed '/^#!/d' ./genomes/Rabbit/Oryctolagus_cuniculus_prefixed.gtf > ./genomes/Rabbit/Oryctolagus_cuniculus_prefixed_noheader.gtf
# sed '/^#!/d' ./genomes/Rat/Rattus_norvegicus_prefixed.gtf > ./genomes/Rat/Rattus_norvegicus_prefixed_noheader.gtf

# # No header GTF files
# no_head_gtf2="./genomes/Chicken/Gallus_gallus_prefixed_noheader.gtf"
# no_head_gtf3="./genomes/Green_Monkey/Chlorocebus_sabaeus_prefixed_noheader.gtf"
# no_head_gtf4="./genomes/Rabbit/Oryctolagus_cuniculus_prefixed_noheader.gtf"
# no_head_gtf5="./genomes/Rat/Rattus_norvegicus_prefixed_noheader.gtf"

# cat ./genomes/Cat/Felis_catus_prefixed.gtf $no_head_gtf2 $no_head_gtf3 $no_head_gtf4 $no_head_gtf5 > five_species_combined_genes.gtf

# # Step 4: Create a combined genome reference using Cell Ranger
# STAR --runThreadN 40 \
#      --runMode genomeGenerate \
#      --genomeDir five_species_reference_genome \
#      --genomeFastaFiles five_species_combined_genome.fa \
#      --genomeSAsparseD 3 \
#      --sjdbGTFfile five_species_combined_genes.gtf \
#      --sjdbOverhang 69 \
# 	 --limitGenomeGenerateRAM 38000000000 \
# 	 --limitSjdbInsertNsj 1062830

# Step 5: Run PIPseeker for FASTQ processing, alignment, and transcript counting
pipseeker full \
	--fastq sequence_fastq/. \
	--star-index-path five_species_reference_genome \
	--hto-fastq hto_fastq/. \
	--hto-tags tags.csv \
	--hto-position 0 \
	--chemistry v4 \
    --output-path pipseeker_output \
    --skip-version-check \
    --resume-last-run