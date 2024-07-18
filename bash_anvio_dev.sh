#!/bin/bash
# This script will run through the Anvio pipline to generate gene clusters.
#  > IMPORTANT - Anvio will not re-run the same analysis, so change the group name or 
#	 delete previously created files (genbank files are the only thing you need to run this).
# Use the group you are using in place of desulfovibrio
# Before running, make sure your Anvio conda environment is running
#   e.g. -->  conda activate anvio-7.1
# Prior to running also make sure to have your references downloaded.
# This script renames "hypothetical protein" to "HP" in order to keep them in gene clusters throughout anvio
# The anvi-pan-genome step has multiple search options, here: --sensitive is DIAMOND sensitive mode.
#   > no flag will be the default DIAMOND fast mode. 


# Setting up the development branch on the cluster: https://anvio.org/install/#5-follow-the-active-development-youre-a-wizard-arry
#  > 15G worked for setup; 2G was not enough memory and setup steps were killed
# ml Miniconda3/4.7.10 
# conda create -y --name anvio-dev python=3.7
# conda activate anvio-dev

# install packages

# link to Anvio codebase:
# cat <<EOF >${CONDA_PREFIX}/etc/conda/activate.d/anvio.sh
# #creating an activation script for the the conda environment for anvi'o
# #development branch so (1) Python knows where to find anvi'o libraries,
# #(2) the shell knows where to find anvi'o programs, and (3) every time
# #the environment is activated it synchronizes with the latest code from
# #active GitHub repository:
# export PYTHONPATH=\$PYTHONPATH:/scratch/hd55218/pangenomes/github_anvio/anvio
# export PATH=\$PATH:/scratch/hd55218/pangenomes/github_anvio/anvio/bin:/scratch/hd55218/pangenomes/github_anvio/anvio/sandbox
# echo -e "\033[1;34mUpdating from anvi'o GitHub \033[0;31m(press CTRL+C to cancel)\033[0m ..."
# cd /scratch/hd55218/pangenomes/github_anvio/anvio && git pull && cd -
# EOF



# for the SAGs, anvio doesn't like the dashes. I replaced them with underscores.
## also, anvio doesn't like periods in file name. change these to underscores

# Installing Anvio on M1 Chip Macs
# Install Miniforge
# Having trouble installing some of the required packages?
#	Try "conda update --all"
#	Try "conda config --set channel_priority flexible"
#	Make sure you have only created ONE environment for anvio
#		check this with "conda info --envs"

# Usage: bash bash_anvio.sh -n 'desulfovibrio'


# subgroup prefix parameter
while getopts "n:" arg; do
  case $arg in
    n) Group=$OPTARG;;
  esac
done

# 
files=$(find ./ -name 'GCA_*.gbff' -maxdepth 1)

# check for missing -n input
if [ -z $Group ] # if you left out the -n parameter
then
	echo 'Add a prefix for your subgroup with -n !'

# check for GCA_ files
elif [ -z $files ] # if the files variable is empty (bc there are no GCA_files)
then

# catting all reference genbanks into one
cat *.gbff > all_"$Group"_refs.gbff

# Replacing 'hypothetical protein' with 'HP in reference files
sed 's/hypothetical protein/HP/g' all_"$Group"_refs.gbff > all_"$Group"_refsPrepped.gbff
mv all_"$Group"_refs.gbff all_"$Group"_refs.gbk

# Parsing the genbank files
echo "Parsing $Group Genbank files..."
anvi-script-process-genbank -i all_"$Group"_refsPrepped.gbff \
 --output-gene-calls all_"$Group"_refs_gene_calls.tsv \
 --output-functions all_"$Group"_refs_functions.tsv \
 --output-fasta all_"$Group"_refs.fa \
 --annotation-source NCBI_PGAP
echo "$Group Genbank files parsed."


# Sometimes there's an error about one bad character in a sequence. This anvio script will ensure those are removed.
anvi-script-reformat-fasta all_"$Group"_refs.fa -o all_"$Group"_refs_reformat.fa --seq-type NT

# Generate a contigs database for all the references
echo "Generating contigs database for $Group..."
anvi-gen-contigs-database -f all_"$Group"_refs_reformat.fa \
-o "$Group"_contigs.db \
-n "$Group"_pangenome \
--external-gene-calls \
all_"$Group"_refs_gene_calls.tsv --ignore-internal-stop-codons
echo "Contigs database generated for $Group."

# Import NCBI PGAP function annotations that came from the Genbank files
echo "Importing $Group NCBI PGAP annotations..."
anvi-import-functions -c "$Group"_contigs.db \
-i all_"$Group"_refs_functions.tsv

# run HMMs (hidden markov models)
echo "Running HMMs on $Group..."
anvi-run-hmms -c "$Group"_contigs.db \
--num-threads 6

# Add in NCBI COG annotations (setup must happen once per anvio environment)
echo "Downloading COG db..."
anvi-setup-ncbi-cogs --reset
# download NCBI COG database
echo "Running NCBI COGs on $Group..."
anvi-run-ncbi-cogs -c "$Group"_contigs.db --num-threads 6

# Add in KEGG annotations (setup must happen once per anvio environment)
echo "Downloading KEGG db..."
anvi-setup-kegg-kofams --reset
# download KEGG database
echo "Running KEGG annotations $Group..."
anvi-run-kegg-kofams -c "$Group"_contigs.db --num-threads 6

# Add in CAZy annotations (setup must happen once per anvio environment)
echo "Downloading CAZy db..."
# download dbCAN database here (check for up-to-date versions)
anvi-setup-cazymes --reset
echo "Running CAZy annotations $Group..."
anvi-run-cazymes -c "$Group"_contigs.db  --num-threads 6
echo "Annotations complete for $Group."

# Now we need to make a collection file with each node and genome associated with that node
echo "Making collection file for $Group..."
for genome in $(ls *.gbff | cut -d . -f -2)
do
  grep "LOCUS" "$genome" | tr -s " " "\t" | cut -f2 > "$genome"_contigs.tmp
  for contig in $(cat "$genome"_contigs.tmp)
  do
    echo "$genome"
  done > "$genome"_name.tmp
  paste "$genome"_contigs.tmp "$genome"_name.tmp > "$genome"_for_cat.tmp
done


# delete the all_group_refs set (3 files)
echo "Cleaning up $Group temp files and prepping its collection..."
rm all_"$Group"_refsPrepped.gbff_contigs.tmp
rm all_"$Group"_refsPrepped.gbff_name.tmp
rm all_"$Group"_refsPrepped.gbff_for_cat.tmp

# concatenate to make the collection file
cat *_for_cat.tmp > "$Group"_collection.tsv
# remove the temp files
rm *.tmp

# find .gbff in collection file and replace with nothing
perl -pi -w -e 's/.gbff//g;' "$Group"_collection.tsv

# we will need a blank profile to store our collection
echo "Creating blank profile and importing $Group collection..."
anvi-profile -c "$Group"_contigs.db \
-o "$Group"_profile \
-S "$Group"_profile \
--blank-profile \
--min-contig-length 0 \
--skip-hierarchical-clustering

# now we can import our collection
anvi-import-collection "$Group"_collection.tsv \
-c "$Group"_contigs.db \
-p "$Group"_profile/PROFILE.db \
-C "$Group"_pangenome \
--contigs-mode
echo "$Group collection imported."

# making an internal genomes file so that anvio knows where to find the profile and contigs databases
echo "Making anvio internal genomes file for $Group..."
echo -e "name\tbin_id\tcollection_id\tprofile_db_path\tcontigs_db_path" > header.tmp
cut -f2 "$Group"_collection.tsv | uniq > name_and_bin_id.tmp
for i in $(cat name_and_bin_id.tmp); do echo ""$Group"_pangenome"; done > collection_id.tmp
for i in $(cat name_and_bin_id.tmp); do echo "$PWD/"$Group"_profile/PROFILE.db"; done > profile_db_path.tmp
for i in $(cat name_and_bin_id.tmp); do echo "$PWD/"$Group"_contigs.db"; done > contigs_db_path.tmp

paste name_and_bin_id.tmp name_and_bin_id.tmp collection_id.tmp profile_db_path.tmp contigs_db_path.tmp > body.tmp
cat header.tmp body.tmp > "$Group"_internal_genomes.tsv
rm *.tmp

# make a genomes storage database using the internal genomes file
echo "Generating $Group genomes storage db from internal genomes file..."
anvi-gen-genomes-storage -i "$Group"_internal_genomes.tsv \
-o "$Group"-GENOMES.db --gene-caller NCBI_PGAP

# run the pangenome analysis. This is where your amino acid sequences are aligned and MCL is used to make gene clusters.
# there are some parameters that can be changed like min-bit and MCL-inflation. I have kept the defaults.
echo "Running pangenome analysis on $Group..."
anvi-pan-genome -g "$Group"-GENOMES.db \
                -n "$Group"-PAN \
                --num-threads 6
echo "Finished pangenome analysis on $Group."

# now we can add a collection without having to open the pangenome circular figs and make bins
echo "Adding $Group collection..."
anvi-script-add-default-collection \
-p "$Group"-PAN/"$Group"-PAN-PAN.db \
-C "$Group" \
-b "$Group"

# get all the amino acid sequences for the gene clusters
echo "Generating $Group fasta file for all gene clusters..."
anvi-get-sequences-for-gene-clusters -g "$Group"-GENOMES.db \
-p "$Group"-PAN/"$Group"-PAN-PAN.db -o "$Group"_GC_fasta

# summarize the pangenome
# You need that "Group"_PAN_gene_clusters_summary.txt file from the SUMMARY file for the next part of the pipeline
echo "Summarizing $Group pangenome..."
anvi-summarize -p "$Group"-PAN/"$Group"-PAN-PAN.db \
-g "$Group"-GENOMES.db -C "$Group"
echo "Finished with $Group. Check for anvio errors before using parsing scripts from Liz to create taxon-specific orthologous gene cluster count tables."
echo "See Steps_after_Anvio.txt for more info."

# unzip
echo "Unzipping pangenome summary file..."
cd 

# PERL SCRIPTS
echo "Running perl scripts: pull NR hits from subgroup alltophits files..."
#perl /Volumes/G-DRIVE\ Thunderbolt\ 3/Helen/META/Diet_Metatranscriptomics/Anvio_pangenome_pipeline/NR_tophit_from_hitcount_file_m8_Anvio_RBH.pl \
#	*_toGenbank_"$Group" -d "$Group" -a /Volumes/G-DRIVE\ Thunderbolt\ 3/Helen/META/Diet_Metatranscriptomics/Anvio_pangenome_pipeline/refs/"$Group"/SUMMARY/"$Group"-PAN_gene_clusters_summary.txt
echo "Finished pulling NR hits from subgroup alltophits files."

echo "Running perl scripts: KEGG counter and hitcounts by GCs..."
#perl /Volumes/G-DRIVE\ Thunderbolt\ 3/Helen/META/Diet_Metatranscriptomics/Anvio_pangenome_pipeline/blast_m8_KEGG_counter_v2_RBH_clusters_forAnvio.pl \
#	./*_"$Group".txt_single_tophit.txt -o "$Group"_by_clusters \
#	-a /Volumes/G-DRIVE\ Thunderbolt\ 3/Helen/META/Diet_Metatranscriptomics/Anvio_pangenome_pipeline/refs/"$Group"/SUMMARY/"$Group"-PAN_gene_clusters_summary.txt
echo "Finished counting KEGG accs and GCs."

#echo "All done! Scroll back up to make sure Anvio steps did not generate any unexpected errors."
#echo "Annotate your hitcount file with TCDB annotations using Cluster_file_annotator_anvio_GCs.pl"


# if the files variable is NOT empty (bc there are GCA_files)
else
	echo 'We found GCA_ files in your directory! Anvio only likes files from RefSeq, and these start with GCF_ .'
fi


## Displaying the pangenome has to be done on personal computer. It has to open up a server and that is not allowed on the desktop.              
#anvi-display-pan -p Acholeplasma-PAN/Acholeplasma-PAN-PAN.db \
#                 -g Acholeplasma-GENOMES.db


# pull subgroup-specific transcripts with blast_m8_taxID_assign_and_get_subgroup_v2.pl
# pick non-redundant tophits from taxon-specific alltophits files
#echo "starting perl script for NR tophit from alltophits files for $Group..."
#
# tally hits by gene cluster, add higher levels of KEGG hierarchy
# The perl script for this command must be edited to properly find the kegg hierarchy file (ko00001.keg)
# > See blast_m8_KEGG_counter_v2_RBH_clusters_forAnvio.pl for more info
#echo "Starting perl script for tallying NR tophits and creating a gene cluster hitcounts table for $Group"
#perl /Volumes/G-DRIVE\ Thunderbolt\ 3/Helen/META/Diet_Metatranscriptomics/Anvio_pangenome_pipeline/blast_m8_KEGG_counter_v2_RBH_clusters_forAnvio.pl \
#	./*_"$Group".txt_single_tophit.txt -o "$Group"_by_clusters \
#	-a /Volumes/G-DRIVE\ Thunderbolt\ 3/Helen/META/Diet_Metatranscriptomics/Anvio_pangenome_pipeline/refs/"$Group"/SUMMARY/"$Group"-PAN_gene_clusters_summary.txt
#echo "All done! Scroll back up to make sure Anvio steps did not generate any unexpected errors."
