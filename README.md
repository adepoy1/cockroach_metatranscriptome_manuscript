# cockroach_metatranscriptome_manuscript

This GitHub repository for all the scripts and data used in a metatranscriptome analysis workflow for manuscript.

## References

First, we need to gather up reference genomes for a taxon of interest. The best place to get these are from NCBI FTP under RefSeq. Download the genbank files (.gbff) for each reference genomes. Make sure that the files that are downloaded from NCBI are GCF_ and no GCA_. Make sure that the files don't have dashes or periods in the file names. Change these if necessary. 

```
# This package found here <https://github.com/kblin/ncbi-genome-download/> can be used to download a lot of genomes at once. It does nest the genomes in folders, so you'll have to move the .gbff files from those nested files for each genome. 
# install the package, following installation instructions from <https://github.com/kblin/ncbi-genome-download/>

pip install ncbi-genome-download

# use this command to download reference genomes. 
# in this case, all genomes from the genus Desulfovibrio

ncbi-genome-download bacteria --genera Desulfovibrio ## 60 references downloaded.

```

```
# the following code can be be used to clean up the file names.
# clean up file names to be just the GCF_ assembly number. 
# this for loop cuts at second underscore. You may have to edit this code depending on how your reference files are named. 

for i in *_genomic.gbff
do
j=`echo $i | cut -d _ -f1-2`
mv $i $j.gbff
done

# for renaming SAG names
for i in *_contigs.gbff
do
j=`echo $i | cut -d _ -f1-3`
mv $i $j.gbff
done

for i in GCF_*
do
mv "$i" "$i.gbff"
done

```

## Pangenome analysis using Anvio

We have a bash script that will run through an entire Anvio pipeline. We relied heavily on Anvio documentation and tutorials to build this pipeline. 

```
bash scripts/bash_anvio_dev.sh -n 'desulfovibrio'

```

## Generating Hit count files using custom perl scripts

Once we have the gene clusters summary file from the Anvio pangenome analysis, we can generate hit count files with perl scripts. In the following perl scripts, the gene clusters summary file will be referred to as the RBH file. 

Files needed to run these scripts include:
- NCBI names file (scripts/names_wSAGs.dmp)
	- names_wSAGs file have SAGs added to the end of the file. Additional SAGs can be added in a similar way. This is only needed if the SAGs are not added to NCBI taxonomy yet. 
- NCBI nodes files (scripts/nodes.dmp)
	- the names and nodes files get updated fairly regularly, so make sure you get the most updated files. 
- all top hit files (alltophits_metatranscriptome)
	- alltophit files are generated using NCBI Blast and the metatranscriptome data. 
	- Files are in m8 blast format - if more information about generating these files is needed contact Dr. Ottesen (ottesen@uga.edu). We generated these files in alltophits/metatranscriptome using custom lab perl scripts. 
- gene cluster summary file from anvio (data/)


The following perl script pulls all the hits from your alltophit files associated with a taxon of interest. In this example, I am using Desulfovibrio.

```
# The directories for the nodes.dmp and names.dmp files are hard coded in the script. Therefore, these directories need to be changed to run the script. Change the directories on lines 6, 36, and 58
# This script can also take an outer limit with the -o flag. This would be used if you want to go beyond the genus level you are looking at. Use this option with caution as many hits will come up from the alltophits file that may not be in the gene cluster summary file from anvio.

# usage:
perl blast_m8_taxID_assign_and_get_subgroup_v2.pl [alltophits files] -s [output prefix] -t [NCBI taxID]

# example:
perl /Volumes/G-DRIVEThunderbolt3/scripts/perl_scripts_anvio_pangenome/blast_m8_taxID_assign_and_get_subgroup_v2.pl /Volumes/G-DRIVEThunderbolt3/alltophits_metaT_genbank/*_toGenbank -s Desulfovibrio_by_taxID -t 872 

```

Now that we have pulled all the hits for a taxon, we are going to get the top hits with priority given to those hits that are found in the gene clusters.

```
# Change line 13 to be the directory that you have your all top hits files for metaT or metaG (if not in present directory otherwise no change to that line will assume that it is in present directory)
# Results are in file [filename]_single_tophit.txt
# options -b bitcut -e eval -d tophit_database (name of tophit file)

# usage:
perl NR_tophit_from_hitcount_file_m8_Anvio_RBH.pl [alltophits with taxon hits pulled (from the previous perl script)] -d name of the all hits database -o [output prefix] -a [gene clusters summary]

# example:
perl /Volumes/G-DRIVEThunderbolt3/scripts/perl_scripts_anvio_pangenome/NR_tophit_from_hitcount_file_m8_Anvio_RBH.pl -o *_toGenbank_Desulfovibrio_by_taxID -d Desulfovibrio_tax_ID_all_hitcounts_database -a desulfovibrio-PAN_gene_clusters_summary.txt

```
Now that we have these single top hit files, we can gets counts for all of our gene clusters. 

```

# usage: 
perl blast_m8_KEGG_counter_v2_RBH_clusters_forAnvio_6.12.2023.pl [single_tophit.txt files (from the previous perl script)] -o [output prefix] -a [gene clusters summary]

# example:
perl blast_m8_KEGG_counter_v2_RBH_clusters_forAnvio.pl *_Desulfovibrio_by_taxID_single_tophit.txt -o Desulfovibrio_by_taxID -a desulfovibrio-PAN_gene_clusters_summary.txt

```

Now we can add in other annotations. The following tutorial will be adding in annotations from the [Transporter database]<https://www.tcdb.org/>

```
# make a diamond database from the transporter database file that was downloaded from their website. then map the gene cluster amino acid sequences to the TCDB diamond database.  desulfovibrio_GC_fasta.txt

#get all top hits 
perl /Volumes/G-DRIVE\ Thunderbolt\ 3/scripts/perl_scripts/blast_m8_extract_all_tophits.pl Dv_GC_tcdb_matches2_1.m8

perl /Volumes/G-DRIVE\ Thunderbolt\ 3/scripts/perl_scripts/blast_m8_extract_all_tophits.pl Ds_GC_tcdb_matches.m8

#get the single top hit
perl /Volumes/G-DRIVE\ Thunderbolt\ 3/scripts/perl_scripts/NR_tophit_from_hitcount_file_m8_all.pl Dv_GC_tcdb_matches2.m8_alltophits -d Dv_GC_TCDB
perl /Volumes/G-DRIVE\ Thunderbolt\ 3/scripts/perl_scripts/NR_tophit_from_hitcount_file_m8_all.pl Ds_GC_matches.m8_alltophits -d Ds_GC_TCDB


add in header to single top hit file
# copy 3rd column and make it a new 4th column (3rd column should be duplicated)
match	function	acc	desc	perc_id	start	end	start	end	len	len	e_val	bit

#get TCDB annotations from single top hit file 
perl /Volumes/G-DRIVE\ Thunderbolt\ 3/scripts/perl_scripts/TCDB_lookup_anvio_GC.pl tcdb_listSuperfamilies.txt tcdb_families.txt Dv_GC_tcdb_matches2.m8_alltophits_single_tophit.txt
perl /Volumes/G-DRIVE\ Thunderbolt\ 3/scripts/perl_scripts/TCDB_lookup_anvio_GC.pl tcdb_listSuperfamilies.txt tcdb_families.txt Ds_GC_matches.m8_alltophits_single_tophit.txt


# put the TCDB annotations into the GC Kegg hit_hitcount file
perl /Volumes/G-DRIVE\ Thunderbolt\ 3/scripts/perl_scripts/Cluster_file_annotator_anvio_GC.pl -h /Volumes/G-DRIVE\ Thunderbolt\ 3/genecluster_groups_analysis/anvio_pangenome/Desulfovibrio/anvio_dev/Desulfovibrio_by_taxID_Kegg_hit_hitcounts.txt -o Desulfovibrio_GC_hits_annotated.txt Dv_GC_tcdb_matches2.m8_alltophits_single_tophit.txt_annotated.txt

```
