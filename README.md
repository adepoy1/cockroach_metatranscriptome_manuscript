# cockroach_metatranscriptome_manuscript

This GitHub repository for all the scripts and data used in a metatranscriptome analysis workflow for manuscript.

## References

First, we need to gather up reference genomes for a taxon of interest. The best place to get these are from NCBI FTP under RefSeq. Download the genbank files (.gbff) for each reference genomes. Make sure that the files that are downloaded from NCBI are GCF_ and no GCA_. Make sure that the files don't have dashes or periods in the file names. Change these if necessary. 
The following code can be be used to clean up the file names.

```
# This package found here <https://github.com/kblin/ncbi-genome-download> can be used to download a lot of genomes at once. It does nest the genomes in folders, so you'll have to move the .gbff files from those nested files for each genome. 
# install the package, following installation from <https://github.com/kblin/ncbi-genome-download>
pip install ncbi-genome-download

# use this command to download reference genomes. 
# in this case, all genomes from the genus Desulfovibrio
ncbi-genome-download bacteria --genera Desulfovibrio ## 60 references downloaded.

```

```
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

```
# Once you have all your reference genomes of interest, you need to concatenate these refs into one file that will be fed into the Anvio Pangenome pipeline. 
cat *.gbff > all_desulfovibrio_refs.gbff

```

```
# This can be used to download a lot of references at once. It does nest the genomes, so you'll have to move the .gbff files from those nested files for each genome. 
pip install ncbi-genome-download
ncbi-genome-download bacteria --genera Desulfovibrio ## 60 references downloaded.


# Once you have all your reference genomes of interest, you need to concatenate these refs into one file that will be fed into the Anvio Pangenome pipeline. 
cat *.gbff > all_desulfovibrio_refs.gbff

```


## Pangenome analysis using Anvio

## Generating Hit count files using custom perl scripts


