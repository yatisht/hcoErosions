# hcoErosions

This repository provides the code to make gene loss predictions (called *hcoErosions*) in a set of species (*species.txt*), given (i) a reference assembly and its annotated gene set (such as *hg38* and *ref_genes/hg38.transcriptCoding_clean.bed*, respectively), (ii) a set of pairwise alignment chains (that can be downloaded from http://hgdownload.cse.ucsc.edu or generated from raw genome assembly using [doBlastzChainNet.pl](https://github.com/ENCODE-DCC/kentUtils/blob/master/src/hg/utils/automation/doBlastzChainNet.pl)) and (iii) a phylogenetic tree for reference and query species (such as *species-tree.nh*), as described in our manuscript [1]. 

The code has been developed using **python2** and tested successfully on an Ubuntu 18.04 LTS Server but should work on other **Linux** distributions as well.

Below, we provide instructions to perform a sample analysis using human (hg38) as the reference assembly to query the mouse (mm10), rat (rn6), rabbit (cavPor3) and armadillo (dasNov3) genomes for high confidence gene losses. 

The sample query species are chosen in such a way so that loss predictions are made primarily for the mouse (mm10) and rat (rn6) genomes and that the computation finishes in <1 day on a 4-core CPU when Step 6 is parallelized. In particular, the sample analysis yields 30 out of the 68 reported gene losses for mouse and rat described in [1]. To yield all 68 predictions, the entire set of 58 placental mammals used in [1] should be provided as input via species.txt. This full analysis may require several days of compute time and may also require some alignment chains not available from http://hgdownload.cse.ucsc.edu to be generated using [doBlastzChainNet.pl](https://github.com/ENCODE-DCC/kentUtils/blob/master/src/hg/utils/automation/doBlastzChainNet.pl).

Note that other researchers can substitute or add to any or all of the species in the provided example to fit their own interests!

*Step 1:* Install prerequisites.
```
sudo apt install python
sudo apt install python-pip
pip install biopython treelib
git clone https://github.com/yatisht/hcoErosions
cd hcoErosions
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/twoBitToFa -P bin/
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/twoBitInfo -P bin/
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bedToBigBed -P bin/
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bedToExons -P bin/
chmod +x bin/*
```

*Step 2:* Download reference assembly (hg38.2bit) from UCSC Genome Browser and extract information.
```
if [ ! -f 2bit/hg38.2bit ]; then
    wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit -P 2bit/
    bin/twoBitInfo 2bit/hg38.2bit 2bit/hg38.chrom.sizes
fi
```

*Step 3:* Download 2bit assembly files for the list of query species in file <span style="font-family:san-serif">species.txt</span> from UCSC Genome Browser. 
```
for species in $(cat species.txt); do
    echo "Downloading 2bit file for ${species} ..."
    if [ ! -f 2bit/${species}.2bit ]; then
        wget http://hgdownload.cse.ucsc.edu/goldenpath/${species}/bigZips/${species}.2bit -P 2bit/
    fi
done
```

*Step 4:* Generate a bed file of potential assembly gaps from the 2bit files of each query species (listed in <span style="font-family:san-serif">species.txt</span>). These gap tracks are used later in step 8 to ensure that a query deletion in the chain file is not a result of an assembly gap. 
```
for species in $(cat species.txt); do
    if [ ! -f gaps/${species}.gap.bed ]; then
        bin/twoBitInfo -nBed  2bit/${species}.2bit stdout | sort -k1,1 -k2,2n > gaps/${species}.gap.bed
    fi
done
```

*Step 5:* Download pairwise alignment chain files between the human (hg38) and each query species in file <span style="font-family:san-serif">species.txt</span> from UCSC Genome Browser.
```
for species in $(cat species.txt); do
    echo "Downloading chain for hg38 and ${species} ..."
    if [ ! -f chains/hg38.${species}.all.chain ]; then
        Species=${species^}
        wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/vs${Species}/hg38.${species}.all.chain.gz -P chains/
        gunzip chains/hg38.${species}.all.chain.gz
    fi
done
```

*Step 6:* For each reference-query pair (where the query species are specified in species.txt), find an orthologous chain for each reference gene in file ref_genes/hg38.transcriptCoding_clean.bed. This step can take several hours (5-10 hrs per reference-query pair) but it is possible to trivially parallelize the following for loop on multiple cores for a large number of query species. 
```
for species in $(cat species.txt); do
    echo "Finding orthologous chain for each reference gene in ${species} alignments ..."
    if [ ! -f chain_ids/${species}.p ]; then
        python pick_ortho_chains.py -reference hg38 -query ${species} -genes ref_genes/hg38.transcriptCoding_clean.bed 
    fi
done
    
```

*Step 7:* For chains selected in step 6 above, we generate bigbed and anydb files to enable faster lookup of alignments in step 8. 
```
for species in $(cat species.txt); do
    echo "Generating bigbed and anydb files for hg38 and ${species} chain ..."
    if [ ! -f chains/hg38.${species}.all.chain.bb ]; then
        tmpfile=$(mktemp --suffix=.bed)
        grep chain chains/hg38.${species}.all.chain | awk '{print $3"\t"$6"\t"$7"\t"$13"\t0\t"$5}' | sort -k1,1 -k2,2n > ${tmpfile}
        bin/bedToBigBed ${tmpfile} 2bit/hg38.chrom.sizes chains/hg38.${species}.all.chain.bb
        rm ${tmpfile}
    fi

    if [ ! -f chains/hg38.${species}.all.chain.db ] | [ ! -f chains/hg38.${species}.all.chain.link.db ]; then
        python chain_to_anydb.py -reference hg38 -query ${species} 
    fi
done
    
```

*Step 8:* Use alignment chains to extract coding sequence of orthologously mapped genes for each species in file <span style="font-family:san-serif">species.txt</span>. 
```
for species in $(cat species.txt); do
    echo "Extracting coding sequence for ${species} ..."
    python extract_coding_seqs.py -reference hg38 -genes ref_genes/hg38.transcriptCoding_clean.bed -queries ${species} -outdir mappings/
done    
```

*Step 9:* For each extracted gene, compute Mahalanobis Distance (md, read [1]) to measure intactness of each gene wrt remaining mapped genes. A single dictionary ( <span style="font-family:san-serif">md/md_dict.p</span>) is created at the end of this step containing md values for all mapped genes for all species in  <span style="font-family:san-serif">species.txt</span>. 
```
python generate_md_values.py -reference hg38  -speciesList species.txt
```

*Step 10:* Apply a phylogenetic filer (read [1]) for each gene to make hcoErosions predictions. The phylogenetic tree is provided as input in Newick format (<span style="font-family:san-serif">species-tree.nh</span>). At the end of this step, a file  <span style="font-family:san-serif">hcoErosions.csv</span> is created that stores all hcoErosions predictions.
```
python phylo_filter.py -reference hg38 -speciesList species.txt -phyloTree species-tree.nh -geneTranscriptIds ref_genes/hg38.canonical_gene_transcript_ids.txt -outFile hcoErosions.csv
```


## Citation

[1] Yatish Turakhia*, Heidi Chen*, Amir Marcovitz*, and Gill Bejerano. "Loss of critical developmental and human disease-causing genes in 58 mammals." bioRxiv (2019): 819169 (* equal contributors).
