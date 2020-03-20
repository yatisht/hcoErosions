# hcoErosions

*Step 1:* Install prerequisites
```
sudo apt install python
sudo apt install python-pip
pip install biopython treelib
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/twoBitToFa -P bin/
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/twoBitInfo -P bin/
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bedToBigBed -P bin/
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bedToExons -P bin/
chmod +x bin/*
```

*Step 2:* Download reference assembly (hg38.2bit) from UCSC Genome Browser
```
if [ ! -f 2bit/hg38.2bit ]; then
    wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit -P 2bit/
    bin/twoBitInfo 2bit/h38.2bit 2bit/hg38.chrom.sizes
fi
```

*Step 3:* Download 2bit assembly files for the list of query species in file <span style="font-family:san-serif">species.list</span> from UCSC Genome Browser 
```
for species in $(cat species.list); do
    echo "Downloading 2bit file for ${species} ..."
    if [ ! -f 2bit/${species}.2bit ]; then
        wget http://hgdownload.cse.ucsc.edu/goldenpath/${species}/bigZips/${species}.2bit -P 2bit/
    fi
done
```

*Step 4:* Generate a bed file of potential assembly gaps from the 2bit files of each query species (listed in <span style="font-family:san-serif">species.list</span>). These gap tracks are used later in step 8 to ensure that a query deletion in the chain file is not a result of an assembly gap. 
```
for species in $(cat species.list); do
    if [ ! -f gaps/${species}.gap.bed ]; then
        tmpfile=$(mktemp --suffix=.fa)
        bin/twoBitToFa 2bit/${species}.2bit ${tmpfile}
        python generate_gap_bed.py ${tmpfile} > gaps/${species}.gap.bed
        rm ${tmpfile}
    fi
done
```

*Step 5:* Download pairwise alignment chain files between the human (hg38) and each query species in file <span style="font-family:san-serif">species.list</span> from UCSC Genome Browser 
```
for species in $(cat species.list); do
    echo "Dowloading chain for hg38 and ${species} ..."
    if [ ! -f chains/hg38.${species}.all.chain ]; then
        Species=${species^}
        wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/vs${Species}/hg38.${species}.all.chain.gz -P chains/
        gunzip chains/hg38.${species}.all.chain.gz
    fi
done
```

*Step 6:* For each reference-query pair (list of query species specified in  <span style="font-family:san-serif">species.list</span>), find an orthologous chain for each reference gene in file <span style="font-family:san-serif">ref_genes/hg38.transcriptCoding_clean.bed</span>. This step is can take several hours (5-10 per reference-query pair) but it is possible to trivially parallelize the for loop below across a large number of query species.  
```
for species in $(cat species.list); do
    echo "Finding orthologous chain for each reference gene in ${species} alignments ..."
    if [ ! -f chain_ids/${species}.p ]; then
        python pickChains.py -r hg38 -q ${species} -b ref_genes/hg38.transcriptCoding_clean.bed 
    fi
done
    
```

*Step 7:* For chains selected in step 6 above, we generate bigbed and anydb files to enable faster lookup of alignments in step 8. 
```
for species in $(cat species.list); do
    echo "Generating bigbed and anydb files for hg38 and ${species} chain ..."
    if [ ! -f chains/hg38.${species}.all.chain.bb ]; then
        tmpfile=$(mktemp --suffix=.bed)
        grep chain chains/hg38.${species}.all.chain | awk '{print $3"\t"$6"\t"$7"\t"$13"\t0\t"$5}' | sort -k1,1 -k2,2n > ${tmpfile}
        bin/bedToBigBed ${tmpfile} 2bit/hg38.chrom.sizes chains/hg38.${species}.all.chain.bb
        rm ${tmpfile}
    fi

    if [ ! -f chains/hg38.${species}.all.chain.db ] | [ ! -f chains/hg38.${species}.all.chain.link.db ]; then
        python chainToDb.py -r hg38 -q ${species} 
    fi
done
    
```

*Step 8:* Use alignment chains to extract coding sequence of orthologously mapped genes for each species in file <span style="font-family:san-serif">species.list</span>. 
```
for species in $(cat species.list); do
    python geneDiffs.py -reference hg38 -genes ref_genes/hg38.transcriptCoding_clean.bed -queries ${species} -outdir mappings/
done    
```

*Step 9:* For each extracted gene, compute Mahalanobis Distance (md, read [1]) to measure intactness of each gene wrt remaining mapped genes. A single dictionary ( <span style="font-family:san-serif">md/md_dict.p</span>) is created at the end of this step containing md values for all mapped genes for all species in  <span style="font-family:san-serif">species.list</span>. 
```
python generate_md_values.py -reference hg38  -speciesList species.list
```

*Step 10:* Apply a phylogenetic filer (read [1]) for each gene to make hcoErosions predictions. The phylogenetic tree is provided as input in Newick format (<span style="font-family:san-serif">species-tree.nh</span>). At the end of this step, a file  <span style="font-family:san-serif">hcoErosions.csv</span> is created that stores all hcoErosions predictions.
```
python phylo_filter.py -reference hg38 -speciesList species.list -phyloTree species-tree.nh -phyloTree ref_genes/hg38.canonical_gene_transcript_ids.txt -outFile hcoErosions.csv 
```


## Citation

[1] Yatish Turakhia, Heidi Chen, Amir Marcovitz, and Gill Bejerano. "Loss of critical developmental and human disease-causing genes in 58 mammals." bioRxiv (2019): 819169
