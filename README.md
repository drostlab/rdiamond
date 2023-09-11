![Visitors](https://api.visitorbadge.io/api/visitors?path=https%3A%2F%2Fgithub.com%2Fdrostlab%2Frdiamond&label=VISITORS&countColor=%23263759&style=flat)


# rdiamond

## Seamless Integration of [DIAMOND2](https://github.com/bbuchfink/diamond)  Sequence Searches in R

### Motivation 

We are excited to introduce [DIAMOND2](https://www.nature.com/articles/s41592-021-01101-x), a cutting-edge pairwise protein aligner tailored to meet the extensive demands of the [Earth BioGenome Project](https://www.earthbiogenome.org/) and other expansive genomics initiatives. [DIAMOND2](https://github.com/bbuchfink/diamond) is a groundbreaking software solution designed to accelerate `BLAST` searches by an factor of up to 10,000x. To offer researchers even more flexibility and integration, we provide `rdiamond`, a dedicated interface package that allows programmatic handling of [DIAMOND2](https://github.com/bbuchfink/diamond) sequence searches directly through R. 

The `rdiamond` package offers streamlined interface functions, enabling users to seamlessly run [DIAMOND2](https://github.com/bbuchfink/diamond) directly within R. Notably, it's designed to handle vast outputs, processing terabytes of [DIAMOND2](https://github.com/bbuchfink/diamond) hit files directly from the disk on a local machine, bypassing memory limitations.

Furthermore, when paired with the [biomartr](https://github.com/ropensci/biomartr) R package, users have the convenience of automatically fetching large-scale genomic data and subsequently searching through it using rdiamond."

This version emphasizes the utility and integration capabilities of the rdiamond package while maintaining clarity.

### Install `rdiamond`

### For Linux Users:

Please install the `libpq-dev` library on you linux machine by typing into the terminal:

```
sudo apt-get install libpq-dev
```

### For all systems install `rdiamond` by typing

```r
# install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

# install Biostrings -> see here for different Biostrings verions:
# http://bioconductor.org/about/release-announcements/
BiocManager::install(c("Biostrings"))

# install.packages("devtools")
# install the current version of rdiamond on your system
devtools::install_github("drostlab/rdiamond", build_vignettes = TRUE, dependencies = TRUE)
```

### Citation

This R package is not formally published yet, but please cite the following paper when using this software for your research:

> Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). [doi:10.1038/s41592-021-01101-x](https://www.nature.com/articles/s41592-021-01101-x)



### Quick start
 
```r
# run diamond assuming that the diamond executable is available
# via the system path ('diamond_exec_path = NULL') and using
# sensitivity_mode = "ultra-sensitive"
diamond_example <- rdiamond::diamond_protein_to_protein(
              query   = system.file('seqs/qry_aa.fa', package = 'rdiamond'),
              subject = system.file('seqs/sbj_aa.fa', package = 'rdiamond'),
              sensitivity_mode = "ultra-sensitive",
              output_path = tempdir(),
              use_arrow_duckdb_connection  = FALSE)

# look at DIAMOND results
diamond_example
```

```
Running diamond with 'diamond blastp  with query: /library/rdiamond/seqs/qry_aa.fa and subject: /library/rdiamond/seqs/sbj_aa.fa using 2 core(s) ...


Alignment Sensitivity Mode:  --ultra-sensitive and max-target-seqs: 500


Masking of low-complexity regions: TANTAN


diamond v2.0.4.142 (C) Max Planck Society for the Advancement of Science
Documentation, support and updates available at http://www.diamondsearch.org

#CPU threads: 4
Scoring parameters: (Matrix=BLOSUM62 Lambda=0.267 K=0.041 Penalties=11/1)
Database input file: /library/rdiamond/seqs/sbj_aa.fa
Opening the database file...  [0.001s]
Loading sequences...  [0s]
Masking sequences...  [0.003s]
Writing sequences...  [0s]
Hashing sequences...  [0s]
Loading sequences...  [0s]
Writing trailer...  [0s]
Closing the input file...  [0s]
Closing the database file...  [0s]
Database hash = e7cd8e84df51b22dd27f3c01d5765fe1
Processed 20 sequences, 9444 letters.
Total time = 0.006s
diamond v2.0.4.142 (C) Max Planck Society for the Advancement of Science
Documentation, support and updates available at http://www.diamondsearch.org

#CPU threads: 2
Scoring parameters: (Matrix=BLOSUM62 Lambda=0.267 K=0.041 Penalties=11/1)
Temporary directory: /var/folders/yn/mgwl8_b56hz4v2c2vlfxj07w00076j/T//RtmpUCSsbL
Opening the database...  [0.002s]
#Target sequences to report alignments for: 500
Reference = /var/folders/yn/mgwl8_b56hz4v2c2vlfxj07w00076j/T//RtmpUCSsbL/sbj_aa.fa.dmnd
Sequences = 20
Letters = 9444
Block size = 400000000
Opening the input file...  [0s]
Opening the output file...  [0s]
Loading query sequences...  [0s]
Masking queries...  [0.002s]
Building query seed set... Algorithm: Double-indexed
 [0s]
Building query histograms...  [0.008s]
Allocating buffers...  [0s]
Loading reference sequences...  [0s]
Masking reference...  [0.001s]
Initializing temporary storage...  [0.013s]
Building reference histograms...  [0.007s]
Allocating buffers...  [0s]
Processing query block 1, reference block 1/1, shape 1/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0.001s]
Building seed filter...  [0.001s]
Searching alignments...  [0.002s]

... 


Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.033s]
Searching alignments...  [0.001s]
Deallocating buffers...  [0s]
Clearing query masking...  [0s]
Computing alignments...  [0.065s]
Deallocating reference...  [0s]
Loading reference sequences...  [0s]
Deallocating buffers...  [0s]
Deallocating queries...  [0s]
Loading query sequences...  [0s]
Closing the input file...  [0s]
Closing the output file...  [0s]
Closing the database file...  [0s]
Deallocating taxonomy...  [0s]
Total time = 0.797s
Reported 20 pairwise alignments, 20 HSPs.
20 queries aligned.


DIAMOND search finished successfully! 


The DIAMOND output file was imported into the running R session. The DIAMOND output file has been stored at: /var/folders/RtmpUCSsbL/qry_aa_sbj_aa_blastp_eval_0.001.blast_tbl
```

```
A tibble: 20 x 20
   query_id subject_id perc_identity num_ident_match…
   <chr>    <chr>              <dbl>            <int>
 1 333554|… AT1G01010…          73.2              347
 2 470181|… AT1G01020…          91.1              224
 3 470180|… AT1G01030…          93.3              335
 4 333551|… AT1G01040…          93.4             1840
 5 909874|… AT1G01050…         100                213
 6 470177|… AT1G01060…          87.5              567
 7 918864|… AT1G01070…          92.6              339
 8 909871|… AT1G01080…          89.3              268
 9 470171|… AT1G01090…          96.8              420
10 333544|… AT1G01110…          87.7              463
11 918858|… AT1G01120…          99.2              525
12 470161|… AT1G01140…          98.5              446
13 918855|… AT1G01150…          72.6              207
14 918854|… AT1G01160…          78.8              141
15 311317|… AT1G01170…          85.6               83
16 909860|… AT1G01180…          92.6              287
17 311315|… AT1G01190…          94.2              502
18 470156|… AT1G01200…          95.8              228
19 311313|… AT1G01210…          95.3              102
20 470155|… AT1G01220…          96.5             1019
# … with 16 more variables: alig_length <int>,
#   mismatches <int>, gap_openings <int>, n_gaps <int>,
#   pos_match <int>, ppos <dbl>, q_start <int>,
#   q_end <int>, q_len <int>, qcovhsp <dbl>, s_start <int>,
#   s_end <int>, s_len <int>, evalue <dbl>,
#   bit_score <dbl>, score_raw <dbl>
```

or 

```r
dplyr::glimpse(diamond_example)
```

```
Rows: 20
Columns: 20
$ query_id          <chr> "333554|PACid:16033839", "470181|PACid:16064328", "470180|PAC…
$ subject_id        <chr> "AT1G01010.1", "AT1G01020.1", "AT1G01030.1", "AT1G01040.1", "…
$ perc_identity     <dbl> 73.2, 91.1, 93.3, 93.4, 100.0, 87.5, 92.6, 89.3, 96.8, 87.7, …
$ num_ident_matches <int> 347, 224, 335, 1840, 213, 567, 339, 268, 420, 463, 525, 446, …
$ alig_length       <int> 474, 246, 359, 1969, 213, 648, 366, 300, 434, 528, 529, 453, …
$ mismatches        <int> 75, 22, 20, 58, 0, 71, 23, 25, 8, 65, 4, 6, 68, 30, 0, 20, 30…
$ gap_openings      <int> 8, 0, 2, 7, 0, 5, 2, 2, 3, 0, 0, 1, 3, 2, 1, 1, 1, 0, 0, 0
$ n_gaps            <int> 52, 0, 4, 71, 0, 10, 4, 7, 6, 0, 0, 1, 10, 8, 14, 3, 1, 0, 0,…
$ pos_match         <int> 369, 231, 338, 1870, 213, 587, 342, 275, 425, 475, 527, 448, …
$ ppos              <dbl> 77.8, 93.9, 94.2, 95.0, 100.0, 90.6, 93.4, 91.7, 97.9, 90.0, …
$ q_start           <int> 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 16, 2, 4, 1, 6, 1, 1, 1
$ q_end             <int> 466, 246, 355, 1963, 213, 640, 362, 299, 433, 528, 529, 453, …
$ q_len             <int> 466, 246, 355, 1963, 213, 640, 362, 299, 433, 528, 529, 453, …
$ qcovhsp           <dbl> 100.0, 100.0, 100.0, 99.9, 100.0, 100.0, 100.0, 100.0, 100.0,…
$ s_start           <int> 1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 5, 4, 2, 1, 1, 1, 1, 1
$ s_end             <int> 430, 246, 359, 1910, 213, 646, 366, 294, 429, 528, 529, 452, …
$ s_len             <int> 430, 246, 359, 1910, 213, 646, 366, 294, 429, 528, 529, 452, …
$ evalue            <dbl> 1.83e-212, 2.79e-157, 4.66e-209, 0.00e+00, 4.51e-158, 0.00e+0…
$ bit_score         <dbl> 584, 428, 568, 3541, 427, 1041, 613, 499, 816, 841, 1029, 866…
$ score_raw         <dbl> 1506, 1100, 1464, 9181, 1098, 2691, 1581, 1284, 2108, 2172, 2…
```

## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions provided in this package.

Furthermore, in case you find some bugs, need additional (more flexible) functionality of parts of this package, or want to contribute to this project please let me know:

https://github.com/drostlab/rdiamond/issues
