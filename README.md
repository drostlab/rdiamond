# rdiamond

## An easy-to-use framework to perform massive [DIAMOND](http://www.diamondsearch.org/index.php) sequence searches with R

### Motivation 

Homepage: drostlab.github.io/rdiamond/

[DIAMOND](http://www.diamondsearch.org/index.php) is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data. The key features are:

- Pairwise alignment of proteins and translated DNA at 100x-20,000x speed of BLAST.
- Frameshift alignments for long read analysis.
- Low resource requirements and suitable for running on standard desktops or laptops.
- Various output formats, including BLAST pairwise, tabular and XML, as well as taxonomic classification.

The `rdiamond` package provides easy-to-use interface functions between R and the standalone command line tool [DIAMOND](http://www.diamondsearch.org/index.php). 

`rdiamond` is designed to enable a new level of data-driven genomics
research by providing the computational tools and data science standards needed
to perform reproducible genomics research at scale.

The exponentially growing number of available sequences in biological databases
revolutionizes the way modern life science research is conducted. Approximately
hundred thousand genomic sequences spanning diverse species from the tree of life
are currently publicly available and free to access. It is now possible to
access and retrieve this data automatically using the R package [biomartr](https://github.com/ropensci/biomartr)
and the next step is to harness this wealth of sequence diversity to explore
and detect novel patterns of evolvability, variation, and disease emergence using `rdiamond`.

The R package [biomartr](https://github.com/ropensci/biomartr)
__solves the problem of retrieving this vast amount of biological sequence data__ in a standardized and computationally reproducible way and the `rdiamond` package aims to __solve the problem of performing massive 
sequence searches on a tree-of-life scale__ in a standardized and computationally reproducible way. 

Both packages, `biomartr` and `rdiamond` are designed to complement
each other seamlessly to provide users with a tool set to automatically
retrieve thousands of biological sequences (thousands of genomes, proteomes, annotations, etc)
and to use these sequences to perform massive sequence searches with [DIAMOND](http://www.diamondsearch.org/index.php) to
extract novel patterns of similarity and divergence between large sets
of species.

### Install `rdiamond`

### For Linux Users:

Please install the `libpq-dev` library on you linux machine by typing into the terminal:

```
sudo apt-get install libpq-dev
```

### For all systems install `rdiamond` by typing

```r
# install.packages("devtools")
# install the current version of rdiamond on your system
devtools::install_github("drostlab/rdiamond", build_vignettes = TRUE, dependencies = TRUE)
```

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
              db_import  = FALSE)

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
Processing query block 1, reference block 1/1, shape 2/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 3/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0.001s]
Building seed filter...  [0s]
Searching alignments...  [0.003s]
Processing query block 1, reference block 1/1, shape 4/64.
Building reference seed array...  [0.001s]
Building query seed array...  [0s]
Computing hash join...  [0.001s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 5/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.003s]
Processing query block 1, reference block 1/1, shape 6/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0.001s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 7/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 8/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 9/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 10/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 11/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 12/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.004s]
Processing query block 1, reference block 1/1, shape 13/64.
Building reference seed array...  [0.001s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 14/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 15/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 16/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 17/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 18/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 19/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 20/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 21/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 22/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 23/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 24/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 25/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 26/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.001s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 27/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0.001s]
Building seed filter...  [0.001s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 28/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 29/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 30/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 31/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 32/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.001s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 33/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.001s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 34/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.001s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 35/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.001s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 36/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.002s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 37/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.001s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 38/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.001s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 39/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.001s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 40/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.001s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 41/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.001s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 42/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.002s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 43/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.003s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 44/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.003s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 45/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.004s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 46/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.006s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 47/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.007s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 48/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.008s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 49/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.007s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 50/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.007s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 51/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.01s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 52/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.013s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 53/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.014s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 54/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.027s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 55/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.027s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 56/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.029s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 57/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.03s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 58/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.03s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 59/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.03s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 60/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.031s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 61/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.033s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 62/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.033s]
Searching alignments...  [0.001s]
Processing query block 1, reference block 1/1, shape 63/64.
Building reference seed array...  [0s]
Building query seed array...  [0s]
Computing hash join...  [0s]
Building seed filter...  [0.033s]
Searching alignments...  [0.002s]
Processing query block 1, reference block 1/1, shape 64/64.
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


## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions provided in this package.

Furthermore, in case you find some bugs, need additional (more flexible) functionality of parts of this package, or want to contribute to this project please let me know:

https://github.com/drostlab/rdiamond/issues
