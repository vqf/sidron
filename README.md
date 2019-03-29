# Sidron
Scripts and example to understand Sidr&oacute;n. Instructions can be found in the publication [Genome Sequencing and Analysis Methods in Chronic Lymphocytic Leukemia](https://www.ncbi.nlm.nih.gov/pubmed/30350214). It should be noted that, after the installation of BLAT, we need to tell the scripts where the BLAT server is located with the environment variable BLATPORT. For instance,
```
export BLATPORT=9006
gfServer start localhost $BLATPORT <path_to_genome_fasta_file>
```
To avoid complications, the path `path_to_genome_fasta_file` should be absolute.
