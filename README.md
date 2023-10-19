## binf_align_denovo

### Assembling the given reads as fasta file into contigs for denovo alignment

`Functions`:

1. Main function:
- assemble_reads: Function to get contigs from a file of sequence reads 

2. Helper functions
- _get_reads: function to read sequences from the given file
- _get_contigs: function to merge reads with specific overlap
- _check_contig: function to check if reads can be mergerd with same k/overlap again
- _get_con:function to get final contig list
