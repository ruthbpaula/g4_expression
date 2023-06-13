#!/bin/bash

## get_g4_at_gene_locations.sh
## Written by: Ruth De Paula

# Formatting refFlat file for intersection (downloaded straight from UCSC genome browser)
grep -v NR_ refFlat.txt | cut -f 1,3,5,6 | gawk '{print $2 "\t" $3 "\t" $4 "\t" $1}' | grep -v "_" | grep -v chrM | sort -u > refFlat_to_intersect

# Formatting exon table for intersection (downloaded straight from Ensembl BioMart). Format: Chromosome/scaffold name, Exon region start (bp), Exon region end (bp), Gene name. Example:
# MT      577     647     MT-TF
grep -v CHR_ all_exon_coord_biomart.txt | grep -v -P "MT\t" | grep -v GL000 | grep -v KI270 | grep -v Chromosome > all_exon_biomart_clean
gawk '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4}' all_exon_biomart_clean > all_exon_biomart_ok

# Formatting g4 table for intersection
# all_g4_toIntersect: Generated from C++ Quad. Format: chromosome, g4_start, g4_end. Example:
# chr10   10005   10091
sed 's/chr0/chr/g' all_g4_toIntersect | sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g' > all_g4_toIntersect_ok

# Run bedtools intersect to find coordinates that overlap
bedtools intersect -a all_g4_toIntersect_ok -b refFlat_to_intersect -wo | sort -u > all_g4_transcripts
bedtools intersect -a all_g4_toIntersect_ok -b refFlat_to_intersect | sort -u > all_g4_transcripts_toIntersect

# Count number of genes with G4 at transcripts
cut -f 7 all_g4_transcripts | sort -u | wc -l  # 16,648 entries

# Save gene names
cut -f 7 all_g4_transcripts | sort -u > all_g4_transcripts_gene_names

# ---

# Run bedtools intersect to find coordinates that DO NOT OVERLAP transcripts (entries in A that do not overlap with B)
bedtools intersect -a refFlat_to_intersect -b all_g4_transcripts_toIntersect -v -wo | sort -u > all_genes_without_g4

# Count number of genes without G4 at transcripts
cut -f 4 all_genes_without_g4 | sort -u | wc -l  # 3,237 entries

# Save gene names
cut -f 4 all_genes_without_g4 | sort -u > all_genes_without_g4_gene_names

# ---

# Run bedtools intersect to find genes with G4 at exons
bedtools intersect -a all_g4_toIntersect_ok -b all_exon_biomart_ok -wo | sort -u > all_g4_exon

# Run bedtools intersect to find refFlat gene names
bedtools intersect -a all_g4_exon -b refFlat_to_intersect -wo | sort -u > all_g4_exon_genes

# Count number of genes with G4 at exons
cut -f 7 all_g4_exon_genes | sort -u | wc -l  # 10,966 entries

# Save gene names
cut -f 7 all_g4_exon_genes | sort -u > all_g4_exon_gene_names

# ---

# Formatting tables for intersection
bedtools intersect -a all_g4_transcripts_toIntersect -b refFlat_to_intersect -wo | sort -u > all_g4_transcripts_genes
cut -f 1,2,3,7 all_g4_transcripts_genes | sort -u > all_g4_transcripts_toIntersect_ok

# Run bedtools intersect to find coordinates that DO NOT OVERLAP (entries in A that do not overlap with B)
bedtools intersect -a all_g4_transcripts_toIntersect_ok -b all_exon_toIntersect_ok -v -wo | sort -u > all_g4_intron

# Count number of genes with G4 at introns
cut -f 4 all_g4_intron | sort -u | wc -l  # 15,103 entries

# Save gene names
cut -f 4 all_g4_intron | sort -u > all_g4_intron_gene_names

# ---

# Split original refFlat file (downloaded straight from UCSC genome browser) into + and -: done on spliced 3utr code
grep -v NR_ refFlat.txt | grep "+" > refFlat_+
grep -v NR_ refFlat.txt | grep -v "+" > refFlat_-

# For + strand, get 5'utr coord
gawk '{print $3 "\t" $5 "\t" $7 "\t" $1}' refFlat_+ > refFlat_unspliced_5utr_+  # chr, txStart, cdsStart, gene name

# For - strand, get 5'utr coord
gawk '{print $3 "\t" $8 "\t" $6 "\t" $1}' refFlat_- > refFlat_unspliced_5utr_-  # chr, cdsEnd, txEnd, gene name

# Merge + and - 5'utr coord
cat refFlat_unspliced_5utr_+ refFlat_unspliced_5utr_- | grep -v "_" | grep -v chrM | sort -u > refFlat_unspliced_5utr

# Run bedtools intersect to find common coordinates
bedtools intersect -a all_g4_exon -b refFlat_unspliced_5utr -wo | sort -u > all_g4_spliced_5utr

# Count number of genes with G4 at 5'utr
cut -f 7 all_g4_spliced_5utr | sort -u | wc -l  # 5,200 entries

# Save gene names
cut -f 7 all_g4_spliced_5utr | sort -u > all_g4_spliced_5utr_gene_names

# ---

# For + strand, get 3'utr coord
gawk '{print $3 "\t" $8 "\t" $6 "\t" $1}' refFlat_+ > refFlat_unspliced_3utr_+  # chr, cdsEnd, txEnd, gene name

# For - strand, get 3'utr coord
gawk '{print $3 "\t" $5 "\t" $7 "\t" $1}' refFlat_- > refFlat_unspliced_3utr_-  # chr, txStart, cdsStart, gene name

# Merge + and - 3'utr coord
cat refFlat_unspliced_3utr_+ refFlat_unspliced_3utr_- | grep -v "_" | grep -v chrM | sort -u > refFlat_unspliced_3utr

# Run bedtools intersect to find common coordinates
bedtools intersect -a all_g4_exon -b refFlat_unspliced_3utr -wo | sort -u > all_g4_spliced_3utr

# Count number of genes with G4 at 3'utr
cut -f 7 all_g4_spliced_3utr | sort -u | wc -l  # 4,757 entries

# Save gene names
cut -f 7 all_g4_spliced_3utr | sort -u > all_g4_spliced_3utr_gene_names

# ---

# For + strand, get TSS coord
gawk '{print $3 "\t" $5 "\t" $5 "\t" $1}' refFlat_+ > refFlat_TSS_+  # chr, txStart, txStart, gene name

# For - strand, get TSS coord
gawk '{print $3 "\t" $6 "\t" $6 "\t" $1}' refFlat_- > refFlat_TSS_-  # chr, txEnd, txEnd, gene name

# Merge + and - TSS coord
cat refFlat_TSS_+ refFlat_TSS_- | grep -v "_" | grep -v chrM | sort -u > refFlat_TSS

# Run bedtools intersect to find common coordinates
bedtools intersect -a all_g4_toIntersect_ok -b refFlat_TSS -wo | sort -u > all_g4_TSS

# Count number of genes with G4 at TSS
cut -f 7 all_g4_TSS | sort -u | wc -l  # 808 entries

# Save gene names
cut -f 7 all_g4_TSS | sort -u > all_g4_TSS_gene_names

# ---

# For + strand, get TSS -45 coord
gawk '{print $3 "\t" $5-45 "\t" $5-45 "\t" $1}' refFlat_+ > refFlat_TSS_-45_+  # chr, txStart-45, txStart-45, gene name

# For - strand, get TSS -45 coord (actually +45)
gawk '{print $3 "\t" $6+45 "\t" $6+45 "\t" $1}' refFlat_- > refFlat_TSS_-45_-  # chr, txEnd, txEnd, gene name

# Merge + and - TSS -45 coord
cat refFlat_TSS_-45_+ refFlat_TSS_-45_- | grep -v "_" | grep -v chrM | sort -u > refFlat_TSS_-45

# Run bedtools intersect to find common coordinates
bedtools intersect -a all_g4_toIntersect_ok -b refFlat_TSS_-45 -wo | sort -u > all_g4_TSS_-45

# Count number of genes with G4 at TSS -45
cut -f 7 all_g4_TSS_-45 | sort -u | wc -l  # 2,340 entries

# Save gene names
cut -f 7 all_g4_TSS_-45 | sort -u > all_g4_TSS_-45_gene_names
