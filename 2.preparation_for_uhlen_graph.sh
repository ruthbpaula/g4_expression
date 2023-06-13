#!/bin/bash

## preparation_for_uhlen_graph.sh
## Written by: Ruth De Paula and Albino Bacolla

# 1. Convert genes to ensembl ids and grep with Uhlen table - https://biit.cs.ut.ee/gprofiler/convert:
grep -w -f exon_genes_list_ensembl ../uhlen_tableS18.txt > uhlen_exon_+.txt
grep -w -f intron_genes_list_ensembl ../uhlen_tableS18.txt > uhlen_intron_+.txt
grep -w -f spliced_3utr_genes_list_ensembl ../uhlen_tableS18.txt > uhlen_spliced_3utr_+.txt
grep -w -f spliced_5utr_genes_list_ensembl ../uhlen_tableS18.txt > uhlen_spliced_5utr_+.txt
grep -w -f TSS_-45_genes_list_ensembl ../uhlen_tableS18.txt > uhlen_TSS_-45_+.txt
grep -w -f TSS_genes_list_ensembl ../uhlen_tableS18.txt > uhlen_TSS_+.txt
grep -w -f transcripts_genes_list_ensembl ../uhlen_tableS18.txt > uhlen_transcripts_+.txt
grep -w -f genes_without_g4_genes_list_ensembl ../uhlen_tableS18.txt > uhlen_genes_without_g4_-.txt  ## will be the ultimate control

# 2. Run log2 of the mean of expression of each gene in Uhlen subtables:
outfile="uhlen_tableS18_mean"

gawk 'NR == 1 { next }
    { T = 0
      for(N = 3; N <= NF; N++) T+= $N;
      T/=(NF - 2)
      print $1 "\t" $2 "\t" log(T + 1)/log(2) "\t1.Exonic_G4" }' uhlen_exon_+.txt > ${outfile}

gawk '{ T = 0
        for(N = 3; N <= NF; N++) T+= $N;
        T/=(NF - 2)
        print $1 "\t" $2 "\t" log(T + 1)/log(2) "\t2.Intronic_G4" }' uhlen_intron_+.txt >> ${outfile}

gawk '{ T = 0
        for(N = 3; N <= NF; N++) T+= $N;
        T/=(NF - 2)
        print $1 "\t" $2 "\t" log(T + 1)/log(2) "\t3.5UTR_G4" }' uhlen_spliced_5utr_+.txt >> ${outfile}

gawk '{ T = 0
        for(N = 3; N <= NF; N++) T+= $N;
        T/=(NF - 2)
        print $1 "\t" $2 "\t" log(T + 1)/log(2) "\t4.3UTR_G4" }' uhlen_spliced_3utr_+.txt >> ${outfile}

gawk '{ T = 0
        for(N = 3; N <= NF; N++) T+= $N;
        T/=(NF - 2)
        print $1 "\t" $2 "\t" log(T + 1)/log(2) "\t5.TSS-45_G4" }' uhlen_TSS_-45_+.txt >> ${outfile}

gawk '{ T = 0
        for(N = 3; N <= NF; N++) T+= $N;
        T/=(NF - 2)
        print $1 "\t" $2 "\t" log(T + 1)/log(2) "\t6.TSS_G4" }' uhlen_TSS_+.txt >> ${outfile}

gawk '{ T = 0
        for(N = 3; N <= NF; N++) T+= $N;
        T/=(NF - 2)
        print $1 "\t" $2 "\t" log(T + 1)/log(2) "\t7.All_G4" }' uhlen_transcripts_+.txt >> ${outfile}
gawk '{ T = 0
        for(N = 3; N <= NF; N++) T+= $N;
        T/=(NF - 2)
        print $1 "\t" $2 "\t" log(T + 1)/log(2) "\t8.No_G4" }' uhlen_genes_without_g4_-.txt >> ${outfile}
		
# 3. Make "uhlen_table_toplot_ok" with columns 3 and 4 from "uhlen_tableS18_mean" and column 3 switched with column 4:
cut -f 3,4 uhlen_tableS18_mean | gawk '{print $2 "\t" $1}' > uhlen_table_toplot_ok