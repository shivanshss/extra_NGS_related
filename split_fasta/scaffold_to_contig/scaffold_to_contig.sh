cat seqs.fa                \
    | seqkit fx2tab        \
    | cut -f 2             \
    | sed -r 's/n+/\n/gi'  \
    | cat -n               \
    | seqkit tab2fx        \
    | seqkit replace -p "(.+)" -r "Contig{nr}" > seq_contig.fa
