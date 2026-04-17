#!/usr/bin/zsh

awk -v OFS='\t' 'BEGIN{smap["gene_id"]="symbol"} NR==FNR{smap[$1]=$2} NR>FNR{print smap[$1], $0}' \
    references/human_ref/ensid_to_symbol.tsv $1 | \
    awk -v g=$2 'BEGIN{ 
        n = split(g, t, ",");
        for(i=1; i<=n; i++){ gc[t[i]]="Y" }
    }
    NR==1||gc[$1]=="Y"' | batcat -l tsv
