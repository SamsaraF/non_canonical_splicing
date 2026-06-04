#!/usr/bin/env zsh
# Important: Must be runned under zsh interpreter

cd references

# gencode annotations def
typeset -A ref_urls=(
    [mouse_ref]="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37"
    [human_ref]="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49"
)

typeset -A mouse_ref=(
    [genome]="GRCm39.genome.fa"
    [transcripts]="gencode.vM37.transcripts.fa"
    [gtf]="gencode.vM37.annotation.gtf"
    [gff3]="gencode.vM37.annotation.gff3"
)

typeset -A human_ref=(
    [genome]="GRCh38.p14.genome.fa"
    [transcripts]="gencode.v49.transcripts.fa"
    [gtf]="gencode.v49.annotation.gtf"
    [gff3]="gencode.v49.annotation.gff3"
)


for ref_name base_url in "${(@kv)ref_urls}"; do
    if ! [ -d $ref_name ]; then
        mkdir $ref_name
    fi
    cd $ref_name

    # download references and annotations
    local -A files=("${(@kvP)ref_name}")
    for file in "${(@v)files}"; do
        if ! [ -f "$file.gz" ]; then
            wget "$base_url/$file.gz"
        fi
        if ! [ -f $file ]; then
            pigz -p 8 -d -c "$file.gz" > $file
        fi
    done

    # convert gtf to bed, using UCSC tools
    local bed_file="${files[gtf]:r}.bed"
    if ! [ -f "${files[gtf]:r}.bed" ]; then
        gtfToGenePred "${files[gtf]}" "${files[gtf]:r}.genePred"
        genePredToBed "${files[gtf]:r}.genePred" "${files[gtf]:r}.bed"
        # rm "${files[gtf]:r}.genePred"
    fi

    # STAR build index
    if ! [ -f "${files[genome]}.fai" ]; then
        samtools faidx "${files[genome]}"
    fi
    if ! [ -d "./star_index" ]; then
        STAR --runMode genomeGenerate --genomeDir ./star_index --runThreadN 32 --genomeFastaFiles "${files[genome]}" --sjdbGTFfile "${files[gtf]}"
    fi

    # Salmon transctipt index
    if ! [ -d "./salmon_index" ]; then
        grep "^>" "${files[genome]}" | cut -d " " -f 1 > decoys.txt
        sed -i -e 's/>//g' decoys.txt
        cat "${files[transcripts]}.gz" "${files[genome]}.gz"> gentrome.fa.gz
        salmon index -t gentrome.fa.gz -d decoys.txt -p 16 -i salmon_index --gencode
    fi

    cd ..
done
