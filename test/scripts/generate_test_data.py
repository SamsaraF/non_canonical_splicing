import sys
import os
import numpy as np


NP_RANDOM_SEED = 42
np.random.seed(NP_RANDOM_SEED)

# Splicing sites definitions:
## U2: +GT/AG -CT/AC; U2 weak: +GC/AG -CT/GC; U12: +AT/AC -GT/AT;
AC_MOTIFS_MAP = {
    '+' : {'U2': 'AG', 'U2w': 'AG', 'U12': 'AC'},
    '-' : {'U2': 'AC', 'U2w': 'GC', 'U12': 'AT'}
}
DN_MOTIFS_MAP = {
    '+' : {'U2': 'GT', 'U2w': 'GC', 'U12': 'AT'},
    '-' : {'U2': 'CT', 'U2w': 'CT', 'U12': 'GT'}
}
all_motifs_set = set(i+j for i in 'ATCG' for j in 'ATCG')
# if not belong to any strand U2/U2w/U12, assign non-canonical motif
NC_AC_MOTIFS_SET = all_motifs_set - set(AC_MOTIFS_MAP['+'].values()) - set(AC_MOTIFS_MAP['-'].values())
NC_DN_MOTIFS_SET = all_motifs_set - set(DN_MOTIFS_MAP['+'].values()) - set(DN_MOTIFS_MAP['-'].values())

# Splicing site usage ratio defs:
## Non-U2 ratio are extradicted 10 times for generating enough data, to test these rare events.
U2_RATIO = 0.7
U2W_RATIO = 0.1
U12_RATIO = 0.1
NC_RATIO = 0.1


def add_splice_by_deletion():
    pass



def create_test_genome(length=10000):
    """Create a random genome sequence."""
    bases = ['A', 'T', 'C', 'G']
    return ''.join(np.random.choice(bases, length))


def add_splice_site(exons, genome, strand='+'):
    for e in range(1, len(exons) -1):
        # randomly choose splice type based on defined ratios
        splice_type = np.random.choice(
            ['U2', 'U2w', 'U12', 'NC'],
            p=[U2_RATIO, U2W_RATIO, U12_RATIO, NC_RATIO]
        )
        # determine motif sequence
        if splice_type == 'NC':
            ac_motif = np.random.choice(list(NC_AC_MOTIFS_SET))
            dn_motif = np.random.choice(list(NC_DN_MOTIFS_SET))
        else:
            ac_motif = AC_MOTIFS_MAP[strand][splice_type]
            dn_motif = DN_MOTIFS_MAP[strand][splice_type]

        genome = genome[:exons[e-1][1]] + ac_motif + genome[exons[e-1][1] + 2:]
        genome = genome = genome[:exons[e][0]-2] + dn_motif + genome[exons[e][0]:]


class Annotation:
    def __init__(self):
        self.genes = {}

    def add_gene(self, gene_id, gene_obj):
        if gene_id not in self.genes:
            self.genes[gene_id] = gene_obj
        else:
            raise ValueError(f"Gene {gene_id} already exists.")

class Gene:
    def __init__(self, gene_id, chrom, strand):
        self.gene_id = gene_id
        self.chrom = chrom
        self.strand = strand
        self.transcripts = {}

    def add_transcript(self, transcript_id, transcript_obj):
        if transcript_id not in self.transcripts:
            self.transcripts[transcript_id] = transcript_obj
        else:
            raise ValueError(f"Transcript {transcript_id} already exists in gene {self.gene_id}.")

class Transcript:
    def __init__(self):
        pass
    pass
    

def create_fastq_reads(genome, annotation):
    coverage = np.zeros(len(genome), dtype=int)
    for gene in annotation.genes.values():
        pass




def write_fasta(genome, output_path):
    """Write the genome sequence to a FASTA file."""
    with open(output_path, 'w') as f:
        f.write(">chr1\n")
        for i in range(0, len(genome), 60):
            f.write(genome[i:i+60] + '\n')