rule make_adapter:
    output: 'results/hela_tm_tg_ont/pychopper/adapter.fa'
    params:
        VNP = 'TTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGC',
        SSP = 'GCAATACGTAACTGAACGAAGTACATT',
    shell: 'echo ">SNP\n{params.VNP}\n>SSP\n{params.SSP}" > {output}'


rule pychopper:
    input:
        adapter = rules.make_adapter.output,
    threads: 4
    params:

    shell: 'pychopper -t {threads} -b {input.adapter} -m edlib '