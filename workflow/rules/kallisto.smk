# For now ignore stranded options

rule kallisto_index:
    input:
        fasta="raw_data/{genome}/genome.fasta"
    output:
        index="raw_data/{genome}/{genome}.idx"
    params:
        extra=""
    log:
        "logs/kallisto_index_{genome}.log"
    threads:
        8
    conda:
        "envs/kallisto.yml"
    shell:
        """
        kallisto index --threads {threads} --index {output.index} {input.fasta} > {log} 2>&1
        """

rule kallisto_se:
    input:
        fastq="raw_data/fastq/{accession}_1.fastq",
        index=rules.kallisto_index.output.index
    output:
        ouput_dir=directory("processed_data/{genome}/{accession}/")
    params:
        extra=""
    log:
        "logs/kallisto_{accession}{genome}.log"
    threads:
        8
    conda:
        "envs/kallisto.yml"
    shell:
        """
        kallisto quant --threads {threads} --index {input.index} --output-dir {output.output_dir} {input.fastq} > {log} 2>&1
        """

rule kallisto_pe:
    input:
        fastq=["raw_data/fastq/{accession}/{accession}_1.fastq", "raw_data/fastq/{accession}/{accession}_2.fastq"],
        index=rules.kallisto_index.output.index
    output:
        ouput_dir=directory("processed_data/{genome}/{accession}/")
    params:
        extra=""
    log:
        "logs/kallisto_{accession}{genome}.log"
    threads:
        8
    conda:
        "envs/kallisto.yml"
    shell:
        """
        kallisto quant --threads {threads} --index {input.index} --output-dir {output.output_dir} {input.fastq} > {log} 2>&1
        """