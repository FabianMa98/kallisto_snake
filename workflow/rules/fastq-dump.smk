# Look for SRR .txt in directory
import glob

accessions_list=[]
for filename in (glob.glob("../SRR/SRR_*.txt")):
    with open(filename, 'r') as f:
        for line in f:
            accessions_list.append(line.strip())

rule get_fastq_split_files:
    #input:
    #    acc=expand("{accession}", accession=accessions_list)
    output:
        output_dir=directory("../raw_data/fastq/{accession}/"),
        fastq="../raw_data/fastq/{accession}/{accession}.fastq.log"
    log:
        "logs/fasterq_dump/{accession}_fasterq.log"
    threads:
        8
    conda:
        "envs/SRA.yml"
    shell:
        """
        fasterq-dump --threads {threads} --split-files {wildcards.accession} -O {output.output_dir} > {log} 2&>1
        """

#rule get_fastq:
#    input:
#        acc=expand("{accession}", accession=accessions_list)
#    output:
#        output_dir=directory("../raw_data/fastq/")
#    log:
#        "logs/fasterq_dump/{accession}_fasterq.log"
#    threads:
#        8
#    conda:
#        "envs/SRA.yml"
#    shell:
#        """
#        fasterq-dump --threads {threads} --split-files {input.acc} -O {output.output_dir} -t tmp/ > {log} 2&>1
#        """