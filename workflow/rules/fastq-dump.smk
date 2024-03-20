# Look for SRR .txt in directory
import glob
import pandas as pd
import numpy as np
from pathlib import Path

accessions_list=[]
for filename in (glob.glob("../SRR/SRR_*.txt")):
    with open(filename, 'r') as f:
        for line in f:
            accessions_list.append(line.strip())

def get_table(SRA_run_info: Path) -> pd.Data:
    """
    read SraRunInfo.csv  
    -------------------
    params:
        SRA_run_info:
            SraRunInfo.csv from NCBI SRA run selector
            with information on LibraryLayout and Run
    -------------------
    returns:
        data
            pd.DataFrame object with all columns from
            params.SRA_run_info
            NOT VERIFIED!
        
    """
    data = pd.read_csv(SRA_run_info)

    return data

def get_SE(SRA_data: pd.DataFrame) -> pd.DataFrame:
    """
    get single end libaries from SRA_run_info
    -------------------
    params:
        SRA_data
            pd.DataFrame object from previous function
    -------------------
    returns:
        SE_files
            pd.DataFrme object containg only single-end
            library
    """
    return SRA_data.loc[SRA_data["LibraryLayout"] == "SINGLE"]

def get_PE(SRA_data: pd.DataFrame) -> pd.DataFrame:
    """
    get single end libaries from SRA_run_info
    -------------------
    params:
        SRA_data
            pd.DataFrame object from previous function
    -------------------
    returns:
        PE_files
            pd.DataFrme object containg only single-end
            library
    """
    return SRA_data.loc[SRA_data["LibraryLayout"] == "PAIRED"]

# VARIABLES:
DOWNLOAD_PATH=config["SRA_download"]
SRA_data=get_table(DOWNLOAD_PATH)
SRA_SE=get_SE(SRA_data)
SRA_SE_acc=SRA_SE["run"].tolist()
SRA_PE=get_PE(SRA_data)
SRA_PE_acc=SRA_PE["run"].tolist()


rule get_SE_fastq:
    output:
        fastq="raw_data/SE/{accession}_1.fastq"
        fastq_gz="raw_data/SE/{accession}_1.fastq.gz"
    log:
    threads:
    params:
        download_folder=config["SRA_download"]
        accession=lamba wildcards: SRA_SE_acc
    conda:
    shell:
        """
        fasterq-dump --threads {threads} --split-files {params.accession} -O {output.output_dir}
        gzip {output.fastq} 
        """

rule get_PE_fastq:
    output:
        fastq_1="raw_data/PE/{accession}_1."
        fastq_gz=
    log:
    threads:
    conda:
    shell:




#rule get_fastq_split_files:
#    #input:
#    #    acc=expand("{accession}", accession=accessions_list)
#    output:
#        output_dir=directory("../raw_data/fastq/{accession}/"),
#        fastq="../raw_data/fastq/{accession}/{accession}.fastq.log"
#    log:
#        "logs/fasterq_dump/{accession}_fasterq.log"
#    threads:
#        8
#    conda:
#        "envs/SRA.yml"
#    shell:
#        """
#        fasterq-dump --threads {threads} --split-files {wildcards.accession} -O {output.output_dir} > {log} 2&>1
#        """

rule download_fastq_single:
    output:
        output_dir=directory("../raw_data/fastq/{accession}/"),
        singleFastq = lambda wildcards:
    log:
        "logs/fasterq_dump/{accession}_fasterq.log"
    threads:
        8
    conda:
        "envs/SRA.yml"
    shell:
        """
        fasterq-dump --threads {threads} --split-files {wildcards.accession} -O {output.output_dir}
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