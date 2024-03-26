# Look for SRR .txt in directory
import glob
import pandas as pd
import numpy as np
from pathlib import Path

# FUNCTIONS
def get_table(SRA_run_info: Path) -> pd.DataFrame:
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
SRA_info=get_table(config["SraRunInfo"])
SRA_SE=get_SE(SRA_info)
SRA_SE_acc=SRA_SE["Run"].tolist()
SRA_PE=get_PE(SRA_info)
SRA_PE_acc=SRA_PE["Run"].tolist()


rule download_SE:
    output:
        fastq=temp("raw_data/SE/{accession}_1.fastq")
    log:
        "logs/fasterq_dump/{accession}_fasterq.log"
    threads:
        8
    params:
        download_folder=config["SRA_download_SE"],
        accession=lambda wildcards: SRA_SE_acc
    conda:
        "../envs/SRA.yml"
    shell:
        """
        fasterq-dump --threads {threads} --split-files {params.accession} -O {params.download_folder} > {log} 2&>1
        """

rule gzip_SE:
    input:
        SE_fastq=rules.download_SE.output.fastq
    output:
        fastq_gz="raw_data/SE/{accession}_1.fastq.gz"
    threads:
        8
    params:
        extra=""
    shell:
        """
        gzip {input.SE_fastq}
        """

rule download_PE:
    output:
        fastqs=temp(expand("raw_data/PE/{accession}_{PE}.fastq", PE = [1, 2], accession = SRA_PE_acc)),
    threads:
        8
    params:
        download_folder=config["SRA_download_PE"],
        accessions=lambda wildcards: SRA_PE_acc
    conda:
        "../envs/SRA.yml"
    shell:
        """
        fasterq-dump --threads {threads} --split-files {params.accessions} -O {params.download_folder} > {log} 2&>1 
        """

rule gzip_PE:
    input:
        fastqs=rules.download_PE.output.fastqs
    output:
        fastqs_gz=expand("raw_data/PE/{accession}_{PE}.fastq.gz", PE = [1, 2], accession = SRA_PE_acc)
    threads:
        8
    params:
        extra=""
    shell:
        """
        gzip {input.fastqs}
        """