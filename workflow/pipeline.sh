#!/bin/bash

snakemake --keep-going --cores 12 --use-conda --rerun-incomplete --debug-dag
