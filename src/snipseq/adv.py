#!/usr/bin/env python3 #only activate this if you want to run using python on command line.

import sys
import time
import tomli
import typer
import snipseq
import logging
import subprocess
import numpy as np 
import regex as re 
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq

from .logger import setup_logger, log, log_error
from .utils import est_chunk_size

def demux(pair, seq):
    m = re.match(pair, seq)
    return m.groupdict() if m else None
    
def adv(
    input_path: Path = typer.Option(..., "-i", "--input", help="Input file: 'fastq' or 'pod5'"),
    arrangement_path: Path = typer.Option(..., "-a", "--arrangement", help="Arrangement TOML file"),
    sequencing_type: str = typer.Option(..., "-s", "--sequencing_type", help="Sequencing type: 'illumina' or 'nanopore'"),
    read_type: str = typer.Option(None, "-r", "--read_type", help="Read type: 'simplex' or 'duplex' [Required for Nanopore]"),
    output_path: Path = typer.Option(..., "-o","--output", help="Output file (e.g., 'SRR30861017_fwd_ss_assign.csv')")
):

    start_time = time.time()
    
    log_file = output_path.with_suffix(".log")
    setup_logger(log_file)
    
    log("Command: " + " ".join(sys.argv))
    log(f"snipseq version {snipseq.__version__}")

    log('Validating sequencing type...')
    sequencing_type = sequencing_type.lower() #this makes sure it is lowercase

    if sequencing_type not in ["illumina", "nanopore"]:
        log_error("Error: --sequencing_type must be 'illumina' or 'nanopore'.")
        raise typer.Exit(code=1)

    with open(arrangement_path, mode="rb") as fp:
        config_arrangement = tomli.load(fp)

    config_arrangement = pd.DataFrame.from_dict(config_arrangement['arrangements'], orient='index', columns=['value']).reset_index()
    config_arrangement.columns = ['features', 'value']
    print(config_arrangement)

    features = ['(?x)']  
    
    for _, row in config_arrangement.iterrows():
        feat = row['features']
        val = row['value']
    
        if isinstance(val, str) and val.isdigit(): # if input is fixed number
            feature = fr"(?P<{feat}>.{{{val}}})"
    
        elif isinstance(val, str) and ',' in val and all(part.strip().isdigit() for part in val.split(',')): # if input is range
            feature = fr"(?P<{feat}>.{{{val}}})"
    
        elif isinstance(val, str): # if input is sequence, fuzzy match with {e<=1}
            feature = fr"(?P<{feat}>{val}){{e<=1}}"
    
        else:
            raise ValueError(f"Unsupported value format for field: {feat} -> {val}")
    
        features.append(feature)
    
    full_features = '\n'.join(features)
    #print(full_features)

    if sequencing_type == "nanopore":
        if not read_type:
            log_error("Error: --read_type is required for nanopore.")
            raise typer.Exit(code=1)
            
        read_type = read_type.lower()

        if read_type == "simplex":

            log('Running Dorado basecalling in simplex mode...')
            subprocess.run(f"dorado basecaller hac {input_path} > simplex.bam", shell=True, check=True)
            subprocess.run("samtools view -F2304 simplex.bam | "
                           "awk '{print $1 \"\\t\" $10}' | "
                           "gzip > simplex.txt.gz", shell=True, check=True)

            df_sample = pd.read_csv("simplex.txt.gz", sep="\t", header=None)
            # df_sample = pd.read_csv(input_path, sep=",", header=None) # for developer: unhash if testing
        
        elif read_type == "duplex":
            log('Running Dorado basecalling in duplex mode...')
            subprocess.run(f"dorado duplex sup {input_path} > duplex.bam", shell=True, check=True)
            subprocess.run(
                "samtools view -F2304 duplex.bam | "
                "awk '/dx:i:1/ || /dx:i:0/ {print $1 \"\\t\" $10}' | "
                "gzip > duplex.txt.gz", shell=True, check=True)

            df_sample = pd.read_csv("duplex.txt.gz", sep="\t", header=None)
            # df_sample = pd.read_csv(input_path, sep=",", header=None) # for developer: unhash if testing

    elif sequencing_type == "illumina":

        subprocess.run(
            f"""zcat {input_path} |
            awk 'NR%4==1 {{id=$1; sub(/^@/, "", id)}} 
                 NR%4==2 {{print id "\\t" $0}}' | 
            gzip > reads.txt.gz""", shell=True, check=True)
        
        df_results = []
        df_counts_list = []

        log("Loading input file by chunk")
        
        df_sample = pd.read_csv("reads.txt.gz", sep="\t", header=None)
        # df_sample = pd.read_csv(input_path, sep=",", header=None) # for developer: unhash if testing


    df_results = []

    log("Loading input file by chunk")
    
    ideal_chunk_size = est_chunk_size(df_sample)
    chunk_id = 1
    for chunk in pd.read_csv('simplex.txt.gz', chunksize=ideal_chunk_size, sep='\t', header=None, names=['read_id','seq']):
    # for chunk in pd.read_csv(input_path, chunksize=ideal_chunk_size, sep='\t', header=None, names=['read_id','seq']): # for developer: unhash if testing
        log(f"Chunk {chunk_id} processed")
        chunk_id += 1

        chunk = chunk.rename(columns={0:'read_id',1:'seq'})
        chunk['retseq'] = chunk['seq'].apply(lambda x: demux(full_features, x)) #retrieve the sequences of each feature
        df_retseq = chunk[chunk['retseq'].notna()].reset_index(drop=True)
        df_retseq_ext = pd.DataFrame(df_retseq['retseq'].values.tolist()) # expand the column that contains the extracted sequences.
        df_retseq_com = pd.concat([df_retseq, df_retseq_ext], axis=1) # paste the read_id info back into the dataframe

        log("Next chunk...be patient my friend")

        df_results.append(df_retseq_com)

    final_df_result = pd.concat(df_results, ignore_index=True)
    final_df_result.to_csv(output_path, index=False)
    log(f"Demultiplexing complete. Main output saved to {output_path}")
        
    end_time = time.time()

    total_reads = len(final_df_result)
    runtime = time.time() - start_time
    
    log(f"Processed {total_reads:,} reads in {runtime:.2f}s")
    log(f"Speed: {total_reads/runtime:,.0f} reads/sec")