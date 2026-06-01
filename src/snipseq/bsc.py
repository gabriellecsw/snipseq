#!/usr/bin/env python3 #only activate this if you want to run using python on command line.

import sys
import time
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
from .utils import validate_position, validate_metadata, rev_complement, est_chunk_size, build_pattern
from . import demux 

def bsc(
    input_path: Path = typer.Option(..., "-i", "--input", help="Input file: 'fastq' or 'pod5'"),
    metadata_path: Path = typer.Option(..., "-m", "--metadata", help="Metadata CSV file"),
    sequencing_type: str = typer.Option(..., "-s", "--sequencing_type", help="Sequencing type: 'illumina' or 'nanopore'"),
    read_type: str = typer.Option(None, "-r", "--read_type", help="Read type: 'simplex' or 'duplex' [Required for Nanopore]"),
    retrieve: bool = typer.Option(False, "-ret", "--retrieve", help="Retrieve sequence of interest [Required if wanting to retrieve sequence of interest]"),
    position: str = typer.Option(None, "-p", "--position", help="Position of sequence of interest [Required if wanting to retrieve sequence of interest] Please choose from: {'P1', 'P2', 'P3', 'P4', 'P5', 'P6'}"),
    feature_length: int = typer.Option(None, "-fl", "--feature_length", help="Length of sequence of interest [Required if wanting to retrieve sequence of interest]"),
    space: int = typer.Option(None, "-sp", "--space", help="Spacing between the barcodes and the sequence of interest [Required if wanting to retrieve sequence of interest and if position is P1, P4, P6]"),
    output_path: Path = typer.Option(..., "-o","--output", help="Output file (e.g., 'SRR30861017_fwd_ss_assign.csv')")
):
    start_time = time.time()
    
    log_file = output_path.with_suffix(".log")
    setup_logger(log_file)
    
    log("Command: " + " ".join(sys.argv))
    log(f"snipseq version {snipseq.__version__}")
    
    if retrieve:
        log('Validating features positions...')
        validate_position(position, space)

    log('Validating sequencing type...')
    sequencing_type = sequencing_type.lower() #this makes sure it is lowercase

    log("Validating metadata...")
    metadata = pd.read_csv(metadata_path, sep=',', header=0)
    validate_metadata(metadata)

    if sequencing_type not in ["illumina", "nanopore"]:
        log_error("Error: --sequencing_type must be 'illumina' or 'nanopore'.")
        raise typer.Exit(code=1)

    if sequencing_type == "nanopore":
        if not read_type:
            log_error("Error: --read_type is required for nanopore.")
            raise typer.Exit(code=1)
            
        read_type = read_type.lower() #this makes sure it is lowercase
        
        df_results = []
        df_counts_list = []

        if read_type == "simplex":

            log('Running Dorado basecalling in simplex mode...')
            subprocess.run(f"dorado basecaller hac {input_path} > simplex.bam", shell=True, check=True)
            subprocess.run("samtools view -F2304 simplex.bam | "
                           "awk '{print $1 \"\\t\" $10}' | "
                           "gzip > simplex.txt.gz", shell=True, check=True)

            df_sample = pd.read_csv("simplex.txt.gz", sep="\t", header=None)
            # df_sample = pd.read_csv(input_path, sep=",", header=None) # for developer: unhash if testing

            ideal_chunk_size = est_chunk_size(df_sample)
            chunk_id = 1
            for chunk in pd.read_csv('simplex.txt.gz', chunksize=ideal_chunk_size, sep='\t', header=None, names=['read_id','seq']):
            # for chunk in pd.read_csv(input_path, chunksize=ideal_chunk_size, sep=',', header=None, names=['read_id','seq']): # for developer: unhash if testing
                log(f"Chunk {chunk_id} processed")
                chunk_id += 1
            
                if retrieve:
                    df_result, df_counts = demux.simplex_assign_retrieve(chunk, metadata_path, position, feature_length, space)

                else:
                    df_result, df_counts = demux.simplex_assign(chunk, metadata_path)

                df_results.append(df_result)
                df_counts_list.append(df_counts)
            
        elif read_type == "duplex":
            log('Running Dorado basecalling in duplex mode...')
            subprocess.run(f"dorado duplex sup {input_path} > duplex.bam", shell=True, check=True)
            subprocess.run(
                "samtools view -F2304 duplex.bam | "
                "awk '/dx:i:1/ || /dx:i:0/ {print $1 \"\\t\" $10}' | "
                "gzip > duplex.txt.gz", shell=True, check=True)

            df_sample = pd.read_csv("duplex.txt.gz", sep="\t", header=None)
            # df_sample = pd.read_csv(input_path, sep=",", header=None) # for developer: unhash if testing
            
            ideal_chunk_size = est_chunk_size(df_sample)
            chunk_id = 1

            for chunk in pd.read_csv('duplex.txt.gz', chunksize=ideal_chunk_size, sep='\t', header=None, names=['read_id','seq']):
            # for chunk in pd.read_csv(input_path, chunksize=ideal_chunk_size, sep=',', header=None, names=['read_id','seq']): # for developer: unhash if testing

                log(f"Chunk {chunk_id} processed")
                chunk_id += 1

                if retrieve:
                    df_result, df_counts = demux.duplex_assign_retrieve(chunk, metadata_path, position, feature_length, space)

                else:
                    df_result, df_counts = demux.duplex_assign(chunk, metadata_path)

                df_results.append(df_result)
                df_counts_list.append(df_counts)

        else:
            log_error("Error: read_type must be 'simplex' or 'duplex'.")
            raise typer.Exit(code=1)

        # Concatenate all results and counts after processing all chunks
        final_df_result = pd.concat(df_results, ignore_index=True)
        final_df_counts = pd.concat(df_counts_list, ignore_index=True)
        
        #as_index=False makes sure that the pairs column stays as a column instead of an index
        #this step is to group the pairs together because they were separated when reading in by chunks
        final_df_counts = final_df_counts.groupby(['pair','front_barcode_name', 'rear_barcode_name',
                                                   'front_barcode_seq','rear_rc_barcode_seq'], as_index=False)['count'].sum()
        
    # If input is short read, this command will be executed.
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

        ideal_chunk_size = est_chunk_size(df_sample)
        
        chunk_id = 1
        for chunk in pd.read_csv("reads.txt.gz", chunksize=ideal_chunk_size, sep='\t', header=None, names=['read_id','seq']):
        # for chunk in pd.read_csv(input_path, chunksize=ideal_chunk_size, sep=',', header=None, names=['read_id','seq']): # for developer: unhash if testing
            log(f"Chunk {chunk_id} processed")
            chunk_id += 1
            
            if retrieve:
                df_result, df_counts = demux.assign_retrieve_short_reads(chunk, metadata_path, position, feature_length, space)

            else:
                df_result, df_counts = demux.assign_short_reads(chunk, metadata_path)
            
            log("Next chunk...be patient my friend")

            df_results.append(df_result)
            df_counts_list.append(df_counts)
            
        final_df_result = pd.concat(df_results, ignore_index=True)
        final_df_counts = pd.concat(df_counts_list, ignore_index=True).groupby(['pair', 'front_barcode_name', 'rear_barcode_name',
                                                                                'front_barcode_seq', 'rear_barcode_seq'], as_index=False)['count'].sum()
        
    final_df_result.to_csv(output_path, index=False)
    log(f"Demultiplexing complete. Main output saved to {output_path}")

    suffix = output_path.suffix or ".csv"  # default to .csv if missing
    counts_name = output_path.with_name(output_path.stem + "_counts" + suffix)

    final_df_counts.to_csv(counts_name, index=False)
    log(f"Demultiplexing complete. Counts output saved to {counts_name}")
        
    end_time = time.time()

    total_reads = len(final_df_result)
    runtime = time.time() - start_time
    
    log(f"Processed {total_reads:,} reads in {runtime:.2f}s")
    log(f"Speed: {total_reads/runtime:,.0f} reads/sec")