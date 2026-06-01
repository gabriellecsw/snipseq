"""
@author      Gabrielle C.
@create date 2026-03-12
@desc        Helper functions for Snipseq
"""

import typer
import numpy as np
import regex as re
import pandas as pd
from Bio.Seq import Seq
from .logger import log, log_error

def validate_metadata(metadata, cli_mode=True):
    REQUIRED_COLUMNS = {
        "pair",
        "front_barcode_name",
        "front_barcode",
        "rear_barcode_name",
        "rear_barcode"
    }

    missing = REQUIRED_COLUMNS - set(metadata.columns)
    
    # check for missing columns
    if missing:
        msg = f"Metadata file is missing required columns: {', '.join(sorted(missing))}"
        if cli_mode:
            log_error(msg)
            raise typer.Exit(code=1)
        else:
            raise ValueError(msg)

    # check for empty values
    for col in REQUIRED_COLUMNS:
        if metadata[col].isna().any():
            msg = f"Metadata column '{col}' contains missing values"
            if cli_mode:
                log_error(msg)
                raise typer.Exit(code=1)
            else:
                raise ValueError(msg)
            
def validate_position(position, space, cli_mode=True):
    if position is None:
        msg = "position cannot be None"
        if cli_mode:
            log_error(msg)
            raise typer.Exit(code=1)
        else:
            raise ValueError(msg)
        
    position = position.upper()
    VALID_POSITIONS = {"P1", "P2", "P3", "P4", "P5", "P6"}

    if position not in VALID_POSITIONS:
        msg = f"position must be one of {sorted(VALID_POSITIONS)}, got: {position}"
        if cli_mode:
            log_error(msg)
            raise typer.Exit(code=1)
        else:
            raise ValueError(msg)

    SPACED_POSITIONS = {"P1", "P4", "P6"}
    if position in SPACED_POSITIONS and space is None:
        msg = "Error: --space is required for spaced patterns."
        if cli_mode:
            log_error(msg)
            raise typer.Exit(code=1)
        else:
            raise ValueError(msg)

def rev_complement(metadata):

    log('Generating reverse complement barcodes.')

    metadata['front_barcode_rc'] = metadata['front_barcode'].apply(lambda x: str(Seq(x).reverse_complement()))
    metadata['rear_barcode_rc'] = metadata['rear_barcode'].apply(lambda x: str(Seq(x).reverse_complement()))

    return metadata

def est_chunk_size(df_sample, sample_n = 100000):
    """
    df_sample: input file after preprocessing.
    sample_n: sample size required to estimate chunk size. Larger sample gives better estimate, adjust if needed.
    """

    log('Estimating chunk size.')

    # df_sample = pd.read_csv(input_path, nrows=sample_n, sep="\t", header=None)
    
    if len(df_sample) == 0:
        raise ValueError("Input file appears empty.")
    
    avg_per_row_bytes = df_sample.memory_usage(deep=True).sum() / len(df_sample) # Estimate memory per row
        
    mem_limit_bytes = 500 * 1024 * 1024 # Set memory budget for each chunk (e.g., 500 MB)
    ideal_chunk_size = max(1, int(mem_limit_bytes / avg_per_row_bytes))

    log(f"Ideal chunk size estimated: {ideal_chunk_size:,}")

    return ideal_chunk_size

def build_pattern(front_barcode,
                  rear_barcode,
                  position,
                  feature_length,
                  space=None,
                  space_range=None,
                  lower=None,
                  upper=None):
    """
    Regex pattern builders for Snipseq barcode extraction
    """

    if feature_length is None:
        raise ValueError("feature_length is required for all positions")

    if position in {"P1", "P6"} and space is None:
        raise ValueError(f"{position} requires space")

    if position == "P4":
        if None in (space, lower, upper):
            raise ValueError("P4 requires space, lower, and upper")

    else:
        if space_range is None:
            raise ValueError(f"{position} requires space_range")

    log(f"Building regex pattern for position {position}")

    if position == "P2":
        pattern = fr'''(?x)
        (?P<feature>.{{{feature_length}}})
        (?P<front>{front_barcode}){{e<=1}}
        (?P<inter>.{{{space_range}}})
        (?P<rear>{rear_barcode}){{e<=1}}
        '''

    elif position == "P3":
        pattern = fr'''(?x)
        (?P<front>{front_barcode}){{e<=1}}
        (?P<feature>.{{{feature_length}}})
        (?P<rear>{rear_barcode}){{e<=1}}
        '''

    elif position == "P5":
        pattern = fr'''(?x)
        (?P<front>{front_barcode}){{e<=1}}
        (?P<inter>.{{{space_range}}})
        (?P<rear>{rear_barcode}){{e<=1}}
        (?P<feature>.{{{feature_length}}})
        '''

    elif position == "P1":
        pattern = fr'''(?x)
        (?P<feature>.{{{feature_length}}})
        (?P<up_space>.{{{space}}})
        (?P<front>{front_barcode}){{e<=1}}
        (?P<inter>.{{{space_range}}})
        (?P<rear>{rear_barcode}){{e<=1}}
        '''

    elif position == "P4":
        user_space = space
        lower_inter = max(0, lower - (user_space + feature_length))
        upper_inter = max(lower_inter, upper - (user_space + feature_length)) # ensure upper >= lower

        pattern = fr'''(?x)
        (?P<front>{front_barcode}){{e<=1}}
        (?P<inter_user>.{{{user_space}}})
        (?P<feature>.{{{feature_length}}})
        (?P<inter>.{{{lower_inter},{upper_inter}}})
        (?P<rear>{rear_barcode}){{e<=1}}
        '''

    elif position == "P6":
        pattern = fr'''(?x)
        (?P<front>{front_barcode}){{e<=1}}
        (?P<inter>.{{{space_range}}})
        (?P<rear>{rear_barcode}){{e<=1}}
        (?P<down_space>.{{{space}}})
        (?P<feature>.{{{feature_length}}})
        '''

    else:
        raise ValueError(f"Unsupported position: {position}")

    return re.compile(pattern)