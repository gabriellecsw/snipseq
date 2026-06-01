"""
@author      Gabrielle C.
@create date 2026-03-12
@desc        Demultiplexing pipeline for different sequencing data or conditions based on users' specifications
"""

import typer
import pandas as pd
import regex as re
import numpy as np
from .logger import log, log_error
from .utils import build_pattern, rev_complement

def compute_space_range(df, metadata, use_rc=False, sample_size=20000):
    """
    Estimate the spacing between front and rear barcodes. 
    
    Parameters
    ----------
    df : pd.Dataframe = Dataframe with a 'seq' column.
    metadata : pd.Dataframe = Dataframe containing front/rear barcode information.
    use_rc : bool = If True, use reverse-complement column for rear barcodes (Nanopore / long reads). If False, use regular rear barcodes (short reads).
    sample_size : int = Number of sequences to sample for estimation. Adjust if needed.

    Returns
    -------
    median_front : float = The median distance from the start of the read to the front barcode.
    median_rear : float = The median distance from the start of the read to the rear barcode.
    """

    front_barcodes = metadata['front_barcode'].unique()
    rear_col = 'rear_barcode_rc' if use_rc else 'rear_barcode'
    rear_barcodes  = metadata[rear_col].unique()

    # Compile regex patterns
    front_patterns = [re.compile(bc) for bc in front_barcodes]
    rear_patterns  = [re.compile(bc) for bc in rear_barcodes]

    log('Sampling input to estimate barcode positions')
    sampled_seqs = df['seq'].dropna().sample(n=min(sample_size, len(df)), random_state=42)

    log('Estimating barcode positions from sampled sequences')
    front_distances = []
    rear_distances  = []

    for seq in sampled_seqs:
        if pd.isna(seq):
            continue

        # ---- FRONT ----
        front_pos = [m.start() for pat in front_patterns if (m := pat.search(seq))]
        if front_pos:
            front_distances.append(min(front_pos))

        # ---- REAR ----
        rear_pos = [m.start() for pat in rear_patterns if (m := pat.search(seq))]
        if rear_pos:
            rear_distances.append(min(rear_pos))

    median_front = np.median(front_distances)
    median_rear  = np.median(rear_distances)

    lengths = []
    
    for i in front_barcodes:
        length = len(i)
        lengths.append(length)
        
    median_fb_length = np.median(lengths)
    
    log('Computing space range between two barcodes.')
    
    space_len = median_rear - median_front - median_fb_length

    if np.isnan(space_len):
        raise ValueError("Unable to estimate barcode spacing — no matches found.")
        
    # to prevent negative spac range eg: {-5,5}
    space_len = max(0, int(round(space_len)))
        
    lower = max(0, space_len - 5)
    upper = space_len + 5
        
    if lower == upper:
        space_range = f"{lower}"
    else:
        space_range = f"{lower},{upper}"
        
    log(f"Using space range: {space_range}")

    return lower, upper, space_range

def extract_match_compiled(seq, compiled_regex):
    if pd.isna(seq):
        return None
    m = compiled_regex.search(seq)
    return m.groupdict() if m else None

def assign_short_reads(df, metadata_path):
    extracted_dfs = []
    count_records = []

    metadata = pd.read_csv(metadata_path, sep=',', header=0)
    metadata = rev_complement(metadata)

    lower, upper, space_range = compute_space_range(df, metadata, use_rc=False)

    log('Demultiplexing in progress...')
    for _, row in metadata.iterrows():
        pair_name = row['pair']
        front_barcode = row['front_barcode']
        rear_barcode = row['rear_barcode']
        
        # Build regex pattern for this barcode pair
        pattern = fr'''(?x)
        (?P<front>{front_barcode}){{e<=1}}
        (?P<space>.{{{space_range}}})
        (?P<rear>{rear_barcode}){{e<=1}}
        '''
    
        regex = re.compile(pattern)
    
        matches = df['seq'].apply(lambda x: extract_match_compiled(x, regex))
        
        # Filter and expand extracted results
        df_valid = df[matches.notna()].copy()
        extracted_cols = matches[matches.notna()].apply(pd.Series)
        
        df_combined = pd.concat([df_valid.reset_index(drop=True),
                                    extracted_cols.reset_index(drop=True)],
                                axis=1)
    
        # Add metadata info
        df_combined['pair'] = pair_name
        df_combined['front_barcode_name'] = row['front_barcode_name']
        df_combined['rear_barcode_name'] = row['rear_barcode_name']
    
        extracted_dfs.append(df_combined)
    
        # Track counts
        count_records.append({
            'pair': pair_name,
            'front_barcode_name': row['front_barcode_name'],
            'rear_barcode_name': row['rear_barcode_name'],
            'front_barcode_seq': front_barcode,
            'rear_barcode_seq': rear_barcode,
            'count': len(df_combined)
        })
    
    # Final outputs
    df_result = pd.concat(extracted_dfs).reset_index(drop=True) if extracted_dfs else pd.DataFrame()
    df_counts = pd.DataFrame(count_records)
    
    return df_result, df_counts

def assign_retrieve_short_reads(df, metadata_path, position, feature_length, space):
    extracted_dfs = []
    count_records = []

    metadata = pd.read_csv(metadata_path, sep=',', header=0)
    metadata = rev_complement(metadata)

    lower, upper, space_range = compute_space_range(df, metadata, use_rc=False)

    log('Demultiplexing & retrieving features of interest in progress...')
    for _, row in metadata.iterrows():
        pair_name = row['pair']
        front_barcode = row['front_barcode']
        rear_barcode = row['rear_barcode']

        regex = build_pattern(
                front_barcode=front_barcode,
                rear_barcode=rear_barcode,
                position=position,
                feature_length=feature_length,
                space=space,
                space_range=space_range,
                lower=lower,
                upper=upper
            )
    
        matches = df['seq'].apply(lambda x: extract_match_compiled(x, regex))
        
        # Filter and expand extracted results
        df_valid = df[matches.notna()].copy()
        extracted_cols = matches[matches.notna()].apply(pd.Series)
        
        df_combined = pd.concat([df_valid.reset_index(drop=True), extracted_cols.reset_index(drop=True)], axis=1)
        
        # Add metadata info
        df_combined['pair'] = pair_name
        df_combined['front_barcode_name'] = row['front_barcode_name']
        df_combined['rear_barcode_name'] = row['rear_barcode_name']
    
        extracted_dfs.append(df_combined)
    
        # Track counts
        count_records.append({
            'pair': pair_name,
            'front_barcode_name': row['front_barcode_name'],
            'rear_barcode_name': row['rear_barcode_name'],
            'front_barcode_seq': front_barcode,
            'rear_barcode_seq': rear_barcode,
            'count': len(df_combined)
        })
    
    # Final outputs
    df_result = pd.concat(extracted_dfs).reset_index(drop=True) if extracted_dfs else pd.DataFrame()
    df_counts = pd.DataFrame(count_records)
    
    return df_result, df_counts

def simplex_assign(df, metadata_path):
    extracted_dfs = []
    count_records = []

    metadata = pd.read_csv(metadata_path, sep=',', header=0)
    metadata = rev_complement(metadata)

    lower, upper, space_range = compute_space_range(df, metadata, use_rc=True)

    log('Demultiplexing in progress...')

    for _, row in metadata.iterrows():
        pair_name = row['pair']
        front_barcode = row['front_barcode']
        front_barcode_rc = row['front_barcode_rc']
        rear_barcode = row['rear_barcode']
        rear_barcode_rc = row['rear_barcode_rc']
            
        # Build regex pattern for this barcode pair
        pattern_fwd = fr'''(?x)
        (?P<front>{front_barcode}){{e<=1}}
        (?P<space>.{{{space_range}}})
        (?P<rear>{rear_barcode_rc}){{e<=1}}
        '''
        
        pattern_rev = fr'''(?x)
        (?P<front>{front_barcode_rc}){{e<=1}}
        (?P<space>.{{{space_range}}})
        (?P<rear>{rear_barcode}){{e<=1}}
        '''
        regex_fwd = re.compile(pattern_fwd)
        regex_rev = re.compile(pattern_rev)

        matches_fwd = df['seq'].apply(lambda x: extract_match_compiled(x, regex_fwd))
        matches_rev = df['seq'].apply(lambda x: extract_match_compiled(x, regex_rev))

        matches = matches_fwd.combine_first(matches_rev)

        # Keep only rows with matches
        mask = matches.notna()
        df_valid = df.loc[mask].copy()   # original sequences, keep 'seq' column
        df_matches = matches[mask].apply(pd.Series)  # expand front/space/rear

        # Combine original sequences with extracted info
        df_combined = pd.concat([df_valid.reset_index(drop=True), df_matches.reset_index(drop=True)], axis=1)

        # Add metadata
        df_combined['pair'] = row['pair']
        df_combined['front_barcode_name'] = row['front_barcode_name']
        df_combined['rear_barcode_name'] = row['rear_barcode_name']
    
        extracted_dfs.append(df_combined)
    
        # 7. Track counts
        count_records.append({
            'pair': pair_name,
            'front_barcode_name': row['front_barcode_name'],
            'rear_barcode_name': row['rear_barcode_name'],
            'front_barcode_seq': front_barcode,
            'rear_rc_barcode_seq': rear_barcode_rc,
            'count': len(df_combined)
            })
        # Final outputs
    df_result = pd.concat(extracted_dfs).reset_index(drop=True) if extracted_dfs else pd.DataFrame()
    df_counts = pd.DataFrame(count_records)
    
    return df_result, df_counts

def simplex_assign_retrieve(df, metadata_path, position, feature_length, space):
    extracted_dfs = []
    count_records = []

    metadata = pd.read_csv(metadata_path, sep=',', header=0)
    metadata = rev_complement(metadata)

    lower, upper, space_range = compute_space_range(df, metadata, use_rc=True)

    log('Demultiplexing & retrieving features of interest in progress...')

    for _, row in metadata.iterrows():
        pair_name = row['pair']
        front_barcode = row['front_barcode']
        front_barcode_rc = row['front_barcode_rc']
        rear_barcode = row['rear_barcode']
        rear_barcode_rc = row['rear_barcode_rc']
        
        regex_fwd = build_pattern(
                front_barcode=front_barcode,
                rear_barcode=rear_barcode_rc,
                position=position,
                feature_length=feature_length,
                space=space,
                space_range=space_range,
                lower=lower,
                upper=upper
            )

        regex_rev = build_pattern(
                front_barcode=front_barcode_rc,
                rear_barcode=rear_barcode,
                position=position,
                feature_length=feature_length,
                space=space,
                space_range=space_range,
                lower=lower,
                upper=upper
            )
        
        matches_fwd = df['seq'].apply(lambda x: extract_match_compiled(x, regex_fwd))
        matches_rev = df['seq'].apply(lambda x: extract_match_compiled(x, regex_rev))

        matches = matches_fwd.combine_first(matches_rev)
        mask = matches.notna()

        df_subset = df.loc[mask]

        df_subset = df_subset.join(matches[mask].apply(pd.Series))

        # print("Columns extracted:" , df_subset.columns.tolist())
        df_subset = df_subset.dropna(subset=['front','feature','rear']).reset_index(drop=True)
        
        # Add metadata
        df_subset['pair'] = pair_name
        df_subset['front_barcode_name'] = row['front_barcode_name']
        df_subset['rear_barcode_name'] = row['rear_barcode_name']

        extracted_dfs.append(df_subset)

        # 7. Track counts
        count_records.append({
            'pair': pair_name,
            'front_barcode_name': row['front_barcode_name'],
            'rear_barcode_name': row['rear_barcode_name'],
            'front_barcode_seq': front_barcode,
            'rear_rc_barcode_seq': rear_barcode_rc,
            'count': len(df_subset)
            })

    # Final outputs
    df_result = pd.concat(extracted_dfs).reset_index(drop=True) if extracted_dfs else pd.DataFrame()
    df_counts = pd.DataFrame(count_records)
    
    return df_result, df_counts

def duplex_assign(df, metadata_path):
    extracted_dfs = []
    count_records = []
    
    metadata = pd.read_csv(metadata_path, sep=',', header=0)
    metadata = rev_complement(metadata)

    lower, upper, space_range = compute_space_range(df, metadata, use_rc=True)

    log('Demultiplexing in progress...')
    for _, row in metadata.iterrows():
        pair_name = row['pair']
        front_barcode = row['front_barcode']
        rear_barcode_rc = row['rear_barcode_rc']
            
        # Build regex pattern for this barcode pair
        pattern = fr'''(?x)
        (?P<front>{front_barcode}){{e<=1}}
        (?P<space>.{{{space_range}}})
        (?P<rear>{rear_barcode_rc}){{e<=1}}
        '''
        
        regex = re.compile(pattern)
        
        # 4. Apply pattern to sequences
        col_name = f"extracted_{pair_name}"
        
        matches = df['seq'].apply(lambda x: extract_match_compiled(x, regex))
            
        # 5. Filter and expand extracted results
        df_valid = df[matches.notna()].copy()
        extracted_cols = matches[matches.notna()].apply(pd.Series)
            
        df_combined = pd.concat([df_valid.reset_index(drop=True),
                                    extracted_cols.reset_index(drop=True)],
                                    axis=1)
        
        # 6. Add metadata info
        df_combined['pair'] = pair_name
        df_combined['front_barcode_name'] = row['front_barcode_name']
        df_combined['rear_barcode_name'] = row['rear_barcode_name']
        
        extracted_dfs.append(df_combined)
        
        # 7. Track counts
        count_records.append({
                'pair': pair_name,
                'front_barcode_name': row['front_barcode_name'],
                'rear_barcode_name': row['rear_barcode_name'],
                'front_barcode_seq': front_barcode,
                'rear_rc_barcode_seq': rear_barcode_rc,
                'count': len(df_combined)
        })
        
    # Final outputs
    df_result = pd.concat(extracted_dfs).reset_index(drop=True) if extracted_dfs else pd.DataFrame()
    df_counts = pd.DataFrame(count_records)
        
    return df_result, df_counts

def duplex_assign_retrieve(df, metadata_path, position, feature_length, space):
    extracted_dfs = []
    count_records = []
    
    metadata = pd.read_csv(metadata_path, sep=',', header=0)
    metadata = rev_complement(metadata)

    lower, upper, space_range = compute_space_range(df, metadata, use_rc=True)

    log('Demultiplexing & retrieving features of interest in progress...')

    for _, row in metadata.iterrows():
        pair_name = row['pair']
        front_barcode = row['front_barcode']
        rear_barcode_rc = row['rear_barcode_rc']
        
        regex_fwd = build_pattern(
                front_barcode=front_barcode,
                rear_barcode=rear_barcode_rc,
                position=position,
                feature_length=feature_length,
                space=space,
                space_range=space_range,
                lower=lower,
                upper=upper
            )

        matches = df['seq'].apply(lambda x: extract_match_compiled(x, regex_fwd))
                
        mask = matches.notna()

        df_subset = df.loc[mask]

        df_subset = df_subset.join(matches[mask].apply(pd.Series))

        # print("Columns extracted:" , df_subset.columns.tolist())
        df_subset = df_subset.dropna(subset=['front','feature','rear']).reset_index(drop=True)
        
        # Add metadata
        df_subset['pair'] = pair_name
        df_subset['front_barcode_name'] = row['front_barcode_name']
        df_subset['rear_barcode_name'] = row['rear_barcode_name']

        extracted_dfs.append(df_subset)

        # Track counts
        count_records.append({
            'pair': pair_name,
            'front_barcode_name': row['front_barcode_name'],
            'rear_barcode_name': row['rear_barcode_name'],
            'front_barcode_seq': front_barcode,
            'rear_rc_barcode_seq': rear_barcode_rc,
            'count': len(df_subset)
            })

    # Final outputs
    df_result = pd.concat(extracted_dfs).reset_index(drop=True) if extracted_dfs else pd.DataFrame()
    df_counts = pd.DataFrame(count_records)
    
    return df_result, df_counts