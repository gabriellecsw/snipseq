# Snipseq: A Demultiplexing and Sequence Retrieval Tool for User-Defined Barcode

## Key Features
- Flexible Barcode Usage: Accepts user-defined barcode pairs, enabling researchers to repurpose barcodes already available in the lab.
- Compatible with Various Sequencing Services (e.g., Plasmidsaurus) and Platforms: Eliminating the need for custom library prep kits. 
- Mapping and Target Sequence Extraction: Reduce the need for separate preprocessing steps. 
- Cost-Efficiency: Enables multiplexing of replicates and conditions, allowing pooled sequencing of multiple samples in one run. Users are able to then link each read to its experimental condition via barcode mapping.
- Flexible Input Files: Accepts `fastq` and `pod5` files as input.

Snipseq is designed to remain accessible even to users with limited bioinformatics experience. At its simplest, only the sequencing reads and a metadata file describing the barcode pairs are required.

![overview](figures/snipseq_pipeline.png)

## INSTALLATION
```r
git clone https://github.com/gabriellecsw/snipseq.git
cd snipseq
pip install .
```

## INPUT FILES
To support a wide range of sequencing platforms, `Snipseq` currently accepts raw sequencing reads generated from Nanopore (long-read) and Illumina (short-read) technologies.

For Illumina raw read files, `Snipseq` expects a `fastq` input file. It extracts the read identifiers and sequences and organises them into a dataframe for downstream processing.

For Nanopore raw read files, `Snipseq` expects a `pod5` input file. It will first perform basecalling using either `simplex` or `duplex` mode from the `Dorado` package, depending on the user’s specification, followed by extraction of read identifiers and sequences.


## HOW TO START
To accommodate different experimental designs and analysis goals, Snipseq provides two operating modes: `basic` mode and `advanced` mode.

Both modes support demultiplexing and sequence retrieval, allowing users to choose the level of complexity appropriate for their experiment.

![modes](figures/snipseq_example_usages.png)

### Basic mode
The `basic` mode is designed for simple experimental setups involving a single feature of interest.

In this mode, users may perform demultiplexing alone or demultiplexing with feature sequence retrieval using the `-ret` argument (`False` for demultiplexing only and `True` to demultiplex and retrieve the feature sequence). 

#### Demultiplex only
Below is the example to run snipseq on the basic mode with `-ret False`:

```r
snipseq bsc
-i SRR30861029_2.fastq.gz # Your sample file that requires demultiplexing.
-m metadata_ime.csv # A metadata file containing your paired barcode sequences.
-s illumina # The sequencing platform used.
-ret False # Whether to retrieve sequence of interest. The default is False, if True, please refer to next subsection of the basic mode for more information.
-o SRR30861029_2_ss_assign_ret.csv # Your output directory.
```

**Note**: If the sequencing type is `nanopore`, the `read_type` parameter must be specified. Users can choose either `simplex` or `duplex`.

#### Demultiplex and retrieve sequence of interest
When running the basic mode with `-ret True`, there are several options to choose from depending on the postion of your sequence of interest:
![overview](figures/snipseq_bsc_position_options.png)

Below is an example to run snipseq on the basic mode with `-ret True`:

```r
snipseq bsc
-i SRR30861029_2.fastq.gz # Your sample file that requires demultiplexing.
-m metadata_ime.csv # A metadata file containing your paired barcode sequences.
-s illumina # The sequencing platform used.
-ret True # Whether to retrieve sequence of interest.
-p P4 # The position of your sequence of interest on your sequenced fragement.
-fl 18 # The length of your sequence of interest.
-sp 0 # See notes above for sp.
-o SRR30861029_2_ss_assign_ret.csv # Your output directory
```

**Note**: If you are using P1, P4 and P6, you must specify the spacing accordingly using the `sp` parameter. 

### Advance mode
The `advanced` mode is intended for more complex experiments involving multiple features or more elaborate read architectures.

In this mode, `Snipseq` automatically performs both demultiplexing and sequence retrieval.

#### Preparing the `.toml` arrangement file
This `.toml` arrangement file contains the arrangement of specific features within the read structure. Refer to [Figure(B)](#HOW-TO-START) on how to prepare a `.toml` arrangement file.

```r
snipseq adv
-i SRR30861029_2.fastq.gz # Your sample file that requires demultiplexing
-a arrangement.toml # A .toml file containing the arrangement of features you would like to extract
-s illumina # The sequencing platform used.
-o SRR30861029_2_ss_assign_ret.csv # Output directory
```

## OUTPUT FILES
Snipseq generates three primary output files:
- The `main_output.csv` file contains the read identifiers and their corresponding barcode assignments.
- The `counts_summary.csv` file summarises the number of reads assigned to each barcode pair, allowing users to quickly evaluate barcode distribution.
- The `log_summary.csv` file records run statistics, including runtime information and the number of assigned and unassigned reads. # will be available soon.