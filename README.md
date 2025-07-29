# Snipseq: A Demultiplexing and Sequence Retrieval Tool for User-Defined Barcode in Reporter Assays
# Key Features
- Flexible Barcode Usage: Accepts user-defined barcode pairs, enabling researchers to repurpose barcodes already available in the lab.
- Compatible with Various Sequencing Services (e.g., Plasmidsaurus) and Platforms: Eliminating the need for custom library prep kits. 
- Mapping, Target Sequence Extraction and Adaptor Trimming: Reduce the need for separate preprocessing steps. 
- Cost-Efficiency: Enables multiplexing of replicates and conditions, allowing pooled sequencing of multiple samples in one run. Users are able to then link each read to its experimental condition via barcode mapping.
- Standard input files: Accepts fastq files or raw sequencing files as input.

# Snipseq is broken down into two pipelines, snipseqMatch and snipseqAssign. Please see the links below for detailed documentation:
- [SnipseqMatch](https://github.com/gabriellecsw/snipseq/tree/main/snipseqMatch): Generation of Oligo-Barcode Library, the input file for snipseqAssign.
- [SnipseqAssign](https://github.com/gabriellecsw/snipseq/tree/main/snipseqAssign): Demultiplexing and target sequence retrieval.
