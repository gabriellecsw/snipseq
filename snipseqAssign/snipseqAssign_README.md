# SnipseqAssign
## Overview Flowchart
## Pipeline Description
- snipseqAssign only demultiplex reads. Please see snipseqMatch if you wish to only retrieve sequences/target/region of interest within your reads or snipseqCore if you wish to carry out both demultiplexing and sequence retrieval. 
- Users are required to prepare the following files to run snipseqAssign:
    - snipseq_assign_metadata.csv: A file containing the barcodes or barcode pairs.

## Input options for ```snipseq assign run```
Basic usage: ```snipseq assign run -i <olb.txt.gz file> -m <metadata file> -sp <spacing between front and rear barcode> -pad <spacing before front barcode>```

| Argument     | Alias  | Description                                                                        | Required |
| ------------ | ------ | ---------------------------------------------------------------------------------- | -------- |
| `--input`    | `-i`   | Path to input file.                                                                | Yes      |
| `--metadata` | `-m`   | Path to metadata CSV file with barcode pairs.                                      | Yes      |
| `--space`    | `-sp`  | Fixed number or range of bases between front and rear barcodes (e.g. 30,68 or 30). | Yes      |
| `--padding`  | `-pad` | Fixed number or range of bases before the front barcode (e.g. 30,68 or 30).        | Yes      |

## Outputs from ```snipseq assign run```
Users will have a 'main' and 'counts' csv files saved as ```snipseq_assign.csv``` and ```snipseq_assign_counts.csv```. Below are the description of each column for each csv file:
### Main output file
| Column Name         | Description                                      |
| ------------------- | ------------------------------------------------ |
| reads               | The read id.                                     |
| seq                 | The sequences of the read.                       |
| space               | The sequence between the front and rear barcode. |
| padding             | The sequence before the front barcode.           |
| front_barcode       | The sequence of the front barcode.               |
| rear_barcode        | The sequence of the rear barcode.                |
| pair                | The pair that each read has been assigned to.    |
| front_barcode_name  | The name of the front barcode.                   |
| rear_barcode_name   | The name of the rear barcode.                    |

### Counts output file
| Column Name         | Description                                      |
| ------------------  | ------------------------------------------------ |
| pair                | The pair that each read has been assigned to.    |
| front_barcode_name  | The name of the front barcode.                   |
| rear_barcode_name   | The name of the rear barcode.                    |
| counts              | The total number of reads for each pair.         |

## An example usage of snipseqAssign with paired barcodes
> snipseq assign run -i mut_lib_duplex_mapped.txt.gz -sp 701,711 -pad 38,65 -m snipseq_assign_metadata.csv
### In this example, snipseq will start finding the front barcode 38 to 65 bases after the start of the read and look for the rear barcode at 739 to 776 bases after the start of the read.
### An example of a ```snipseq_assign_metadata.csv``` file with paired barcodes separated by comma. The sequences should be specified in a 5'to 3' manner. 

```
pair,front_barcode_name,front_barcode,rear_barcode_name,rear_barcode
pair11,S515_F,TTCTAGCT,N710_R,CAGCCTCG
pair12,S520_F,AAGGCTAT,N715_R,CCTGAGAT
pair13,S520_F,AAGGCTAT,N710_R,CAGCCTCG
pair14,S515_F,TTCTAGCT,N715_R,CCTGAGAT
```
- In this example, snipseq assign will first get the reverse complement sequence for both front and rear barcodes. Next, reads with 'TTCTAGCT' and 'CGAGGCTG' (reverse complement of 'CAGCCTCG') will be assigned as pair11. This process will be repeated for each pair.
