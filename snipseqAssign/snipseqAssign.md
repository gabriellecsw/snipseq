# SnipseqAssign

## Reference Diagram 

A double-ended barcode with different flanks and barcode sequences for front and rear barcodes is described here.

> 5′--- padding --- front_barcode--- trail_frt --- UMI_fwd --- fwd_primer --- target_flank1 --- target --- target_flank2 --- rev_primer --- UMI_rev --- trail_rev --- rear_barcode ---3′ 

## Reads Arrangement Options

The table below describes the arrangement options in more detail. Users are free to add or drop the options according to their experimental design and what they want to extract from the sequencing data. For each option, users can specify either the exact sequence, single postition, or a range of position to match.

| Option        | Description                                                                                                                       |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| padding       | The region upstream of the front custom barcode.                                                                                  |
| front_barcode | The front custom barcode.                                                                                                         |
| trail_frt     | The trailing flank of the front custom barcode.                                                                                   |
| UMI_fwd       | The length of the UMI used in the forward primer.                                                                                 |
| fwd_primer    | The sequence or length of the forward primer.                                                                                     |
| target_flank1 | The region between the forward primer and the target.                                                                             |
| target        | The region of interest/ variation. To extract the sequences of the variable region, use the length of the sequence as input.      |
| target_flank2 | The region between the target and the reverse primer.                                                                             |
| rev_primer    | The sequence or length of the forward primer. If using sequences, please use the reverse complement of the reverse primer.        |
| UMI_rev       | The length of the UMI used in the reverse primer.                                                                                 |
| trail_rev     | The trailing flank of the rear custom barcode. If using sequences, please use the reverse complement of the rear custom barcode.  |
| rear_barcode  | The rear custom barcode. If using sequences, please use the reverse complement of the rear custom barcode.                        |

## Example Usage
In this example, the target is 'NTCATGNN'. By leveraging the flexibility of snipseq, we split the ```target``` into 3 parts, ``` up_v ```, ```cons_v``` and ``` down_v ```.

```
pair11_e1 = r"""(?x)
(?P<padding>.{30,68})
(?P<front_barcode>TTCTAGCT){e<=1}
(?P<trail_frt>TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG){e<=1}
(?P<UMI_fwd>.{2})
(?P<fwd_primer>GTGTTCACTAGCAACCTCAAAC){e<=1}
(?P<target_flank1>ATTC){e<=1}
(?P<up_v>.{1})
(?P<cons_v>TCATG){e<=1}
(?P<down_v>.{2})
(?P<target_flank2>.{575,585})
(?P<rev_primer>GAGAGAGTCACCACATACGAAG){e<=1}
(?P<UMI_rev>.{2})
(?P<trail_rev>CTGTCTCTTATACACATCTCCGAGCCCACGAGAC){e<=1}
(?P<rear_barcode>ATCTCAGG){e<=1}
"""

```
## Description for the options in the example usage:
- Assigning an exact sequence           ```(?P<front_barcode>TTCTAGCT){e<=1}``` : Find these sequences in the read after ```padding``` by allowing up to one mismatch and assign them as the ```front_barcode```.
- Assigning a range of positions        ```(?P<padding>.{30,68})```             : Assign sequences from the start of the read up to 30bp or 68bp as ```padding```.
- Assigning a fixed number of positions ```(?P<UMI_fwd>.{2})```                 : Assign 2 bases after ```trail_frt``` as ```UMI_fwd```.
