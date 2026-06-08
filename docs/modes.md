# Modes

PERREO includes three execution modes selected according to sequencing technology.

<img width="2437" height="875" alt="PERREO mode logos" src="https://github.com/user-attachments/assets/5240009a-2645-4fd1-9fff-7e5ced4cf3b5" />

| Requirement | SR-SE | SR-PE | LR |
| --- | --- | --- | --- |
| Single-end FASTQ files | Yes | No | Yes |
| Paired-end FASTQ files | No | Yes | No |
| Reference genome FASTA | Yes | Yes | Yes |
| Genome GTF | Yes | Yes | No |
| Repeats GTF | Yes | Yes | Yes |
| Batch correction setting | Yes | Yes | Yes |
| Differential expression method | Yes | Yes | Yes |
| Duplicate removal option | Yes | Yes | No |

## Default Parameters

| Argument | Default |
| --- | --- |
| `-threads` | `8` |
| `-ram` | `32` |
| `-trimming` | `simple` |
| `-mismatch_align` | `0.05` |
| `-trimming_quality_threshold` | `30` |
| `-min_length_trim` | `16` |
| `-initial_trim_read1` | `0` |
| `-initial_trim_read2` | `0` |
| `-initial_trim_read` | `0` |
| `-k_num` | `2` |
| `-log2FC` | `1` |
| `-FDR` | `0.05` |
| `-prediction_model` | `no` |
| `-batch` / `-batch_effect` | `no` |

Long-read mode does not use short-read trimming parameters because adapter and barcode processing should be handled during basecalling.
