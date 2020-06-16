# demux_i7

Demultiplex fasta based on the header. This version uses i7 index from paired
barcode.

You can set edit distance between observed barcode and barcodes specified in the BARCODES file.

```
  -h, --help            show this help message and exit
  -f FASTQ_FOLDER, --fastq_folder FASTQ_FOLDER
                        Folder containing fastq.gz files for demultiplex
  -o OUTPUT, --output OUTPUT
                        Output folder
  -b BARCODES, --barcodes BARCODES
                        File with barcodes that should be demultiplexed and
                        samle names. Other barcodes will be discarded
  -d DISTANCE, --distance DISTANCE
                        Distance between observed and target barcode
```