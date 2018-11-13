
# Bisulfite Treatment of RNA for Determination of 5-methylcytosine

For information on the process, see [here](https://link.springer.com/protocol/10.1007%2F978-1-4939-6807-7_8). The basic idea is that bisulfite treatment of RNA modifies cytosine:

![bisulphite chemistry](/assets/20181012_bisulphite/bisulfite_chemistry.gif)

Since it converts unmethylated cytosine to uracil faster than it convert 5-methylcytosine to thymine, we can selectively convert only the unmethylated cytosine.

This means that when we sequence it, cytosines that *remain* cytosine are methylated, whereas those converted to uracil (seen as thymine) are not!

# Notes on our experiment from Nara

Seven separate experiments are being sequenced simultaneously (barcodes 1, 2, 4, 5, 6, 7, & 12 from [NEB E7335](https://www.neb.com/-/media/catalog/datacards-or-manuals/manuale7335.pdf)) through multiplexing. Fortunately, the iSeq is able to demultiplex each one on its own (we give it indices in the local run manager corresponding to the barcodes below). If it fails to demultiplex, e.g. we see variance in sample quantity when they should all be relatively equal, you need to go into the unknown reads and manually separate by barcode.

| Sample Name | Barcode 
| :--- | -------
| *in vivo* 1 | ATCACG 
| *in vivo* 2 | CGATGT 
| *in vitro* 1 | TGACCA
| *in vitro* 2 | ACAGTG 
| CTRL | GCCAAT
| NSUN2-kd | CAGATC
| DMNT2-kd | CTTGTA


<u>Descriptions:</u>

The *in vitro* replicates are **negative controls**, since they were not grown in cells with the machinery necessary for methylation.

CTRL = **transfection control**: scrambled siRNA that was transfected into the cell and shouldn't have any impact on methylation. This is to confirm that the mere process of transfection didn't alter methylation in some way.

NSUN2-kd = **NSUN2 knockdown** (a methyltransferase that catalyzes m5C formation, specifically on tRNA). We want to see if NSUN2 is the methyltransferase involved in alkylating EBER1.

DMNT2-kd = **DNMT2 knockdown** (A DNA & RNA methyltransferase that might methylate ncRNA). Again, we want to see if DNMT2 is related to EBER1's methylation.

# Setup

## User-Defined Variables

Change these before proceeding.

**Note:** The HITS-CLIP adapters, RL5 and RL3, were used to pull out the *in vivo* and *in vitro* reads only. The EBER1s from the siRNA transfection control and methyltransferase knockdowns were obtained simply through PCR, and so should not have adapters.


```python
# Note: do NOT include an end "/"
output_dir = '.'

# Define separately for each sample
paired_end = [True, True, True, True, True, True, True]
run_fastqc = [True, True, True, True, True, True, True]
has_adapters = [True, True, True, True, False, False, False]

# This will not only serve as a convenient means for indicating the input
# files, but will be used as file prefixes for output files
file_prefixes = ['Vitro-1_S2', 'Vitro-2_S4', 'Vivo-1_S1', 'Vivo-2_S3',
             'CTRL_S5', 'NSUN2_S6', 'DMNT2_S7']

# Main list from which files will be processed.
# Single-End should be in the format:
#  ['f1', 'f2', ..., 'fn']
# Paired-End should be in the format:
#  [('f1r1', 'f1r2'), ('f2r1', f2r2'), ..., ('fnr1', 'fnr2')]
input_files = []

# If paired-end, this is a somewhat clean way to do it.
lee_lab_dir = '/Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/'
r1_suffix = '_L001_R1_001.fastq.gz'
r2_suffix = '_L001_R2_001.fastq.gz'

for i in (file_prefixes):
    input_files += [(lee_lab_dir + i + r1_suffix, lee_lab_dir + i + r2_suffix)]
    
    
# Testing to make sure valid input
assert len(paired_end) == len(run_fastqc)
assert len(paired_end) == len(has_adapters)
assert len(paired_end) == len(file_prefixes)
assert len(paired_end) == len(input_files)
```

## Static Variables


```python
# Adapters - same as the ones used in HITS-CLIP.
RL5 = 'AGGGAGGACGATGCGG' # 5' adapter
RL3 = 'GTGTCAGTCACTTCCAGCGG' # 3' adapter
```


```python
eber1_pos = (6629, 6795)

# getSeqByCoords.py -s 6629 -e 6795 -f EBV_Mod_B98-5 -c "AJ507799.2 Human herpesvirus 4 complete wild type genome"
eber1_seq = 'AGGACCTACGCTGCCCTAGAGGTTTTGCTAGGGAGGAGACGTGTGTGGCTGTAGCCACCCGTCCCGGGTACAAGTCCCGGGTGGTGAGGACGGTGTCTGTGGTTGTCTTCCCAGACTCTGCTTTCTGCCGTCTTCGGTCAAGTACCAGCTGGTGGTCCGCATGTTTT'
c_sites_eber1 = [x for x in range(0, len(eber1_seq)) if eber1_seq[x] == 'C']
```

## Load Python Libraries


```python
import os
import matplotlib.pyplot as plt
import numpy as np
import csv

from matplotlib.patches import Rectangle
from Bio.Seq import Seq
```


```python
!mkdir -p $output_dir
```

# Pre-process Data

## Get total Raw Read Counts


```python
total_raw_reads = []
for i in range(0, len(input_files)):
    in_file = input_files[i]
    if (paired_end[i]):
        in_file = in_file[0]

    raw_read_count = 0
    if in_file[-2:] == 'gz':
        raw_read_count = !gzcat $in_file | wc -l
        raw_read_count = raw_read_count[0]
    else:
        raw_read_count = !wc -l $in_file
        raw_read_count = raw_read_count[0].split()[0]
    
    if ('.fasta' in in_file or '.fa.' in in_file or in_file[-2:] == 'fa'):
        total_raw_reads.append(int(int(raw_read_count) / 2))
    else:
        total_raw_reads.append(int(int(raw_read_count) / 4))
        
print(total_raw_reads)
```

    [294330, 216369, 241087, 284886, 217289, 236392, 227075]


## FastQC


```python
def runFastqc(idx, fastqc_folder):
    
    # Get the input file(s) using the parameter index
    in_file = input_files[idx]
    # If it's paired-end, convert it to a string
    print(in_file)
    if (paired_end[idx]):
        in_file = in_file[0] + ' ' + in_file[1]
    print(in_file)
    prefix = file_prefixes[idx]
    
    # Run fastqc
    !fastqc -o $fastqc_folder --noextract $in_file
```


```python
fastqc_out = output_dir + '/fastqc'
!mkdir -p $fastqc_out
for i in range(0, len(input_files)):
    if (run_fastqc[i]):
        runFastqc(i, fastqc_out)
```


### View the FastQC files


```python
!open $output_dir"/fastqc/"*.html
```

### Clean up non-html files


```python
!rm $output_dir"/fastqc/"*.zip
```

## PEAR


```python
def pear(r1, r2, out_file):
    !pear -f $r1 -r $r2 -o $out_file
```


```python
# All files after the PEAR process will be added to this list,
# even if they don't get paired.
combined_files = []

!mkdir -p $output_dir"/pear"

for i in range(0, len(input_files)):
    if (paired_end[i]):
        r1 = input_files[i][0]
        r2 = input_files[i][1]
        
        # Get a prefix for the PEAR files
        out_file_name = (output_dir + '/pear/' + file_prefixes[i] + '_paired')
        
        # Run PEAR
        pear(r1, r2, out_file_name)
        
        # PEAR puts all the assembled reads into a ".assembled" file
        paired_file_name = out_file_name + '.assembled.fastq'
        combined_files.append(paired_file_name)
        
    else:
        # If not PEAR'd, just add the file to the new list so we can continue
        # processing with the new list.
        combined_files.append(input_files[i])
```

     ____  _____    _    ____ 
    |  _ \| ____|  / \  |  _ \
    | |_) |  _|   / _ \ | |_) |
    |  __/| |___ / ___ \|  _ <
    |_|   |_____/_/   \_\_| \_\
    
    PEAR v0.9.11 [Nov 5, 2017]
    
    Citation - PEAR: a fast and accurate Illumina Paired-End reAd mergeR
    Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593
    
    Forward reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/Vitro-1_S2_L001_R1_001.fastq.gz
    Reverse reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/Vitro-1_S2_L001_R2_001.fastq.gz
    PHRED..............................: 33
    Using empirical frequencies........: YES
    Statistical method.................: OES
    Maximum assembly length............: 999999
    Minimum assembly length............: 50
    p-value............................: 0.010000
    Quality score threshold (trimming).: 0
    Minimum read size after trimming...: 1
    Maximal ratio of uncalled bases....: 1.000000
    Minimum overlap....................: 10
    Scoring method.....................: Scaled score
    Threads............................: 1
    
    Allocating memory..................: 200,000,000 bytes
    Computing empirical frequencies....: DONE
      A: 0.244157
      C: 0.232684
      G: 0.292497
      T: 0.230661
      0 uncalled bases
    Assemblying reads: 100%
    
    Assembled reads ...................: 278,139 / 294,330 (94.499%)
    Discarded reads ...................: 0 / 294,330 (0.000%)
    Not assembled reads ...............: 16,191 / 294,330 (5.501%)
    Assembled reads file...............: Testing1/pear/Vitro-1_S2_paired.assembled.fastq
    Discarded reads file...............: Testing1/pear/Vitro-1_S2_paired.discarded.fastq
    Unassembled forward reads file.....: Testing1/pear/Vitro-1_S2_paired.unassembled.forward.fastq
    Unassembled reverse reads file.....: Testing1/pear/Vitro-1_S2_paired.unassembled.reverse.fastq
     ____  _____    _    ____ 
    |  _ \| ____|  / \  |  _ \
    | |_) |  _|   / _ \ | |_) |
    |  __/| |___ / ___ \|  _ <
    |_|   |_____/_/   \_\_| \_\
    
    PEAR v0.9.11 [Nov 5, 2017]
    
    Citation - PEAR: a fast and accurate Illumina Paired-End reAd mergeR
    Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593
    
    Forward reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/Vitro-2_S4_L001_R1_001.fastq.gz
    Reverse reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/Vitro-2_S4_L001_R2_001.fastq.gz
    PHRED..............................: 33
    Using empirical frequencies........: YES
    Statistical method.................: OES
    Maximum assembly length............: 999999
    Minimum assembly length............: 50
    p-value............................: 0.010000
    Quality score threshold (trimming).: 0
    Minimum read size after trimming...: 1
    Maximal ratio of uncalled bases....: 1.000000
    Minimum overlap....................: 10
    Scoring method.....................: Scaled score
    Threads............................: 1
    
    Allocating memory..................: 200,000,000 bytes
    Computing empirical frequencies....: DONE
      A: 0.310145
      C: 0.184589
      G: 0.224564
      T: 0.280701
      0 uncalled bases
    Assemblying reads: 100%
    
    Assembled reads ...................: 213,945 / 216,369 (98.880%)
    Discarded reads ...................: 0 / 216,369 (0.000%)
    Not assembled reads ...............: 2,424 / 216,369 (1.120%)
    Assembled reads file...............: Testing1/pear/Vitro-2_S4_paired.assembled.fastq
    Discarded reads file...............: Testing1/pear/Vitro-2_S4_paired.discarded.fastq
    Unassembled forward reads file.....: Testing1/pear/Vitro-2_S4_paired.unassembled.forward.fastq
    Unassembled reverse reads file.....: Testing1/pear/Vitro-2_S4_paired.unassembled.reverse.fastq
     ____  _____    _    ____ 
    |  _ \| ____|  / \  |  _ \
    | |_) |  _|   / _ \ | |_) |
    |  __/| |___ / ___ \|  _ <
    |_|   |_____/_/   \_\_| \_\
    
    PEAR v0.9.11 [Nov 5, 2017]
    
    Citation - PEAR: a fast and accurate Illumina Paired-End reAd mergeR
    Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593
    
    Forward reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/Vivo-1_S1_L001_R1_001.fastq.gz
    Reverse reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/Vivo-1_S1_L001_R2_001.fastq.gz
    PHRED..............................: 33
    Using empirical frequencies........: YES
    Statistical method.................: OES
    Maximum assembly length............: 999999
    Minimum assembly length............: 50
    p-value............................: 0.010000
    Quality score threshold (trimming).: 0
    Minimum read size after trimming...: 1
    Maximal ratio of uncalled bases....: 1.000000
    Minimum overlap....................: 10
    Scoring method.....................: Scaled score
    Threads............................: 1
    
    Allocating memory..................: 200,000,000 bytes
    Computing empirical frequencies....: DONE
      A: 0.316036
      C: 0.176882
      G: 0.213004
      T: 0.294078
      0 uncalled bases
    Assemblying reads: 100%
    
    Assembled reads ...................: 238,909 / 241,087 (99.097%)
    Discarded reads ...................: 0 / 241,087 (0.000%)
    Not assembled reads ...............: 2,178 / 241,087 (0.903%)
    Assembled reads file...............: Testing1/pear/Vivo-1_S1_paired.assembled.fastq
    Discarded reads file...............: Testing1/pear/Vivo-1_S1_paired.discarded.fastq
    Unassembled forward reads file.....: Testing1/pear/Vivo-1_S1_paired.unassembled.forward.fastq
    Unassembled reverse reads file.....: Testing1/pear/Vivo-1_S1_paired.unassembled.reverse.fastq
     ____  _____    _    ____ 
    |  _ \| ____|  / \  |  _ \
    | |_) |  _|   / _ \ | |_) |
    |  __/| |___ / ___ \|  _ <
    |_|   |_____/_/   \_\_| \_\
    
    PEAR v0.9.11 [Nov 5, 2017]
    
    Citation - PEAR: a fast and accurate Illumina Paired-End reAd mergeR
    Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593
    
    Forward reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/Vivo-2_S3_L001_R1_001.fastq.gz
    Reverse reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/Vivo-2_S3_L001_R2_001.fastq.gz
    PHRED..............................: 33
    Using empirical frequencies........: YES
    Statistical method.................: OES
    Maximum assembly length............: 999999
    Minimum assembly length............: 50
    p-value............................: 0.010000
    Quality score threshold (trimming).: 0
    Minimum read size after trimming...: 1
    Maximal ratio of uncalled bases....: 1.000000
    Minimum overlap....................: 10
    Scoring method.....................: Scaled score
    Threads............................: 1
    
    Allocating memory..................: 200,000,000 bytes
    Computing empirical frequencies....: DONE
      A: 0.316625
      C: 0.176464
      G: 0.212449
      T: 0.294462
      0 uncalled bases
    Assemblying reads: 100%
    
    Assembled reads ...................: 282,861 / 284,886 (99.289%)
    Discarded reads ...................: 0 / 284,886 (0.000%)
    Not assembled reads ...............: 2,025 / 284,886 (0.711%)
    Assembled reads file...............: Testing1/pear/Vivo-2_S3_paired.assembled.fastq
    Discarded reads file...............: Testing1/pear/Vivo-2_S3_paired.discarded.fastq
    Unassembled forward reads file.....: Testing1/pear/Vivo-2_S3_paired.unassembled.forward.fastq
    Unassembled reverse reads file.....: Testing1/pear/Vivo-2_S3_paired.unassembled.reverse.fastq
     ____  _____    _    ____ 
    |  _ \| ____|  / \  |  _ \
    | |_) |  _|   / _ \ | |_) |
    |  __/| |___ / ___ \|  _ <
    |_|   |_____/_/   \_\_| \_\
    
    PEAR v0.9.11 [Nov 5, 2017]
    
    Citation - PEAR: a fast and accurate Illumina Paired-End reAd mergeR
    Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593
    
    Forward reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/CTRL_S5_L001_R1_001.fastq.gz
    Reverse reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/CTRL_S5_L001_R2_001.fastq.gz
    PHRED..............................: 33
    Using empirical frequencies........: YES
    Statistical method.................: OES
    Maximum assembly length............: 999999
    Minimum assembly length............: 50
    p-value............................: 0.010000
    Quality score threshold (trimming).: 0
    Minimum read size after trimming...: 1
    Maximal ratio of uncalled bases....: 1.000000
    Minimum overlap....................: 10
    Scoring method.....................: Scaled score
    Threads............................: 1
    
    Allocating memory..................: 200,000,000 bytes
    Computing empirical frequencies....: DONE
      A: 0.326329
      C: 0.178733
      G: 0.178535
      T: 0.316404
      0 uncalled bases
    Assemblying reads: 100%
    
    Assembled reads ...................: 216,551 / 217,289 (99.660%)
    Discarded reads ...................: 0 / 217,289 (0.000%)
    Not assembled reads ...............: 738 / 217,289 (0.340%)
    Assembled reads file...............: Testing1/pear/CTRL_S5_paired.assembled.fastq
    Discarded reads file...............: Testing1/pear/CTRL_S5_paired.discarded.fastq
    Unassembled forward reads file.....: Testing1/pear/CTRL_S5_paired.unassembled.forward.fastq
    Unassembled reverse reads file.....: Testing1/pear/CTRL_S5_paired.unassembled.reverse.fastq
     ____  _____    _    ____ 
    |  _ \| ____|  / \  |  _ \
    | |_) |  _|   / _ \ | |_) |
    |  __/| |___ / ___ \|  _ <
    |_|   |_____/_/   \_\_| \_\
    
    PEAR v0.9.11 [Nov 5, 2017]
    
    Citation - PEAR: a fast and accurate Illumina Paired-End reAd mergeR
    Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593
    
    Forward reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/NSUN2_S6_L001_R1_001.fastq.gz
    Reverse reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/NSUN2_S6_L001_R2_001.fastq.gz
    PHRED..............................: 33
    Using empirical frequencies........: YES
    Statistical method.................: OES
    Maximum assembly length............: 999999
    Minimum assembly length............: 50
    p-value............................: 0.010000
    Quality score threshold (trimming).: 0
    Minimum read size after trimming...: 1
    Maximal ratio of uncalled bases....: 1.000000
    Minimum overlap....................: 10
    Scoring method.....................: Scaled score
    Threads............................: 1
    
    Allocating memory..................: 200,000,000 bytes
    Computing empirical frequencies....: DONE
      A: 0.325675
      C: 0.177391
      G: 0.183544
      T: 0.313391
      0 uncalled bases
    Assemblying reads: 100%
    
    Assembled reads ...................: 234,812 / 236,392 (99.332%)
    Discarded reads ...................: 0 / 236,392 (0.000%)
    Not assembled reads ...............: 1,580 / 236,392 (0.668%)
    Assembled reads file...............: Testing1/pear/NSUN2_S6_paired.assembled.fastq
    Discarded reads file...............: Testing1/pear/NSUN2_S6_paired.discarded.fastq
    Unassembled forward reads file.....: Testing1/pear/NSUN2_S6_paired.unassembled.forward.fastq
    Unassembled reverse reads file.....: Testing1/pear/NSUN2_S6_paired.unassembled.reverse.fastq
     ____  _____    _    ____ 
    |  _ \| ____|  / \  |  _ \
    | |_) |  _|   / _ \ | |_) |
    |  __/| |___ / ___ \|  _ <
    |_|   |_____/_/   \_\_| \_\
    
    PEAR v0.9.11 [Nov 5, 2017]
    
    Citation - PEAR: a fast and accurate Illumina Paired-End reAd mergeR
    Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593
    
    Forward reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/DMNT2_S7_L001_R1_001.fastq.gz
    Reverse reads file.................: /Volumes/Lee_Lab/Deep_Sequencing_Files/EBV_Bisulphite_20181012/DMNT2_S7_L001_R2_001.fastq.gz
    PHRED..............................: 33
    Using empirical frequencies........: YES
    Statistical method.................: OES
    Maximum assembly length............: 999999
    Minimum assembly length............: 50
    p-value............................: 0.010000
    Quality score threshold (trimming).: 0
    Minimum read size after trimming...: 1
    Maximal ratio of uncalled bases....: 1.000000
    Minimum overlap....................: 10
    Scoring method.....................: Scaled score
    Threads............................: 1
    
    Allocating memory..................: 200,000,000 bytes
    Computing empirical frequencies....: DONE
      A: 0.328972
      C: 0.175439
      G: 0.174914
      T: 0.320676
      0 uncalled bases
    Assemblying reads: 100%
    
    Assembled reads ...................: 225,836 / 227,075 (99.454%)
    Discarded reads ...................: 0 / 227,075 (0.000%)
    Not assembled reads ...............: 1,239 / 227,075 (0.546%)
    Assembled reads file...............: Testing1/pear/DMNT2_S7_paired.assembled.fastq
    Discarded reads file...............: Testing1/pear/DMNT2_S7_paired.discarded.fastq
    Unassembled forward reads file.....: Testing1/pear/DMNT2_S7_paired.unassembled.forward.fastq
    Unassembled reverse reads file.....: Testing1/pear/DMNT2_S7_paired.unassembled.reverse.fastq


### Remove non-assembled files


```python
!rm $output_dir"/pear/"*.unassembled.* $output_dir"/pear/"*.discarded.*
```

### Get Paired Read Counts


```python
total_peared_reads = []
for i in range(0, len(combined_files)):
    in_file = combined_files[i]
    pear_read_count = 0
    if in_file[-2:] == 'gz':
        pear_read_count = !gzcat $in_file | wc -l
        pear_read_count = pear_read_count[0]
    else:
        pear_read_count = !wc -l $in_file
        pear_read_count = pear_read_count[0].split()[0]
        
    if ('.fasta' in in_file or '.fa.' in in_file or in_file[-2:] == 'fa'):
        total_peared_reads.append(int(int(pear_read_count) / 2))
    else:
        total_peared_reads.append(int(int(pear_read_count) / 4))
        
print(total_peared_reads)
```

    [278139, 213945, 238909, 282861, 216551, 234812, 225836]


## Cutadapt

Remove HITS-CLIP adapters with cutadapt. It will only retains reads with either:
- RL5 ... RL3
- Reverse Complement of RL3 ... Reverse Complement of RL5


It will remove the adapters from those reads and combine them all into a single file for additional processing.


```python
def runCutAdapt(adp5, adp3, in_file, out_trimmed, out_untrimmed, err=0, min_overlap=0):
    e = ''
    O = ''
    if (err != 0):
        e = '-e ' + str(err) + ' '
    if (min_overlap != 0):
        O = '-O ' + str(min_overlap) + ' '
        
    # Anything less than this after removing adapters really isn't worth trying to map
    minimum_length = 15
    
    print('Running cutadapt on \"' + in_file + '\", outputting to \"' +
          out_trimmed + '\"')
    !cutadapt -m $minimum_length --no-indels -g $adp5"..."$adp3 $e $O --untrimmed-o $out_untrimmed -o $out_trimmed $in_file

```


```python
# All files after the cutadapt process will be added to this list,
# even if they don't get trimmed.
trimmed_files = []

!mkdir -p $output_dir"/cutadapt"

revRL3 = ''.join(Seq(RL3).reverse_complement())
revRL5 = ''.join(Seq(RL5).reverse_complement())

for i in range(0, len(combined_files)):
    if (has_adapters[i]):

        # Get an output prefix for the cutadapt files
        out_file_base = output_dir + '/cutadapt/' + file_prefixes[i]
        
        # Get the file extension
        file_ext = 'fastq'
        if ('.fasta' in combined_files[i] or '.fa.' in combined_files[i] or
            combined_files[i][-2:] == 'fa'):
            file_ext = 'fasta'
        
        # Cut in 5'-->3' Direction
        out_trimmed_1 = out_file_base + '_trimmed_1.' + file_ext
        out_untrimmed_1 = out_file_base + '_untrimmed_1.' + file_ext
        runCutAdapt(RL5, RL3, combined_files[i], out_trimmed_1, out_untrimmed_1)
        
        # Cut in the 3'-->5' Direction using the reverse complements
        out_trimmed_2 = out_file_base + '_trimmed_2.' + file_ext
        out_untrimmed_2 = out_file_base + '_untrimmed_2.' + file_ext
        runCutAdapt(revRL3, revRL5, out_untrimmed_1, out_trimmed_2, out_untrimmed_2)
        
        # The final file will be call "mySample_trimmed.fastq"
        out_trimmed_final = out_file_base + '_trimmed.' + file_ext
        !cat $out_trimmed_1 $out_trimmed_2 > $out_trimmed_final
        
        # Remove intermediary files
        !rm $out_trimmed_1 $out_trimmed_2 $out_untrimmed_1 $out_untrimmed_2
        
        trimmed_files.append(out_trimmed_final)
        
    else:
        # If no adapters, just add the file to the new list so we can continue
        # with alignment using the new list
        trimmed_files.append(combined_files[i])
```

    Running cutadapt on "Testing1/pear/Vitro-1_S2_paired.assembled.fastq", outputting to "Testing1/cutadapt/Vitro-1_S2_trimmed_1.fastq"
    This is cutadapt 1.18 with Python 3.7.0
    Command line parameters: -m 15 --no-indels -g AGGGAGGACGATGCGG...GTGTCAGTCACTTCCAGCGG --untrimmed-o Testing1/cutadapt/Vitro-1_S2_untrimmed_1.fastq -o Testing1/cutadapt/Vitro-1_S2_trimmed_1.fastq Testing1/pear/Vitro-1_S2_paired.assembled.fastq
    Processing reads on 1 core in single-end mode ...
    Finished in 3.57 s (13 us/read; 4.68 M reads/minute).
    
    === Summary ===
    
    Total reads processed:                 278,139
    Reads with adapters:                    72,767 (26.2%)
    Reads that were too short:              13,261 (4.8%)
    Reads written (passing filters):        59,506 (21.4%)
    
    Total basepairs processed:    48,013,831 bp
    Total written (filtered):      4,129,129 bp (8.6%)
    
    === Adapter 2 ===
    
    Sequence: AGGGAGGACGATGCGG...GTGTCAGTCACTTCCAGCGG; Type: linked; Length: 16+20; 5' trimmed: 72767 times; 3' trimmed: 72767 times
    
    No. of allowed errors:
    0-9 bp: 0; 10-16 bp: 1
    
    No. of allowed errors:
    0-9 bp: 0; 10-19 bp: 1; 20 bp: 2
 
    Running cutadapt on "Testing1/cutadapt/Vitro-1_S2_untrimmed_1.fastq", outputting to "Testing1/cutadapt/Vitro-1_S2_trimmed_2.fastq"
    This is cutadapt 1.18 with Python 3.7.0
    Command line parameters: -m 15 --no-indels -g CCGCTGGAAGTGACTGACAC...CCGCATCGTCCTCCCT --untrimmed-o Testing1/cutadapt/Vitro-1_S2_untrimmed_2.fastq -o Testing1/cutadapt/Vitro-1_S2_trimmed_2.fastq Testing1/cutadapt/Vitro-1_S2_untrimmed_1.fastq
    Processing reads on 1 core in single-end mode ...
    Finished in 3.49 s (17 us/read; 3.53 M reads/minute).
    
    === Summary ===
    
    Total reads processed:                 205,372
    Reads with adapters:                   192,897 (93.9%)
    Reads that were too short:              10,549 (5.1%)
    Reads written (passing filters):       182,348 (88.8%)
    
    Total basepairs processed:    36,372,962 bp
    Total written (filtered):     22,811,346 bp (62.7%)
    
    === Adapter 2 ===
    
    Sequence: CCGCTGGAAGTGACTGACAC...CCGCATCGTCCTCCCT; Type: linked; Length: 20+16; 5' trimmed: 192897 times; 3' trimmed: 192897 times
    
    No. of allowed errors:
    0-9 bp: 0; 10-19 bp: 1; 20 bp: 2
    
    No. of allowed errors:
    0-9 bp: 0; 10-16 bp: 1
    
    Running cutadapt on "Testing1/pear/Vitro-2_S4_paired.assembled.fastq", outputting to "Testing1/cutadapt/Vitro-2_S4_trimmed_1.fastq"
    This is cutadapt 1.18 with Python 3.7.0
    Command line parameters: -m 15 --no-indels -g AGGGAGGACGATGCGG...GTGTCAGTCACTTCCAGCGG --untrimmed-o Testing1/cutadapt/Vitro-2_S4_untrimmed_1.fastq -o Testing1/cutadapt/Vitro-2_S4_trimmed_1.fastq Testing1/pear/Vitro-2_S4_paired.assembled.fastq
    Processing reads on 1 core in single-end mode ...
    Finished in 2.87 s (13 us/read; 4.47 M reads/minute).
    
    === Summary ===
    
    Total reads processed:                 213,945
    Reads with adapters:                    39,578 (18.5%)
    Reads that were too short:                 348 (0.2%)
    Reads written (passing filters):        39,230 (18.3%)
    
    Total basepairs processed:    41,320,831 bp
    Total written (filtered):      5,375,777 bp (13.0%)
    
    === Adapter 2 ===
    
    Sequence: AGGGAGGACGATGCGG...GTGTCAGTCACTTCCAGCGG; Type: linked; Length: 16+20; 5' trimmed: 39578 times; 3' trimmed: 39578 times
    
    No. of allowed errors:
    0-9 bp: 0; 10-16 bp: 1
    
    No. of allowed errors:
    0-9 bp: 0; 10-19 bp: 1; 20 bp: 2
    
    Running cutadapt on "Testing1/cutadapt/Vitro-2_S4_untrimmed_1.fastq", outputting to "Testing1/cutadapt/Vitro-2_S4_trimmed_2.fastq"
    This is cutadapt 1.18 with Python 3.7.0
    Command line parameters: -m 15 --no-indels -g CCGCTGGAAGTGACTGACAC...CCGCATCGTCCTCCCT --untrimmed-o Testing1/cutadapt/Vitro-2_S4_untrimmed_2.fastq -o Testing1/cutadapt/Vitro-2_S4_trimmed_2.fastq Testing1/cutadapt/Vitro-2_S4_untrimmed_1.fastq
    Processing reads on 1 core in single-end mode ...
    Finished in 3.04 s (17 us/read; 3.44 M reads/minute).
    
    === Summary ===
    
    Total reads processed:                 174,367
    Reads with adapters:                   168,275 (96.5%)
    Reads that were too short:                 208 (0.1%)
    Reads written (passing filters):       168,067 (96.4%)
    
    Total basepairs processed:    34,451,193 bp
    Total written (filtered):     27,221,517 bp (79.0%)
    
    === Adapter 2 ===
    
    Sequence: CCGCTGGAAGTGACTGACAC...CCGCATCGTCCTCCCT; Type: linked; Length: 20+16; 5' trimmed: 168275 times; 3' trimmed: 168275 times
    
    No. of allowed errors:
    0-9 bp: 0; 10-19 bp: 1; 20 bp: 2
    
    No. of allowed errors:
    0-9 bp: 0; 10-16 bp: 1
    
    Running cutadapt on "Testing1/pear/Vivo-1_S1_paired.assembled.fastq", outputting to "Testing1/cutadapt/Vivo-1_S1_trimmed_1.fastq"
    This is cutadapt 1.18 with Python 3.7.0
    Command line parameters: -m 15 --no-indels -g AGGGAGGACGATGCGG...GTGTCAGTCACTTCCAGCGG --untrimmed-o Testing1/cutadapt/Vivo-1_S1_untrimmed_1.fastq -o Testing1/cutadapt/Vivo-1_S1_trimmed_1.fastq Testing1/pear/Vivo-1_S1_paired.assembled.fastq
    Processing reads on 1 core in single-end mode ...
    Finished in 2.87 s (12 us/read; 5.00 M reads/minute).
    
    === Summary ===
    
    Total reads processed:                 238,909
    Reads with adapters:                    28,590 (12.0%)
    Reads that were too short:                  26 (0.0%)
    Reads written (passing filters):        28,564 (12.0%)
    
    Total basepairs processed:    47,510,456 bp
    Total written (filtered):      4,306,641 bp (9.1%)
    
    === Adapter 2 ===
    
    Sequence: AGGGAGGACGATGCGG...GTGTCAGTCACTTCCAGCGG; Type: linked; Length: 16+20; 5' trimmed: 28590 times; 3' trimmed: 28590 times
    
    No. of allowed errors:
    0-9 bp: 0; 10-16 bp: 1
    
    No. of allowed errors:
    0-9 bp: 0; 10-19 bp: 1; 20 bp: 2
    
    Running cutadapt on "Testing1/cutadapt/Vivo-1_S1_untrimmed_1.fastq", outputting to "Testing1/cutadapt/Vivo-1_S1_trimmed_2.fastq"
    This is cutadapt 1.18 with Python 3.7.0
    Command line parameters: -m 15 --no-indels -g CCGCTGGAAGTGACTGACAC...CCGCATCGTCCTCCCT --untrimmed-o Testing1/cutadapt/Vivo-1_S1_untrimmed_2.fastq -o Testing1/cutadapt/Vivo-1_S1_trimmed_2.fastq Testing1/cutadapt/Vivo-1_S1_untrimmed_1.fastq
    Processing reads on 1 core in single-end mode ...
    Finished in 3.52 s (17 us/read; 3.59 M reads/minute).
    
    === Summary ===
    
    Total reads processed:                 210,319
    Reads with adapters:                   201,798 (95.9%)
    Reads that were too short:                  55 (0.0%)
    Reads written (passing filters):       201,743 (95.9%)
    
    Total basepairs processed:    42,167,308 bp
    Total written (filtered):     33,357,622 bp (79.1%)
    
    === Adapter 2 ===
    
    Sequence: CCGCTGGAAGTGACTGACAC...CCGCATCGTCCTCCCT; Type: linked; Length: 20+16; 5' trimmed: 201798 times; 3' trimmed: 201798 times
    
    No. of allowed errors:
    0-9 bp: 0; 10-19 bp: 1; 20 bp: 2
    
    No. of allowed errors:
    0-9 bp: 0; 10-16 bp: 1
    
    Running cutadapt on "Testing1/pear/Vivo-2_S3_paired.assembled.fastq", outputting to "Testing1/cutadapt/Vivo-2_S3_trimmed_1.fastq"
    This is cutadapt 1.18 with Python 3.7.0
    Command line parameters: -m 15 --no-indels -g AGGGAGGACGATGCGG...GTGTCAGTCACTTCCAGCGG --untrimmed-o Testing1/cutadapt/Vivo-2_S3_untrimmed_1.fastq -o Testing1/cutadapt/Vivo-2_S3_trimmed_1.fastq Testing1/pear/Vivo-2_S3_paired.assembled.fastq
    Processing reads on 1 core in single-end mode ...
    Finished in 3.59 s (13 us/read; 4.72 M reads/minute).
    
    === Summary ===
    
    Total reads processed:                 282,861
    Reads with adapters:                    57,702 (20.4%)
    Reads that were too short:                  26 (0.0%)
    Reads written (passing filters):        57,676 (20.4%)
    
    Total basepairs processed:    57,094,634 bp
    Total written (filtered):      9,574,511 bp (16.8%)
    
    === Adapter 2 ===
    
    Sequence: AGGGAGGACGATGCGG...GTGTCAGTCACTTCCAGCGG; Type: linked; Length: 16+20; 5' trimmed: 57702 times; 3' trimmed: 57702 times
    
    No. of allowed errors:
    0-9 bp: 0; 10-16 bp: 1
    
    No. of allowed errors:
    0-9 bp: 0; 10-19 bp: 1; 20 bp: 2
    
    Running cutadapt on "Testing1/cutadapt/Vivo-2_S3_untrimmed_1.fastq", outputting to "Testing1/cutadapt/Vivo-2_S3_trimmed_2.fastq"
    This is cutadapt 1.18 with Python 3.7.0
    Command line parameters: -m 15 --no-indels -g CCGCTGGAAGTGACTGACAC...CCGCATCGTCCTCCCT --untrimmed-o Testing1/cutadapt/Vivo-2_S3_untrimmed_2.fastq -o Testing1/cutadapt/Vivo-2_S3_trimmed_2.fastq Testing1/cutadapt/Vivo-2_S3_untrimmed_1.fastq
    Processing reads on 1 core in single-end mode ...
    Finished in 3.90 s (17 us/read; 3.46 M reads/minute).
    
    === Summary ===
    
    Total reads processed:                 225,159
    Reads with adapters:                   219,036 (97.3%)
    Reads that were too short:                  31 (0.0%)
    Reads written (passing filters):       219,005 (97.3%)
    
    Total basepairs processed:    45,437,367 bp
    Total written (filtered):     36,411,028 bp (80.1%)
    
    === Adapter 2 ===
    
    Sequence: CCGCTGGAAGTGACTGACAC...CCGCATCGTCCTCCCT; Type: linked; Length: 20+16; 5' trimmed: 219036 times; 3' trimmed: 219036 times
    
    No. of allowed errors:
    0-9 bp: 0; 10-19 bp: 1; 20 bp: 2
    
    No. of allowed errors:
    0-9 bp: 0; 10-16 bp: 1
    

### Get Trimmed Read Counts


```python
total_trimmed_reads = []
for i in range(0, len(trimmed_files)):
    in_file = trimmed_files[i]
    trimmed_read_count = 0
    if in_file[-2:] == 'gz':
        trimmed_read_count = !gzcat $in_file | wc -l
        trimmed_read_count = trimmed_read_count[0]
    else:
        trimmed_read_count = !wc -l $in_file
        trimmed_read_count = trimmed_read_count[0].split()[0]
        
    if ('.fasta' in in_file or '.fa.' in in_file or in_file[-2:] == 'fa'):
        total_trimmed_reads.append(int(int(trimmed_read_count) / 2))
    else:
        total_trimmed_reads.append(int(int(trimmed_read_count) / 4))
        
print(total_trimmed_reads)
```

    [241854, 207297, 230307, 276681, 216551, 234812, 225836]


# Align the Data with Bismark

Apparently, most of the bisulphite-aware aligners are pretty similar ([Tsuji, T., Weng, Z., 2016](https://academic.oup.com/bib/article-lookup/doi/10.1093/bib/bbv103)). Given that, and given its success in other experiments ([Wreczycka et al., 2017](https://www.sciencedirect.com/science/article/pii/S0168165617315936#bib0175)), I decided to go with *Bismark* ([Krueger et al., 2011](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btr167)).

Bismark v0.20.0 source code was downloaded from [here](https://github.com/FelixKrueger/Bismark/releases/tag/0.20.0). General pipeline protocol can be found [here](https://www.epigenesys.eu/en/protocols/bio-informatics/483-quality-control-trimming-and-alignment-of-bisulfite-seq-data-prot-57).

Per the protocol, we need to generate a bismark genome index. Since we only care about mapping to EBER1, and since we're actually mapping RNA instead of DNA (i.e. there shouldn't be any reads that extend beyond EBER1's genomic locus), I figured it would make sense to JUST map it against the EBER1 sequence. I therefore prepared a special genome file specifically for EBER1 (see Simulation_Generator notebook).

However, just using the EBER1 sequence was causing errors in Bismark, apparently because "a sequence aligns right to the vert[sic] end (or start for reverse sequences) of a reference genome, and typically this is on the MT. The error arises because Bismark would like to extract 2 bp further downstream of the sequence to determine the cytosine context but this is obviously not possible when there are no further 2bp available." ([Felix Krueger](https://github.com/FelixKrueger/Bismark/issues/68)).

Consequently, **I map to a subset of the EBV genome that goes up- and down-stream 20 bp from the EBER1 gene.**


```python
eber1_subset_genome_folder = '/Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset'
```


```python
sam_files = []

!mkdir -p $output_dir"/bismark"

for i in range(0, len(trimmed_files)):
    
    file_ext = 'fastq'
    if ('.fasta' in trimmed_files[i] or '.fa.' in trimmed_files[i] or
        trimmed_files[i][-2:] == 'fa'):
        file_ext = 'fasta'
    
    file_to_align = trimmed_files[i]
    prefix = file_prefixes[i]

    # -B = file prefix
    # -o = (relative) output folder
    # --non_directional specifies that it could be either same as genome
    #   or in the reverse direction
    # -N = Number of mismatches allowed in a "seed alignment". Lower = more sensitive.
    # -L = Length of "seed substrings". Lower = more sensitive, default is 20, max is 32.
    !bismark -B $prefix -o $output_dir"/bismark" --$file_ext --non_directional --bowtie2 -N 1 -L 10 $eber1_subset_genome_folder/ $file_to_align
    
    bam_file = output_dir + '/bismark/' + prefix + '.bam'
    sam_file = bam_file.replace('.bam', '.sam')
    
    # Create SAM file
    !samtools view $bam_file > $sam_file
    
    sam_files.append(sam_file)

```

    Path to Bowtie 2 specified as: bowtie2
    Bowtie seems to be working fine (tested command 'bowtie2 --version' [2.3.4])
    Output format is BAM (default)
    Alignments will be written out in BAM format. Samtools found here: '/usr/local/bin/samtools'
    Reference genome folder provided is /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/	(absolute path is '/Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/)'
    FastQ format specified
    
    Input files to be analysed (in current folder '/Users/jpk90/Desktop/BisulfiteSequencing/experiments/test'):
    Testing1/cutadapt/Vitro-1_S2_trimmed.fastq
    Library was specified to be not strand-specific (non-directional), therefore alignments to all four possible bisulfite strands (OT, CTOT, OB and CTOB) will be reported
    Output will be written into the directory: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/
    Setting parallelization to single-threaded (default)
    
    Current working directory is: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test
    
    Now reading in and storing sequence information of the genome specified in: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/
    
    chr EBV_Mod_B98-5_EBER1_Subset (207 bp)
    
    Single-core mode: setting pid to 1
    
    Single-end alignments will be performed
    =======================================
    
    Input file is in FastQ format
    Writing a C -> T converted version of the input file Vitro-1_S2_trimmed.fastq to Vitro-1_S2_trimmed.fastq_C_to_T.fastq
    Writing a G -> A converted version of the input file Vitro-1_S2_trimmed.fastq to Vitro-1_S2_trimmed.fastq_G_to_A.fastq
    
    Created C -> T as well as G -> A converted versions of the FastQ file Vitro-1_S2_trimmed.fastq (241854 sequences in total)
    
    Input files are Vitro-1_S2_trimmed.fastq_C_to_T.fastq and Vitro-1_S2_trimmed.fastq_G_to_A.fastq (FastQ)
    
    Now running 4 individual instances of Bowtie 2 against the bisulfite genome of /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/ with the specified options: -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals
    
    Now starting the Bowtie 2 aligner for CTreadCTgenome (reading in sequences from Vitro-1_S2_trimmed.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:2620:1000_1:N:0:2	4	*	0	0	*	*	0	0	TGTTAGTTATTTTTAGTGG	FFFFFFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for CTreadGAgenome (reading in sequences from Vitro-1_S2_trimmed.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:2620:1000_1:N:0:2	4	*	0	0	*	*	0	0	TGTTAGTTATTTTTAGTGG	FFFFFFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadCTgenome (reading in sequences from Vitro-1_S2_trimmed.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:2620:1000_1:N:0:2	4	*	0	0	*	*	0	0	TATCAATCACTTCCAACAA	FFFFFFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadGAgenome (reading in sequences from Vitro-1_S2_trimmed.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:2620:1000_1:N:0:2	4	*	0	0	*	*	0	0	TATCAATCACTTCCAACAA	FFFFFFFFFFFFFFFFFFF	YT:Z:UU
    
    >>> Writing bisulfite mapping results to /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/Vitro-1_S2.bam <<<
    
    
    Reading in the sequence file Testing1/cutadapt/Vitro-1_S2_trimmed.fastq
    241854 reads; of these:
      241854 (100.00%) were unpaired; of these:
        241851 (100.00%) aligned 0 times
        3 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.00% overall alignment rate
    241854 reads; of these:
      241854 (100.00%) were unpaired; of these:
        225013 (93.04%) aligned 0 times
        16841 (6.96%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    6.96% overall alignment rate
    241854 reads; of these:
      241854 (100.00%) were unpaired; of these:
        241846 (100.00%) aligned 0 times
        8 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.00% overall alignment rate
    241854 reads; of these:
      241854 (100.00%) were unpaired; of these:
        148286 (61.31%) aligned 0 times
        93568 (38.69%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    38.69% overall alignment rate
    Processed 241854 sequences in total
    
    
    Successfully deleted the temporary files Vitro-1_S2_trimmed.fastq_C_to_T.fastq and Vitro-1_S2_trimmed.fastq_G_to_A.fastq
    
    Final Alignment report
    ======================
    Sequences analysed in total:	241854
    Number of alignments with a unique best hit from the different alignments:	110412
    Mapping efficiency:	45.7%
    
    Sequences with no alignments under any condition:	131442
    Sequences did not map uniquely:	0
    Sequences which were discarded because genomic sequence could not be extracted:	0
    
    Number of sequences with unique best (first) alignment came from the bowtie output:
    CT/CT:	16834	((converted) top strand)
    CT/GA:	3	((converted) bottom strand)
    GA/CT:	93568	(complementary to (converted) top strand)
    GA/GA:	7	(complementary to (converted) bottom strand)
    
    Final Cytosine Methylation Report
    =================================
    Total number of C's analysed:	4611592
    
    Total methylated C's in CpG context:	67918
    Total methylated C's in CHG context:	84285
    Total methylated C's in CHH context:	134146
    Total methylated C's in Unknown context:	17
    
    Total unmethylated C's in CpG context:	903706
    Total unmethylated C's in CHG context:	1312194
    Total unmethylated C's in CHH context:	2109343
    Total unmethylated C's in Unknown context:	354
    
    C methylated in CpG context:	7.0%
    C methylated in CHG context:	6.0%
    C methylated in CHH context:	6.0%
    C methylated in Unknown context (CN or CHN):	4.6%
    
    
    Bismark completed in 0d 0h 0m 50s
    
    ====================
    Bismark run complete
    ====================
    
    Path to Bowtie 2 specified as: bowtie2
    Bowtie seems to be working fine (tested command 'bowtie2 --version' [2.3.4])
    Output format is BAM (default)
    Alignments will be written out in BAM format. Samtools found here: '/usr/local/bin/samtools'
    Reference genome folder provided is /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/	(absolute path is '/Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/)'
    FastQ format specified
    
    Input files to be analysed (in current folder '/Users/jpk90/Desktop/BisulfiteSequencing/experiments/test'):
    Testing1/cutadapt/Vitro-2_S4_trimmed.fastq
    Library was specified to be not strand-specific (non-directional), therefore alignments to all four possible bisulfite strands (OT, CTOT, OB and CTOB) will be reported
    Output will be written into the directory: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/
    Setting parallelization to single-threaded (default)
    
    Current working directory is: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test
    
    Now reading in and storing sequence information of the genome specified in: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/
    
    chr EBV_Mod_B98-5_EBER1_Subset (207 bp)
    
    Single-core mode: setting pid to 1
    
    Single-end alignments will be performed
    =======================================
    
    Input file is in FastQ format
    Writing a C -> T converted version of the input file Vitro-2_S4_trimmed.fastq to Vitro-2_S4_trimmed.fastq_C_to_T.fastq
    Writing a G -> A converted version of the input file Vitro-2_S4_trimmed.fastq to Vitro-2_S4_trimmed.fastq_G_to_A.fastq
    
    Created C -> T as well as G -> A converted versions of the FastQ file Vitro-2_S4_trimmed.fastq (207297 sequences in total)
    
    Input files are Vitro-2_S4_trimmed.fastq_C_to_T.fastq and Vitro-2_S4_trimmed.fastq_G_to_A.fastq (FastQ)
    
    Now running 4 individual instances of Bowtie 2 against the bisulfite genome of /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/ with the specified options: -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals
    
    Now starting the Bowtie 2 aligner for CTreadCTgenome (reading in sequences from Vitro-2_S4_trimmed.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:2300:1000_1:N:0:4	0	EBV_Mod_B98-5_EBER1_Subset_CT_converted	18	0	23M1D146M	*	0	0	GGGAGGATTTATGTTGTTTTAGAGTTTTGTTAGGGAGGAGATGTGTGTGGTTGTAGTTATTTGTTTTGGGTATAAGTTTTGGGTGGTGAGGATGGTGTTTGTGGTTGTTTTTTTAGATTTTGTTTTTTGTTGTTTTTGGTTAAGTATTAGTTGGTGGTTTGTATGTTTT	FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF	AS:i:-26	XN:i:0	XM:i:3	XO:i:1	XG:i:1	NM:i:4	MD:Z:0T0T0T20^G146	YT:Z:UU
    Now starting the Bowtie 2 aligner for CTreadGAgenome (reading in sequences from Vitro-2_S4_trimmed.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:2300:1000_1:N:0:4	4	*	0	0	*	*	0	0	GGGAGGATTTATGTTGTTTTAGAGTTTTGTTAGGGAGGAGATGTGTGTGGTTGTAGTTATTTGTTTTGGGTATAAGTTTTGGGTGGTGAGGATGGTGTTTGTGGTTGTTTTTTTAGATTTTGTTTTTTGTTGTTTTTGGTTAAGTATTAGTTGGTGGTTTGTATGTTTT	FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadCTgenome (reading in sequences from Vitro-2_S4_trimmed.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:2300:1000_1:N:0:4	4	*	0	0	*	*	0	0	AAAAAAATTTATATTATTTTAAAATTTTATTAAAAAAAAAATATATATAATTATAATTATTTATCCCAAATATAAATCCCAAATAATAAAAACAATATTTATAATTATTTTTTTAAATTTTATTTTTTATTATCTTTAATTAAATATTAATTAATAATTTATATATTTT	FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadGAgenome (reading in sequences from Vitro-2_S4_trimmed.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:2300:1000_1:N:0:4	4	*	0	0	*	*	0	0	AAAAAAATTTATATTATTTTAAAATTTTATTAAAAAAAAAATATATATAATTATAATTATTTATCCCAAATATAAATCCCAAATAATAAAAACAATATTTATAATTATTTTTTTAAATTTTATTTTTTATTATCTTTAATTAAATATTAATTAATAATTTATATATTTT	FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF	YT:Z:UU
    
    >>> Writing bisulfite mapping results to /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/Vitro-2_S4.bam <<<
    
    
    Reading in the sequence file Testing1/cutadapt/Vitro-2_S4_trimmed.fastq
    207297 reads; of these:
      207297 (100.00%) were unpaired; of these:
        192651 (92.93%) aligned 0 times
        14646 (7.07%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    7.07% overall alignment rate
    207297 reads; of these:
      207297 (100.00%) were unpaired; of these:
        207292 (100.00%) aligned 0 times
        5 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.00% overall alignment rate
    207297 reads; of these:
      207297 (100.00%) were unpaired; of these:
        207288 (100.00%) aligned 0 times
        9 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.00% overall alignment rate
    207297 reads; of these:
      207297 (100.00%) were unpaired; of these:
        115203 (55.57%) aligned 0 times
        92094 (44.43%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    44.43% overall alignment rate
    Processed 207297 sequences in total
    
    
    Successfully deleted the temporary files Vitro-2_S4_trimmed.fastq_C_to_T.fastq and Vitro-2_S4_trimmed.fastq_G_to_A.fastq
    
    Final Alignment report
    ======================
    Sequences analysed in total:	207297
    Number of alignments with a unique best hit from the different alignments:	106740
    Mapping efficiency:	51.5%
    
    Sequences with no alignments under any condition:	100557
    Sequences did not map uniquely:	0
    Sequences which were discarded because genomic sequence could not be extracted:	0
    
    Number of sequences with unique best (first) alignment came from the bowtie output:
    CT/CT:	14646	((converted) top strand)
    CT/GA:	4	((converted) bottom strand)
    GA/CT:	92090	(complementary to (converted) top strand)
    GA/GA:	0	(complementary to (converted) bottom strand)
    
    Final Cytosine Methylation Report
    =================================
    Total number of C's analysed:	4493862
    
    Total methylated C's in CpG context:	87260
    Total methylated C's in CHG context:	104997
    Total methylated C's in CHH context:	153551
    Total methylated C's in Unknown context:	31
    
    Total unmethylated C's in CpG context:	858225
    Total unmethylated C's in CHG context:	1247617
    Total unmethylated C's in CHH context:	2042212
    Total unmethylated C's in Unknown context:	619
    
    C methylated in CpG context:	9.2%
    C methylated in CHG context:	7.8%
    C methylated in CHH context:	7.0%
    C methylated in Unknown context (CN or CHN):	4.8%
    
    
    Bismark completed in 0d 0h 0m 51s
    
    ====================
    Bismark run complete
    ====================
    
    Path to Bowtie 2 specified as: bowtie2
    Bowtie seems to be working fine (tested command 'bowtie2 --version' [2.3.4])
    Output format is BAM (default)
    Alignments will be written out in BAM format. Samtools found here: '/usr/local/bin/samtools'
    Reference genome folder provided is /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/	(absolute path is '/Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/)'
    FastQ format specified
    
    Input files to be analysed (in current folder '/Users/jpk90/Desktop/BisulfiteSequencing/experiments/test'):
    Testing1/cutadapt/Vivo-1_S1_trimmed.fastq
    Library was specified to be not strand-specific (non-directional), therefore alignments to all four possible bisulfite strands (OT, CTOT, OB and CTOB) will be reported
    Output will be written into the directory: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/
    Setting parallelization to single-threaded (default)
    
    Current working directory is: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test
    
    Now reading in and storing sequence information of the genome specified in: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/
    
    chr EBV_Mod_B98-5_EBER1_Subset (207 bp)
    
    Single-core mode: setting pid to 1
    
    Single-end alignments will be performed
    =======================================
    
    Input file is in FastQ format
    Writing a C -> T converted version of the input file Vivo-1_S1_trimmed.fastq to Vivo-1_S1_trimmed.fastq_C_to_T.fastq
    Writing a G -> A converted version of the input file Vivo-1_S1_trimmed.fastq to Vivo-1_S1_trimmed.fastq_G_to_A.fastq
    
    Created C -> T as well as G -> A converted versions of the FastQ file Vivo-1_S1_trimmed.fastq (230307 sequences in total)
    
    Input files are Vivo-1_S1_trimmed.fastq_C_to_T.fastq and Vivo-1_S1_trimmed.fastq_G_to_A.fastq (FastQ)
    
    Now running 4 individual instances of Bowtie 2 against the bisulfite genome of /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/ with the specified options: -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals
    
    Now starting the Bowtie 2 aligner for CTreadCTgenome (reading in sequences from Vivo-1_S1_trimmed.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:5650:1000_1:N:0:1	0	EBV_Mod_B98-5_EBER1_Subset_CT_converted	21	42	167M	*	0	0	AGGATTTATGTTGTTTTAGAGGTTTTGTTAGGGAGGAGATGTGTGTGGTTGTAGTTATTTGTTTTGGGTATAAGTTTTGGGTGGTGAGGATGGTGTTTGTGGTTGTTTTTTTAGATTTTGTTTTTTGTTGTTTTTGGTTAAGTATTAGTTGGTGGTTTGTATGTTTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:167	YT:Z:UU
    Now starting the Bowtie 2 aligner for CTreadGAgenome (reading in sequences from Vivo-1_S1_trimmed.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:5650:1000_1:N:0:1	4	*	0	0	*	*	0	0	AGGATTTATGTTGTTTTAGAGGTTTTGTTAGGGAGGAGATGTGTGTGGTTGTAGTTATTTGTTTTGGGTATAAGTTTTGGGTGGTGAGGATGGTGTTTGTGGTTGTTTTTTTAGATTTTGTTTTTTGTTGTTTTTGGTTAAGTATTAGTTGGTGGTTTGTATGTTTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadCTgenome (reading in sequences from Vivo-1_S1_trimmed.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:5650:1000_1:N:0:1	4	*	0	0	*	*	0	0	AAAATTTATATTATTTTAAAAATTTTATTAAAAAAAAAATATATATAATTATAATTATTTATTTTAAATATAAATTTTAAATAATAAAAATAATATTTATAATTATTTTTTTAAATTTTATTTTTTATTATTTTTAATTAAATACCAATTAATAATTTATATATTTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadGAgenome (reading in sequences from Vivo-1_S1_trimmed.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:5650:1000_1:N:0:1	4	*	0	0	*	*	0	0	AAAATTTATATTATTTTAAAAATTTTATTAAAAAAAAAATATATATAATTATAATTATTTATTTTAAATATAAATTTTAAATAATAAAAATAATATTTATAATTATTTTTTTAAATTTTATTTTTTATTATTTTTAATTAAATACCAATTAATAATTTATATATTTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	YT:Z:UU
    
    >>> Writing bisulfite mapping results to /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/Vivo-1_S1.bam <<<
    
    
    Reading in the sequence file Testing1/cutadapt/Vivo-1_S1_trimmed.fastq
    230307 reads; of these:
      230307 (100.00%) were unpaired; of these:
        202816 (88.06%) aligned 0 times
        27491 (11.94%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    11.94% overall alignment rate
    230307 reads; of these:
      230307 (100.00%) were unpaired; of these:
    230307     reads; of these:230307 (100.00%
    ) aligned 0 times
          2303070 ( (0.00%100.00) aligned exactly 1 time%
    ) were unpaired; of these:    
    0     (2303060.00 (%) aligned >1 times100.00
    %) aligned 0 times0.00
    %     overall alignment rate1
     (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.00% overall alignment rate
    230307 reads; of these:
      230307 (100.00%) were unpaired; of these:
        32740 (14.22%) aligned 0 times
        197567 (85.78%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    85.78% overall alignment rate
    Processed 230307 sequences in total
    
    
    Successfully deleted the temporary files Vivo-1_S1_trimmed.fastq_C_to_T.fastq and Vivo-1_S1_trimmed.fastq_G_to_A.fastq
    
    Final Alignment report
    ======================
    Sequences analysed in total:	230307
    Number of alignments with a unique best hit from the different alignments:	225058
    Mapping efficiency:	97.7%
    
    Sequences with no alignments under any condition:	5249
    Sequences did not map uniquely:	0
    Sequences which were discarded because genomic sequence could not be extracted:	0
    
    Number of sequences with unique best (first) alignment came from the bowtie output:
    CT/CT:	27490	((converted) top strand)
    CT/GA:	0	((converted) bottom strand)
    GA/CT:	197567	(complementary to (converted) top strand)
    GA/GA:	1	(complementary to (converted) bottom strand)
    
    Final Cytosine Methylation Report
    =================================
    Total number of C's analysed:	9490995
    
    Total methylated C's in CpG context:	157566
    Total methylated C's in CHG context:	216360
    Total methylated C's in CHH context:	454984
    Total methylated C's in Unknown context:	47
    
    Total unmethylated C's in CpG context:	1838802
    Total unmethylated C's in CHG context:	2653216
    Total unmethylated C's in CHH context:	4170067
    Total unmethylated C's in Unknown context:	794
    
    C methylated in CpG context:	7.9%
    C methylated in CHG context:	7.5%
    C methylated in CHH context:	9.8%
    C methylated in Unknown context (CN or CHN):	5.6%
    
    
    Bismark completed in 0d 0h 1m 29s
    
    ====================
    Bismark run complete
    ====================
    
    Path to Bowtie 2 specified as: bowtie2
    Bowtie seems to be working fine (tested command 'bowtie2 --version' [2.3.4])
    Output format is BAM (default)
    Alignments will be written out in BAM format. Samtools found here: '/usr/local/bin/samtools'
    Reference genome folder provided is /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/	(absolute path is '/Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/)'
    FastQ format specified
    
    Input files to be analysed (in current folder '/Users/jpk90/Desktop/BisulfiteSequencing/experiments/test'):
    Testing1/cutadapt/Vivo-2_S3_trimmed.fastq
    Library was specified to be not strand-specific (non-directional), therefore alignments to all four possible bisulfite strands (OT, CTOT, OB and CTOB) will be reported
    Output will be written into the directory: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/
    Setting parallelization to single-threaded (default)
    
    Current working directory is: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test
    
    Now reading in and storing sequence information of the genome specified in: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/
    
    chr EBV_Mod_B98-5_EBER1_Subset (207 bp)
    
    Single-core mode: setting pid to 1
    
    Single-end alignments will be performed
    =======================================
    
    Input file is in FastQ format
    Writing a C -> T converted version of the input file Vivo-2_S3_trimmed.fastq to Vivo-2_S3_trimmed.fastq_C_to_T.fastq
    Writing a G -> A converted version of the input file Vivo-2_S3_trimmed.fastq to Vivo-2_S3_trimmed.fastq_G_to_A.fastq
    
    Created C -> T as well as G -> A converted versions of the FastQ file Vivo-2_S3_trimmed.fastq (276681 sequences in total)
    
    Input files are Vivo-2_S3_trimmed.fastq_C_to_T.fastq and Vivo-2_S3_trimmed.fastq_G_to_A.fastq (FastQ)
    
    Now running 4 individual instances of Bowtie 2 against the bisulfite genome of /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/ with the specified options: -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals
    
    Now starting the Bowtie 2 aligner for CTreadCTgenome (reading in sequences from Vivo-2_S3_trimmed.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1340:1000_1:N:0:3	0	EBV_Mod_B98-5_EBER1_Subset_CT_converted	18	8	170M	*	0	0	GGGAGGATTTATGTTGTTTTAGAGGTTTTGTTAGGGAGGAGATGTGTGTGGTTGTAGTTATTTGTTTTGGGTATAAGTTTTGGGTGGTGAGGATGGTGTTTGTGGTTGTTTTTTTAGATTTTGTTTTTTGTTGTTTTTGGTTAAGTATTAGTTGGTGGTTTGTATGTTTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	AS:i:-18	XN:i:0	XM:i:3	XO:i:0	XG:i:0	NM:i:3	MD:Z:0T0T0T167	YT:Z:UU
    Now starting the Bowtie 2 aligner for CTreadGAgenome (reading in sequences from Vivo-2_S3_trimmed.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1340:1000_1:N:0:3	4	*	0	0	*	*	0	0	GGGAGGATTTATGTTGTTTTAGAGGTTTTGTTAGGGAGGAGATGTGTGTGGTTGTAGTTATTTGTTTTGGGTATAAGTTTTGGGTGGTGAGGATGGTGTTTGTGGTTGTTTTTTTAGATTTTGTTTTTTGTTGTTTTTGGTTAAGTATTAGTTGGTGGTTTGTATGTTTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadCTgenome (reading in sequences from Vivo-2_S3_trimmed.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1340:1000_1:N:0:3	4	*	0	0	*	*	0	0	AAAAAAATTTATATTATTTTAAAAATTTTATTAAAAAAAAAATATATATAATTATAATTATTTATCTTAAATATAAATCCCAAATAATAAAAATAATATCTATAATTATTTTTTTAAATTTTATTTTTTATTATTTTTAATTAAATATTAATTAATAATTTATATATTTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadGAgenome (reading in sequences from Vivo-2_S3_trimmed.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1340:1000_1:N:0:3	4	*	0	0	*	*	0	0	AAAAAAATTTATATTATTTTAAAAATTTTATTAAAAAAAAAATATATATAATTATAATTATTTATCTTAAATATAAATCCCAAATAATAAAAATAATATCTATAATTATTTTTTTAAATTTTATTTTTTATTATTTTTAATTAAATATTAATTAATAATTTATATATTTT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	YT:Z:UU
    
    >>> Writing bisulfite mapping results to /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/Vivo-2_S3.bam <<<
    
    
    Reading in the sequence file Testing1/cutadapt/Vivo-2_S3_trimmed.fastq
    276681 reads; of these:
      276681 (276681 reads; of these:
    100.00  %) were unpaired; of these:
        276681 (100.00%) aligned 0 times
        2766810 ( (0.00%) aligned exactly 1 time
    100.00    %0) were unpaired; of these: (
    0.00    %276681) aligned >1 times (
    100.000.00%%) aligned 0 times overall alignment rate
    
        0 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.00% overall alignment rate
    276681 reads; of these:
      276681 (100.00%) were unpaired; of these:
        220918 (79.85%) aligned 0 times
        55763 (20.15%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    20.15% overall alignment rate
    276681 reads; of these:
      276681 (100.00%) were unpaired; of these:
        62597 (22.62%) aligned 0 times
        214084 (77.38%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    77.38% overall alignment rate
    Processed 276681 sequences in total
    
    
    Successfully deleted the temporary files Vivo-2_S3_trimmed.fastq_C_to_T.fastq and Vivo-2_S3_trimmed.fastq_G_to_A.fastq
    
    Final Alignment report
    ======================
    Sequences analysed in total:	276681
    Number of alignments with a unique best hit from the different alignments:	269847
    Mapping efficiency:	97.5%
    
    Sequences with no alignments under any condition:	6834
    Sequences did not map uniquely:	0
    Sequences which were discarded because genomic sequence could not be extracted:	0
    
    Number of sequences with unique best (first) alignment came from the bowtie output:
    CT/CT:	55763	((converted) top strand)
    CT/GA:	0	((converted) bottom strand)
    GA/CT:	214084	(complementary to (converted) top strand)
    GA/GA:	0	(complementary to (converted) bottom strand)
    
    Final Cytosine Methylation Report
    =================================
    Total number of C's analysed:	11556648
    
    Total methylated C's in CpG context:	208130
    Total methylated C's in CHG context:	277461
    Total methylated C's in CHH context:	581560
    Total methylated C's in Unknown context:	4
    
    Total unmethylated C's in CpG context:	2211715
    Total unmethylated C's in CHG context:	3217889
    Total unmethylated C's in CHH context:	5059893
    Total unmethylated C's in Unknown context:	707
    
    C methylated in CpG context:	8.6%
    C methylated in CHG context:	7.9%
    C methylated in CHH context:	10.3%
    C methylated in Unknown context (CN or CHN):	0.6%
    
    
    Bismark completed in 0d 0h 1m 53s
    
    ====================
    Bismark run complete
    ====================
    
    Path to Bowtie 2 specified as: bowtie2
    Bowtie seems to be working fine (tested command 'bowtie2 --version' [2.3.4])
    Output format is BAM (default)
    Alignments will be written out in BAM format. Samtools found here: '/usr/local/bin/samtools'
    Reference genome folder provided is /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/	(absolute path is '/Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/)'
    FastQ format specified
    
    Input files to be analysed (in current folder '/Users/jpk90/Desktop/BisulfiteSequencing/experiments/test'):
    Testing1/pear/CTRL_S5_paired.assembled.fastq
    Library was specified to be not strand-specific (non-directional), therefore alignments to all four possible bisulfite strands (OT, CTOT, OB and CTOB) will be reported
    Output will be written into the directory: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/
    Setting parallelization to single-threaded (default)
    
    Current working directory is: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test
    
    Now reading in and storing sequence information of the genome specified in: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/
    
    chr EBV_Mod_B98-5_EBER1_Subset (207 bp)
    
    Single-core mode: setting pid to 1
    
    Single-end alignments will be performed
    =======================================
    
    Input file is in FastQ format
    Writing a C -> T converted version of the input file CTRL_S5_paired.assembled.fastq to CTRL_S5_paired.assembled.fastq_C_to_T.fastq
    Writing a G -> A converted version of the input file CTRL_S5_paired.assembled.fastq to CTRL_S5_paired.assembled.fastq_G_to_A.fastq
    
    Created C -> T as well as G -> A converted versions of the FastQ file CTRL_S5_paired.assembled.fastq (216551 sequences in total)
    
    Input files are CTRL_S5_paired.assembled.fastq_C_to_T.fastq and CTRL_S5_paired.assembled.fastq_G_to_A.fastq (FastQ)
    
    Now running 4 individual instances of Bowtie 2 against the bisulfite genome of /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/ with the specified options: -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals
    
    Now starting the Bowtie 2 aligner for CTreadCTgenome (reading in sequences from CTRL_S5_paired.assembled.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1710:1000_1:N:0:5	0	EBV_Mod_B98-5_EBER1_Subset_CT_converted	21	0	74M1D24M1D30M2D35M	*	0	0	AGGATTTATGTTGTTTTAGAGGTTTTGTTAGGGAGGAGATGTGTGTGGTTGTAGTTATTTGTTTTGGGTATAAGTTTGGGTGGTGAGGATGGTGTTTGGGTTGTTTTTTTAGATTTTGTTTTTTGTTGTTTGGTTAAGTATTAGTTGGTGGTTTGTATGTTTT	FFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFF	AS:i:-27	XN:i:0	XM:i:0	XO:i:3	XG:i:4	NM:i:4	MD:Z:74^T24^T30^TT35	YT:Z:UU
    Now starting the Bowtie 2 aligner for CTreadGAgenome (reading in sequences from CTRL_S5_paired.assembled.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1710:1000_1:N:0:5	4	*	0	0	*	*	0	0	AGGATTTATGTTGTTTTAGAGGTTTTGTTAGGGAGGAGATGTGTGTGGTTGTAGTTATTTGTTTTGGGTATAAGTTTGGGTGGTGAGGATGGTGTTTGGGTTGTTTTTTTAGATTTTGTTTTTTGTTGTTTGGTTAAGTATTAGTTGGTGGTTTGTATGTTTT	FFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadCTgenome (reading in sequences from CTRL_S5_paired.assembled.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1710:1000_1:N:0:5	4	*	0	0	*	*	0	0	AAAATTTATATTATTTTAAAAATTTTATTAAAAAAAAAATATATATAATTATAATTATTTATTTTAAATATAAATTTAAATAATAAAAATAATATTTAAATTATTTTTTTAAATTTTATTTTTTATTATTTAATTAAATACTAATTAATAATTTATATATTTT	FFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadGAgenome (reading in sequences from CTRL_S5_paired.assembled.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1710:1000_1:N:0:5	4	*	0	0	*	*	0	0	AAAATTTATATTATTTTAAAAATTTTATTAAAAAAAAAATATATATAATTATAATTATTTATTTTAAATATAAATTTAAATAATAAAAATAATATTTAAATTATTTTTTTAAATTTTATTTTTTATTATTTAATTAAATACTAATTAATAATTTATATATTTT	FFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFF	YT:Z:UU
    
    >>> Writing bisulfite mapping results to /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/CTRL_S5.bam <<<
    
    
    Reading in the sequence file Testing1/pear/CTRL_S5_paired.assembled.fastq
    216551 reads; of these:
      216551 (100.00%) were unpaired; of these:
        216551 (100.00%) aligned 0 times
        0 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.00% overall alignment rate
    216551 reads; of these:
      216551 (100.00%) were unpaired; of these:
        216551 (100.00%) aligned 0 times
        0 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.00% overall alignment rate
    216551 reads; of these:
      216551 (100.00%) were unpaired; of these:
        167557 (77.38%) aligned 0 times
        48994 (22.62%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    22.62% overall alignment rate
    216551 reads; of these:
      216551 (100.00%) were unpaired; of these:
        72784 (33.61%) aligned 0 times
        143767 (66.39%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    66.39% overall alignment rate
    Processed 216551 sequences in total
    
    
    Successfully deleted the temporary files CTRL_S5_paired.assembled.fastq_C_to_T.fastq and CTRL_S5_paired.assembled.fastq_G_to_A.fastq
    
    Final Alignment report
    ======================
    Sequences analysed in total:	216551
    Number of alignments with a unique best hit from the different alignments:	192761
    Mapping efficiency:	89.0%
    
    Sequences with no alignments under any condition:	23790
    Sequences did not map uniquely:	0
    Sequences which were discarded because genomic sequence could not be extracted:	0
    
    Number of sequences with unique best (first) alignment came from the bowtie output:
    CT/CT:	48994	((converted) top strand)
    CT/GA:	0	((converted) bottom strand)
    GA/CT:	143767	(complementary to (converted) top strand)
    GA/GA:	0	(complementary to (converted) bottom strand)
    
    Final Cytosine Methylation Report
    =================================
    Total number of C's analysed:	8214826
    
    Total methylated C's in CpG context:	59179
    Total methylated C's in CHG context:	265627
    Total methylated C's in CHH context:	283559
    Total methylated C's in Unknown context:	12
    
    Total unmethylated C's in CpG context:	1662983
    Total unmethylated C's in CHG context:	2221819
    Total unmethylated C's in CHH context:	3721659
    Total unmethylated C's in Unknown context:	2030
    
    C methylated in CpG context:	3.4%
    C methylated in CHG context:	10.7%
    C methylated in CHH context:	7.1%
    C methylated in Unknown context (CN or CHN):	0.6%
    
    
    Bismark completed in 0d 0h 1m 20s
    
    ====================
    Bismark run complete
    ====================
    
    Path to Bowtie 2 specified as: bowtie2
    Bowtie seems to be working fine (tested command 'bowtie2 --version' [2.3.4])
    Output format is BAM (default)
    Alignments will be written out in BAM format. Samtools found here: '/usr/local/bin/samtools'
    Reference genome folder provided is /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/	(absolute path is '/Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/)'
    FastQ format specified
    
    Input files to be analysed (in current folder '/Users/jpk90/Desktop/BisulfiteSequencing/experiments/test'):
    Testing1/pear/NSUN2_S6_paired.assembled.fastq
    Library was specified to be not strand-specific (non-directional), therefore alignments to all four possible bisulfite strands (OT, CTOT, OB and CTOB) will be reported
    Output will be written into the directory: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/
    Setting parallelization to single-threaded (default)
    
    Current working directory is: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test
    
    Now reading in and storing sequence information of the genome specified in: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/
    
    chr EBV_Mod_B98-5_EBER1_Subset (207 bp)
    
    Single-core mode: setting pid to 1
    
    Single-end alignments will be performed
    =======================================
    
    Input file is in FastQ format
    Writing a C -> T converted version of the input file NSUN2_S6_paired.assembled.fastq to NSUN2_S6_paired.assembled.fastq_C_to_T.fastq
    Writing a G -> A converted version of the input file NSUN2_S6_paired.assembled.fastq to NSUN2_S6_paired.assembled.fastq_G_to_A.fastq
    
    Created C -> T as well as G -> A converted versions of the FastQ file NSUN2_S6_paired.assembled.fastq (234812 sequences in total)
    
    Input files are NSUN2_S6_paired.assembled.fastq_C_to_T.fastq and NSUN2_S6_paired.assembled.fastq_G_to_A.fastq (FastQ)
    
    Now running 4 individual instances of Bowtie 2 against the bisulfite genome of /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/ with the specified options: -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals
    
    Now starting the Bowtie 2 aligner for CTreadCTgenome (reading in sequences from NSUN2_S6_paired.assembled.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1170:1000_1:N:0:6	4	*	0	0	*	*	0	0	AAAATATATAAATTATTAATTAATATTTATTAAAAATAATAAAAAATAAAATTTAAAAAAATAATTATAAATATTATTTTTTTATTTAAAATTTATATTTGAAATGGGTAATTATAGTTATATATATTTTTTTTTTAATAAAATTTTTAAAATAATATAAATTTT	FFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for CTreadGAgenome (reading in sequences from NSUN2_S6_paired.assembled.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1170:1000_1:N:0:6	4	*	0	0	*	*	0	0	AAAATATATAAATTATTAATTAATATTTATTAAAAATAATAAAAAATAAAATTTAAAAAAATAATTATAAATATTATTTTTTTATTTAAAATTTATATTTGAAATGGGTAATTATAGTTATATATATTTTTTTTTTAATAAAATTTTTAAAATAATATAAATTTT	FFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadCTgenome (reading in sequences from NSUN2_S6_paired.assembled.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1170:1000_1:N:0:6	16	EBV_Mod_B98-5_EBER1_Subset_CT_converted	21	23	84M1D52M1D29M	*	0	0	AGGATTTATGTTGTTTTAGAGGTTTTGTTAGGGAGGAGATGTGTGTGGTTGTAGTTATTTGTTTTGGGTATAAGTTTTGGGTGGGAGGATGGTGTTTGTGGTTGTTTTTTTAGATTTTGTTTTTTGTTGTTTTTGGTAAGTATTAGTTGGTGGTTTGTATGTTTT	FFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFF	AS:i:-16	XN:i:0	XM:i:0	XO:i:2	XG:i:2	NM:i:2	MD:Z:84^T52^T29	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadGAgenome (reading in sequences from NSUN2_S6_paired.assembled.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1170:1000_1:N:0:6	4	*	0	0	*	*	0	0	AAAACATACAAACCACCAACTAATACTTACCAAAAACAACAAAAAACAAAATCTAAAAAAACAACCACAAACACCATCCTCCCACCCAAAACTTATACCCAAAACAAATAACTACAACCACACACATCTCCTCCCTAACAAAACCTCTAAAACAACATAAATCCT	FFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFF	YT:Z:UU
    
    >>> Writing bisulfite mapping results to /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/NSUN2_S6.bam <<<
    
    
    Reading in the sequence file Testing1/pear/NSUN2_S6_paired.assembled.fastq
    234812 reads; of these:
      234812 (100.00%) were unpaired; of these:
        234812 (100.00%) aligned 0 times
        0 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.00% overall alignment rate
    234812 reads; of these:
      234812 (100.00%) were unpaired; of these:
        234812 (100.00%) aligned 0 times
        0 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.00% overall alignment rate
    234812 reads; of these:
      234812 (100.00%) were unpaired; of these:
        179953 (76.64%) aligned 0 times
        54859 (23.36%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    23.36% overall alignment rate
    234812 reads; of these:
      234812 (100.00%) were unpaired; of these:
        129899 (55.32%) aligned 0 times
        104913 (44.68%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    44.68% overall alignment rate
    Processed 234812 sequences in total
    
    
    Successfully deleted the temporary files NSUN2_S6_paired.assembled.fastq_C_to_T.fastq and NSUN2_S6_paired.assembled.fastq_G_to_A.fastq
    
    Final Alignment report
    ======================
    Sequences analysed in total:	234812
    Number of alignments with a unique best hit from the different alignments:	159772
    Mapping efficiency:	68.0%
    
    Sequences with no alignments under any condition:	75040
    Sequences did not map uniquely:	0
    Sequences which were discarded because genomic sequence could not be extracted:	0
    
    Number of sequences with unique best (first) alignment came from the bowtie output:
    CT/CT:	54859	((converted) top strand)
    CT/GA:	0	((converted) bottom strand)
    GA/CT:	104913	(complementary to (converted) top strand)
    GA/GA:	0	(complementary to (converted) bottom strand)
    
    Final Cytosine Methylation Report
    =================================
    Total number of C's analysed:	6674041
    
    Total methylated C's in CpG context:	90729
    Total methylated C's in CHG context:	127036
    Total methylated C's in CHH context:	269415
    Total methylated C's in Unknown context:	6
    
    Total unmethylated C's in CpG context:	1313705
    Total unmethylated C's in CHG context:	1903253
    Total unmethylated C's in CHH context:	2969903
    Total unmethylated C's in Unknown context:	3507
    
    C methylated in CpG context:	6.5%
    C methylated in CHG context:	6.3%
    C methylated in CHH context:	8.3%
    C methylated in Unknown context (CN or CHN):	0.2%
    
    
    Bismark completed in 0d 0h 1m 19s
    
    ====================
    Bismark run complete
    ====================
    
    Path to Bowtie 2 specified as: bowtie2
    Bowtie seems to be working fine (tested command 'bowtie2 --version' [2.3.4])
    Output format is BAM (default)
    Alignments will be written out in BAM format. Samtools found here: '/usr/local/bin/samtools'
    Reference genome folder provided is /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/	(absolute path is '/Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/)'
    FastQ format specified
    
    Input files to be analysed (in current folder '/Users/jpk90/Desktop/BisulfiteSequencing/experiments/test'):
    Testing1/pear/DMNT2_S7_paired.assembled.fastq
    Library was specified to be not strand-specific (non-directional), therefore alignments to all four possible bisulfite strands (OT, CTOT, OB and CTOB) will be reported
    Output will be written into the directory: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/
    Setting parallelization to single-threaded (default)
    
    Current working directory is: /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test
    
    Now reading in and storing sequence information of the genome specified in: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/
    
    chr EBV_Mod_B98-5_EBER1_Subset (207 bp)
    
    Single-core mode: setting pid to 1
    
    Single-end alignments will be performed
    =======================================
    
    Input file is in FastQ format
    Writing a C -> T converted version of the input file DMNT2_S7_paired.assembled.fastq to DMNT2_S7_paired.assembled.fastq_C_to_T.fastq
    Writing a G -> A converted version of the input file DMNT2_S7_paired.assembled.fastq to DMNT2_S7_paired.assembled.fastq_G_to_A.fastq
    
    Created C -> T as well as G -> A converted versions of the FastQ file DMNT2_S7_paired.assembled.fastq (225836 sequences in total)
    
    Input files are DMNT2_S7_paired.assembled.fastq_C_to_T.fastq and DMNT2_S7_paired.assembled.fastq_G_to_A.fastq (FastQ)
    
    Now running 4 individual instances of Bowtie 2 against the bisulfite genome of /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/ with the specified options: -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals
    
    Now starting the Bowtie 2 aligner for CTreadCTgenome (reading in sequences from DMNT2_S7_paired.assembled.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1220:1000_1:N:0:7	4	*	0	0	*	*	0	0	AAAATATATAAATTATTAATTGGTATTTAATTAAAAATAATAAAAAATAAAATTTAAAAAAATAATTATAAATATTATTTTTATTATTTAAAATTTATATTTAAAATAAATAATTATAATTATATATATTTTTTTTTTAATAAAATTTTTAAAATAATATAAATTTT	FFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for CTreadGAgenome (reading in sequences from DMNT2_S7_paired.assembled.fastq_C_to_T.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1220:1000_1:N:0:7	4	*	0	0	*	*	0	0	AAAATATATAAATTATTAATTGGTATTTAATTAAAAATAATAAAAAATAAAATTTAAAAAAATAATTATAAATATTATTTTTATTATTTAAAATTTATATTTAAAATAAATAATTATAATTATATATATTTTTTTTTTAATAAAATTTTTAAAATAATATAAATTTT	FFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFF	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadCTgenome (reading in sequences from DMNT2_S7_paired.assembled.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --nofw)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/CT_conversion/BS_CT
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1220:1000_1:N:0:7	16	EBV_Mod_B98-5_EBER1_Subset_CT_converted	21	42	167M	*	0	0	AGGATTTATGTTGTTTTAGAGGTTTTGTTAGGGAGGAGATGTGTGTGGTTGTAGTTATTTGTTTTGGGTATAAGTTTTGGGTGGTGAGGATGGTGTTTGTGGTTGTTTTTTTAGATTTTGTTTTTTGTTGTTTTTGGTTAAGTATTAGTTGGTGGTTTGTATGTTTT	FFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFF	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:167	YT:Z:UU
    Now starting the Bowtie 2 aligner for GAreadGAgenome (reading in sequences from DMNT2_S7_paired.assembled.fastq_G_to_A.fastq with options -q -N 1 -L 10 --score-min L,0,-0.2 --ignore-quals --norc)
    Using Bowtie 2 index: /Volumes/Lee_Lab/Jack_Files/Genomes/EBV/EBV_Mod_B98-5_EBER1_Subset/Bisulfite_Genome/GA_conversion/BS_GA
    
    Found first alignment:	FS10000252:5:BNT40308-1613:1:1101:1220:1000_1:N:0:7	4	*	0	0	*	*	0	0	AAAACATACAAACCACCAACTAATACTTAACCAAAAACAACAAAAAACAAAATCTAAAAAAACAACCACAAACACCATCCTCACCACCCAAAACTTATACCCAAAACAAATAACTACAACCACACACATCTCCTCCCTAACAAAACCTCTAAAACAACATAAATCCT	FFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFF	YT:Z:UU
    
    >>> Writing bisulfite mapping results to /Users/jpk90/Desktop/BisulfiteSequencing/experiments/test/Testing1/bismark/DMNT2_S7.bam <<<
    
    
    Reading in the sequence file Testing1/pear/DMNT2_S7_paired.assembled.fastq
    225836 reads; of these:225836 reads; of these:
      
    225836   (225836 (100.00%100.00) were unpaired; of these:%
    ) were unpaired; of these:    
    225836     (225836100.00 (%100.00) aligned 0 times%
    ) aligned 0 times    
    0     (00.00 (%0.00) aligned exactly 1 time%
    ) aligned exactly 1 time    
    0     (00.00 (%0.00) aligned >1 times%
    ) aligned >1 times0.00
    %0.00 overall alignment rate%
     overall alignment rate
    225836 reads; of these:
      225836 (100.00%) were unpaired; of these:
        154299 (68.32%) aligned 0 times
        71537 (31.68%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    31.68% overall alignment rate
    225836 reads; of these:
      225836 (100.00%) were unpaired; of these:
        89668 (39.70%) aligned 0 times
        136168 (60.30%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    60.30% overall alignment rate
    Processed 225836 sequences in total
    
    
    Successfully deleted the temporary files DMNT2_S7_paired.assembled.fastq_C_to_T.fastq and DMNT2_S7_paired.assembled.fastq_G_to_A.fastq
    
    Final Alignment report
    ======================
    Sequences analysed in total:	225836
    Number of alignments with a unique best hit from the different alignments:	207705
    Mapping efficiency:	92.0%
    
    Sequences with no alignments under any condition:	18131
    Sequences did not map uniquely:	0
    Sequences which were discarded because genomic sequence could not be extracted:	0
    
    Number of sequences with unique best (first) alignment came from the bowtie output:
    CT/CT:	71537	((converted) top strand)
    CT/GA:	0	((converted) bottom strand)
    GA/CT:	136168	(complementary to (converted) top strand)
    GA/GA:	0	(complementary to (converted) bottom strand)
    
    Final Cytosine Methylation Report
    =================================
    Total number of C's analysed:	8872929
    
    Total methylated C's in CpG context:	14425
    Total methylated C's in CHG context:	202400
    Total methylated C's in CHH context:	222971
    Total methylated C's in Unknown context:	11
    
    Total unmethylated C's in CpG context:	1844120
    Total unmethylated C's in CHG context:	2483132
    Total unmethylated C's in CHH context:	4105881
    Total unmethylated C's in Unknown context:	1645
    
    C methylated in CpG context:	0.8%
    C methylated in CHG context:	7.5%
    C methylated in CHH context:	5.2%
    C methylated in Unknown context (CN or CHN):	0.7%
    
    
    Bismark completed in 0d 0h 1m 20s
    
    ====================
    Bismark run complete
    ====================
    


### Get Aligned Read Counts


```python
total_aligned_reads = []
for i in range(0, len(sam_files)):
    in_file = sam_files[i]
    aligned_read_count = 0
    aligned_read_count = !wc -l $in_file
    aligned_read_count = aligned_read_count[0].split()[0]

    total_aligned_reads.append(int(aligned_read_count))
        
print(total_aligned_reads)
```

    [110412, 106740, 225058, 269847, 192761, 159772, 207705]


# Post-Alignment Analysis

## Get methylated positions from Bismark's XM:Z tag on the SAM file

Each line in the output file should have an "XM:Z:..." tag. Each position in this tag (after the "XM:Z:") corresponds to the nucleotide position mapped on the reference genome. If the nucleotide was A, G, or T on the reference genome, it appears as ".". If it was c and unmethyled, it's a lowercase letter; if it's C and methylated, it's an uppercase letter.

See the [Bismark User Guide](https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf) (p. 11-12) for more info.


```python
def getCPositions(in_sam_file, csv_file_path):
    num_recs  = 0
    s_file = open(in_sam_file, 'r')
    pos = []
    pos_dict = {}
    fixed_seq = ''
    for line in s_file:
        parts = line.strip().split()
        for part in parts:
            if 'XM:Z:' in part:
                methylation_seq = part[5:]
        if methylation_seq is not None:
            num_recs += 1
            #pos += [x for x in range(0, len(methylation_seq)) if methylation_seq[x] in ('ZHXU')]
            for i in range(0, len(methylation_seq)):
                if methylation_seq[i] in ('ZHXU'):
                    # Append i+1 for position, since string is 0-based but seq is 1-based
                    pos.append(i+1)
                    if (i+1) in pos_dict:
                        pos_dict[i+1] += 1
                    else:
                        pos_dict[i+1] = 1
    s_file.close()
    
    # Write to csv file
    csv_file = open(csv_file_path, 'w+')
    writer = csv.writer(csv_file)
    for key, value in pos_dict.items():
       writer.writerow([key, value])
    
    return [num_recs, pos, pos_dict]
```

## Heatmap Generation

Not sure how helpful this is, but figured I'd try it to see.


```python
def makeHeatMapFromDict(pos_dict):
    el = len(eber1_seq)
    box = np.zeros((el, el))
    for i in range(0, el):
        if (i in pos_dict):
            box[i,i] = pos_dict[i]
    
    plt.imshow(box, cmap='Greys', interpolation='nearest')
    plt.show()
```

## Histogram of Methylated C Positions


```python
def plotMethylLocations(title, num_recs, positions, img_out_path=''):
    plt.title(title)
    plt.ylim([0, num_recs])
    plt.xlabel('Nucleotide Position in EBER1 Sequence')
    plt.ylabel('Reads with Methylated C')
    plt.hist(positions, bins=range(0, len(eber1_seq)))
    #plt.hist(positions, bins=[x for x in range(0, len(eber1_seq)) if eber1_seq[x] == 'C'])
    if (img_out_path != ''):
        plt.savefig(img_out_path)
    plt.show()
```

## Run the Analysis


```python
hist_folder = output_dir + '/analysis/histograms'
csv_folder = output_dir + '/analysis/methylation_loci'

!mkdir -p $hist_folder
!mkdir -p $csv_folder

pos_dicts = []
for i in range(0, len(sam_files)):

    # Get Methylation loci
    csv_file_name = csv_folder + '/' + file_prefixes[i] + '.csv'
    num_recs, pos, pos_dict = getCPositions(sam_files[i], csv_file_name)
    
    print('Significantly methylated positions (at least 20% of all' +
          ' mapped reads retained a cytosine at this position):\n')
    print([str(x) + ': ' + str(round(pos_dict[x] / num_recs * 100, 2)) + '%' for x in pos_dict if (pos_dict[x] / num_recs) > 0.2])
    
    # Plot Histogram
    title = 'Methylated Nucleotides for ' + file_prefixes[i]
    histogram_file_name = hist_folder + '/' + file_prefixes[i] + '.png'
    plotMethylLocations(title, num_recs, pos, histogram_file_name)
    #makeHeatMapFromDict(pos_dict)
    
    pos_dicts += [pos_dict]
```

    Significantly methylated positions (at least 20% of all mapped reads retained a cytosine at this position):
    
    []



![png](/assets/20181012_bisulphite/output_49_1.png)


    Significantly methylated positions (at least 20% of all mapped reads retained a cytosine at this position):
    
    []



![png](/assets/20181012_bisulphite/output_49_3.png)


    Significantly methylated positions (at least 20% of all mapped reads retained a cytosine at this position):
    
    ['145: 62.68%']



![png](/assets/20181012_bisulphite/output_49_5.png)


    Significantly methylated positions (at least 20% of all mapped reads retained a cytosine at this position):
    
    ['145: 70.6%']



![png](/assets/20181012_bisulphite/output_49_7.png)


    Significantly methylated positions (at least 20% of all mapped reads retained a cytosine at this position):
    
    ['145: 62.86%', '49: 45.86%', '59: 48.04%']



![png](/assets/20181012_bisulphite/output_49_9.png)


    Significantly methylated positions (at least 20% of all mapped reads retained a cytosine at this position):
    
    []



![png](/assets/20181012_bisulphite/output_49_11.png)


    Significantly methylated positions (at least 20% of all mapped reads retained a cytosine at this position):
    
    ['145: 89.01%', '146: 82.08%']



![png](/assets/20181012_bisulphite/output_49_13.png)


# General Stats:

| Sample Name  | Raw Reads | Paired Reads | With Adapters | Aligned EBER1
| :----------  | --------- | ------------ | ------------- | -------------
| *in vivo* 1  | {{total_raw_reads[0]}} | {{total_peared_reads[0]}} | {{total_trimmed_reads[0]}} | {{total_aligned_reads[0]}}
| *in vivo* 2  | {{total_raw_reads[1]}} | {{total_peared_reads[1]}} | {{total_trimmed_reads[1]}} | {{total_aligned_reads[1]}}
| *in vitro* 1 | {{total_raw_reads[2]}} | {{total_peared_reads[2]}} | {{total_trimmed_reads[2]}} | {{total_aligned_reads[2]}}
| *in vitro* 2 | {{total_raw_reads[3]}} | {{total_peared_reads[3]}} | {{total_trimmed_reads[3]}} | {{total_aligned_reads[3]}}
| CTRL         | {{total_raw_reads[4]}} | {{total_peared_reads[4]}} | {{total_trimmed_reads[4]}} | {{total_aligned_reads[4]}}
| NSUN2-kd     | {{total_raw_reads[5]}} | {{total_peared_reads[5]}} | {{total_trimmed_reads[5]}} | {{total_aligned_reads[5]}}
| DMNT2-kd     | {{total_raw_reads[6]}} | {{total_peared_reads[6]}} | {{total_trimmed_reads[6]}} | {{total_aligned_reads[6]}}

## Method to make Overall .csv


```python
def generateFullCsv(pos_dicts, out_csv_file_path):
    
    # Get the maximum position in any of the position dictionaries
    # Should be the length of EBER1, but because of misalignments, it
    # can go over.
    max_pos = 1
    for i in pos_dicts:
        if max_pos < max(i):
            max_pos = max(i)
            
    # Populate the dictionaries with '0' in whatever positions they
    # don't have, up to the maximum position in any of the dictionaries
    for i in range(1, max_pos + 1):
        for j in range(0, len(pos_dicts)):
            if i not in pos_dicts[j]:
                pos_dicts[j][i] = 0
    
    # Create a sorted list of tuples from each dictionary
    pos_lists = []
    for i in pos_dicts:
        pos_lists += [sorted(i.items())]
        
    # Finally, start writing to a .csv file
    csv_file = open(out_csv_file_path, 'w+')
    writer = csv.writer(csv_file)
    
    # Print headers
    headers = ['Position', 'Nucleotide']
    for i in file_prefixes:
        headers += [str(i + ' Num')]
        headers += [str(i + ' %')]
    writer.writerow(headers)
    
    # Print data rows
    for i in range(0, len(pos_lists[0])):
        pos = pos_lists[0][i][0]
        nt = 'N'
        if (pos <= len(eber1_seq)):
            nt = eber1_seq[pos-1]
        output_row = [pos, nt]
        for j in range(0, len(pos_lists)):
            # Number
            num = pos_lists[j][i][1]
            output_row.append(num)
            # Percentage of Total Reads for Sample
            output_row.append(str(round(num / total_aligned_reads[j] * 100, 4)) + '%')
        writer.writerow(output_row)
```


```python
overall_csv = output_dir + '/analysis/overall.csv'
generateFullCsv(pos_dicts, overall_csv)
```

## Forna

Place the following in [forna](http://nibiru.tbi.univie.ac.at/forna/) to get a picture of the nucleotide on EBER1:

Under Add Molecule:

>\>EBER1
AGGACCTACGCTGCCCTAGAGGTTTTGCTAGGGAGGAGACGTGTGTGGCTGTAGCCACCCGTCCCGGGTACAAGTCCCGGGTGGTGAGGACGGTGTCTGTGGTTGTCTTCCCAGACTCTGCTTTCTGCCGTCTTCGGTCAAGTACCAGCTGGTGGTCCGCATGTTTT

Under Colors --> Custom:

>145:red

(based on whichever numbers were returned from the dictionary above). Also, make the figure look as close to the one in [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4378190/figure/F1/) as possible.

# 10/17 Update - New Analysis

Belle was curious as to the distribution of methylation marks per read. Let's find out...


```python
def getCPositions2(in_sam_file):
    num_recs  = 0
    s_file = open(in_sam_file, 'r')
    pos = []
    read_dict = {}
    fixed_seq = ''
    for line in s_file:
        pos_list = ''
        parts = line.strip().split()
        for part in parts:
            if 'XM:Z:' in part:
                methylation_seq = part[5:]
        if methylation_seq is not None:
            num_recs += 1
            #pos += [x for x in range(0, len(methylation_seq)) if methylation_seq[x] in ('ZHXU')]
            for i in range(0, len(methylation_seq)):
                if methylation_seq[i] in ('ZHXU'):
                    # Append i+1 for position, since string is 0-based but seq is 1-based
                    pos.append(i+1)
                    pos_list += str(i+1) + '_'
        if pos_list in read_dict:
            read_dict[pos_list] += 1
        else:
            read_dict[pos_list] = 1
    s_file.close()
    
    return read_dict
```


```python
!samtools view bismark/Vivo-1_S1_L001_paired_trimmed_bismark_bt2.bam > bismark/Vivo-1.sam
```


```python
x = getCPositions2('bismark/Vivo-1.sam')
```


```python
import operator
sorted_x = sorted(x.items(), key=operator.itemgetter(1), reverse=True)
print(sorted_x)
```


```python
def getNumForPos(in_sam_file, has_pos, doesnt_have_pos):
    '''
    positions should be a list, e.g. [144, 145, ..., 2]
    '''
    num_recs  = 0
    s_file = open(in_sam_file, 'r')
    sum_reads = 0
    for line in s_file:
        parts = line.strip().split()
        for part in parts:
            if 'XM:Z:' in part:
                methylation_seq = part[5:]
        if methylation_seq is not None:
            has_all_pos = True
            for i in has_pos:
                if ((i-1 >= len(methylation_seq)) or 
                    (methylation_seq[i-1] not in ('ZHXU'))):
                    has_all_pos = False
            for i in doesnt_have_pos:
                if ((i-1 < len(methylation_seq)) and 
                    (methylation_seq[i-1] in ('ZHXU'))):
                    has_all_pos = False
            if (has_all_pos):
                sum_reads += 1
    s_file.close()
    
    return sum_reads
```


```python
include = [59]
dont_include = [x for x in range(1, 160) if x != 59]

getNumForPos('bismark/Vivo-1.sam', include, dont_include)
```


