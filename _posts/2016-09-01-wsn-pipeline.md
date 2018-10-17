---
layout: post
title:  "WSN HITS-CLIP 2016-09-01"
date:   2018-10-17 15:07:23 -0400
categories: HITS-CLIP
---
# HITS-CLIP Pipeline

For information on the process, see [here](https://www.nature.com/articles/nprot.2014.012). The general idea behind the process is:

1. Cross-link proteins (e.g. NP) to RNA (e.g. influenza vRNA)
2. Isolate the RNA, add adapters, & PCR
3. Sequence the results to determine where RNA your protein of interest is binding to 

![protocol.jpg](attachment:protocol.jpg)

## Notes on implementation

1. Although we do paired-end sequencing, we treat them as single-end.
2. When removing adapters, **we only look for whether the read has a 5' or 3' adapter, not both (should this be the case??)**

# Setup

## User-Defined Variables

Change these before proceeding.


```python
# Note: do NOT include an end "/"
output_dir = 'testing_2'

# This will not only serve as a convenient means for indicating the input
# files, but will be used as file prefixes for output files
file_prefixes = ['WSN']

# Define separately for each sample
paired_end = [True]
run_fastqc = [True]
has_adapters = [True]

# Strain-specific files
chrom_file = '/Volumes/Lee_Lab/Jack_Files/Chrom_Size_Files/WSN.chrom.sizes'
genome_file = '/Volumes/Lee_Lab/Jack_Files/Genomes/Influenza/WSN/WSN.fasta'
```

STAR parameters


```python
num_threads = 1
star_idx_folder = '/Volumes/Lee_Lab/Jack_Files/Genomes/Influenza/WSN/star_indices'
```

## Concatenate Files Together

This should only be done if you have a bunch of reads that you want to process as a single file. When we use the NextSeq, this is typically how it is done (the reads are broken up into multiple files representing different lanes). **Since the RNA for HITS-CLIP is relatively small, we treat them as single-end instead of paired-end, so we combine the R1 and R2 files together into one file (Ask Nara why we do this instead of pairing).**


```python
!mkdir -p $output_dir"/data"

combined_r1_file = output_dir + '/data/' + file_prefixes[0] + '_combined_r1.fastq'
combined_r2_file = output_dir + '/data/' + file_prefixes[0] + '_combined_r2.fastq'
lee_lab_dir = ('/Volumes/Lee_Lab/Deep_Sequencing_Files/Influenza_Project/' +
               'WSN_NP_HITS-CLIP/WSN_NP_HITS-CLIP_090116')

!gzcat $lee_lab_dir/*R1*.gz > $combined_r1_file
!gzcat $lee_lab_dir/*R2*.gz > $combined_r2_file

!gzip $combined_r1_file
!gzip $combined_r2_file

combined_r1_file = combined_r1_file + '.gz'
combined_r2_file = combined_r2_file + '.gz'
```

## Verify user-provided parameters look right


```python
# Main list from which files will be processed.
# Single-End should be in the format:
#  ['f1', 'f2', ..., 'fn']
# Paired-End should be in the format:
#  [('f1r1', 'f1r2'), ('f2r1', f2r2'), ..., ('fnr1', 'fnr2')]
input_files = [(combined_r1_file, combined_r2_file)]

# Testing to make sure valid input
assert len(paired_end) == len(run_fastqc)
assert len(paired_end) == len(has_adapters)
assert len(paired_end) == len(file_prefixes)
assert len(paired_end) == len(input_files)
```

## Static Variables


```python
# Adapters
RL5 = 'AGGGAGGACGATGCGG' # 5' adapter
RL3 = 'GTGTCAGTCACTTCCAGCGG' # 3' adapter
```

## Load Python Libraries


```python
import math
from Bio.Seq import Seq
```

# Pre-process Data

## Get Total Raw Read Counts


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

    [567007]


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

    ('testing_2/data/WSN_combined_r1.fastq.gz', 'testing_2/data/WSN_combined_r2.fastq.gz')
    testing_2/data/WSN_combined_r1.fastq.gz testing_2/data/WSN_combined_r2.fastq.gz
    Started analysis of WSN_combined_r1.fastq.gz
    Approx 5% complete for WSN_combined_r1.fastq.gz
    Approx 10% complete for WSN_combined_r1.fastq.gz
    Approx 15% complete for WSN_combined_r1.fastq.gz
    Approx 20% complete for WSN_combined_r1.fastq.gz
    Approx 25% complete for WSN_combined_r1.fastq.gz
    Approx 30% complete for WSN_combined_r1.fastq.gz
    Approx 35% complete for WSN_combined_r1.fastq.gz
    Approx 40% complete for WSN_combined_r1.fastq.gz
    Approx 45% complete for WSN_combined_r1.fastq.gz
    Approx 50% complete for WSN_combined_r1.fastq.gz
    Approx 55% complete for WSN_combined_r1.fastq.gz
    Approx 60% complete for WSN_combined_r1.fastq.gz
    Approx 65% complete for WSN_combined_r1.fastq.gz
    Approx 70% complete for WSN_combined_r1.fastq.gz
    Approx 75% complete for WSN_combined_r1.fastq.gz
    Approx 80% complete for WSN_combined_r1.fastq.gz
    Approx 85% complete for WSN_combined_r1.fastq.gz
    Approx 90% complete for WSN_combined_r1.fastq.gz
    Approx 95% complete for WSN_combined_r1.fastq.gz
    Approx 100% complete for WSN_combined_r1.fastq.gz
    Analysis complete for WSN_combined_r1.fastq.gz
    Started analysis of WSN_combined_r2.fastq.gz
    Approx 5% complete for WSN_combined_r2.fastq.gz
    Approx 10% complete for WSN_combined_r2.fastq.gz
    Approx 15% complete for WSN_combined_r2.fastq.gz
    Approx 20% complete for WSN_combined_r2.fastq.gz
    Approx 25% complete for WSN_combined_r2.fastq.gz
    Approx 30% complete for WSN_combined_r2.fastq.gz
    Approx 35% complete for WSN_combined_r2.fastq.gz
    Approx 40% complete for WSN_combined_r2.fastq.gz
    Approx 45% complete for WSN_combined_r2.fastq.gz
    Approx 50% complete for WSN_combined_r2.fastq.gz
    Approx 55% complete for WSN_combined_r2.fastq.gz
    Approx 60% complete for WSN_combined_r2.fastq.gz
    Approx 65% complete for WSN_combined_r2.fastq.gz
    Approx 70% complete for WSN_combined_r2.fastq.gz
    Approx 75% complete for WSN_combined_r2.fastq.gz
    Approx 80% complete for WSN_combined_r2.fastq.gz
    Approx 85% complete for WSN_combined_r2.fastq.gz
    Approx 90% complete for WSN_combined_r2.fastq.gz
    Approx 95% complete for WSN_combined_r2.fastq.gz
    Approx 100% complete for WSN_combined_r2.fastq.gz
    Analysis complete for WSN_combined_r2.fastq.gz


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
    
    Forward reads file.................: testing_2/data/WSN_combined_r1.fastq.gz
    Reverse reads file.................: testing_2/data/WSN_combined_r2.fastq.gz
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
      A: 0.239821
      C: 0.225821
      G: 0.323995
      T: 0.210363
      11350 uncalled bases
    Assemblying reads: 100%
    
    Assembled reads ...................: 532,702 / 567,007 (93.950%)
    Discarded reads ...................: 58 / 567,007 (0.010%)
    Not assembled reads ...............: 34,247 / 567,007 (6.040%)
    Assembled reads file...............: testing_2/pear/WSN_paired.assembled.fastq
    Discarded reads file...............: testing_2/pear/WSN_paired.discarded.fastq
    Unassembled forward reads file.....: testing_2/pear/WSN_paired.unassembled.forward.fastq
    Unassembled reverse reads file.....: testing_2/pear/WSN_paired.unassembled.reverse.fastq


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

    [532702]


## Filter by Quality

Uses tool provided by the [CTK](https://github.com/chaolinzhanglab/ctk).


```python
def filterByQuality(fastq_file, output_file):
    # Only consider quality after the length of adapter 5
    adp5_len = len(RL5)
    script = '/usr/local/bin/ctk-master/fastq_filter.pl'
    !perl $script -v -f mean:"$adp5_len"-24:20 -of fastq $fastq_file $output_file
```


```python
filtered_files = []
for i in range(0, len(combined_files)):
    out_file_name = (output_dir + '/data/' + file_prefixes[i] + '_filtered.fastq')
    filterByQuality(combined_files[i], out_file_name)
    !gzip $out_file_name
    out_file_name = out_file_name + '.gz'
    filtered_files.append(out_file_name)
```

    CMD = fastq_filter.pl -v -f mean:16-24:20 -of fastq testing_2/pear/WSN_paired.assembled.fastq testing_2/data/WSN_filtered.fastq
    f=mean:16-24:20
    1 filters detected
    100000 ...
    200000 ...
    300000 ...
    400000 ...
    500000 ...


## Collapse PCR Duplicates

Uses tool provided by the [CTK](https://github.com/chaolinzhanglab/ctk).


```python
def collapsePcrDuplicates(fastq_file, output_file):
    script = '/usr/local/bin/ctk-master/fastq2collapse.pl'
    !perl $script -v $fastq_file $output_file
```


```python
collapsed_files = []
for i in range(0, len(filtered_files)):
    out_file_name = (output_dir + '/data/' + file_prefixes[i] + '_collapsed.fastq')
    collapsePcrDuplicates(filtered_files[i], out_file_name)
    !gzip $out_file_name
    out_file_name = out_file_name + '.gz'
    collapsed_files.append(out_file_name)
```

    gunzip -c testing_2/data/WSN_filtered.fastq.gz | awk '{if(NR%4==1) {print $1}}' > ./fastq2collapse.pl_1539787747_0.546352763268263/id
    gunzip -c testing_2/data/WSN_filtered.fastq.gz | awk '{if(NR%4==2) {print $0}}' > ./fastq2collapse.pl_1539787747_0.546352763268263/seq
    gunzip -c testing_2/data/WSN_filtered.fastq.gz | awk '{if(NR%4==0) {print $0}}' > ./fastq2collapse.pl_1539787747_0.546352763268263/qual
    paste ./fastq2collapse.pl_1539787747_0.546352763268263/id ./fastq2collapse.pl_1539787747_0.546352763268263/qual ./fastq2collapse.pl_1539787747_0.546352763268263/seq| sort -k 3 | uniq -f 2 -c | awk '{print $2"#"$1"\n"$4"\n+\n"$3}' > testing_2/data/WSN_collapsed.fastq


## Get Preprocessed Read Count


```python
total_preprocessed_reads = []
for i in range(0, len(collapsed_files)):
    in_file = collapsed_files[i]
    preprocessed_read_count = 0
    preprocessed_read_count = !gzcat $in_file | wc -l
    preprocessed_read_count = preprocessed_read_count[0]
    total_preprocessed_reads.append(int(int(preprocessed_read_count) / 4))
        
print(total_preprocessed_reads)
```

    [351203]


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
    min_len = 15
    
    print('Running cutadapt on \"' + in_file + '\", outputting to \"' +
          out_trimmed + '\"')
    
    # -a means it will look for anchored 5' adaptor, and if 3' adaptor is present,
    # it will also cut it off.
    #!cutadapt -M $max_len -m $min_len --no-indels -a $adp5"..."$adp3 $e $O --untrimmed-o $out_untrimmed -o $out_trimmed $in_file
    !cutadapt -q 20 -m $min_len --no-indels -a $adp5"..."$adp3 $e $O --untrimmed-o $out_untrimmed -o $out_trimmed $in_file
```


```python
# All files after the cutadapt process will be added to this list,
# even if they don't get trimmed.
trimmed_files = []

!mkdir -p $output_dir"/cutadapt"

revRL3 = ''.join(Seq(RL3).reverse_complement())
revRL5 = ''.join(Seq(RL5).reverse_complement())

for i in range(0, len(collapsed_files)):
    if (has_adapters[i]):

        # Get an output prefix for the cutadapt files
        out_file_base = output_dir + '/cutadapt/' + file_prefixes[i]
        
        # Get the file extension
        file_ext = 'fastq'
        
        # Cut in 5'-->3' Direction
        out_trimmed_1 = out_file_base + '_trimmed_1.fq'
        out_untrimmed_1 = out_file_base + '_untrimmed_1.fq'
        runCutAdapt(RL5, RL3, collapsed_files[i], out_trimmed_1, out_untrimmed_1, err=0.25, min_overlap=10)
        
        # Cut in the 3'-->5' Direction using the reverse complements
        out_trimmed_2 = out_file_base + '_trimmed_2.fq'
        out_untrimmed_2 = out_file_base + '_untrimmed_2.fq'
        runCutAdapt(revRL3, revRL5, out_untrimmed_1, out_trimmed_2, out_untrimmed_2)
        
        # The final file will be call "mySample_trimmed.fastq"
        out_trimmed_final = out_file_base + '_trimmed.fq'
        !cat $out_trimmed_1 $out_trimmed_2 > $out_trimmed_final
        
        # Remove intermediary files
        !rm $out_trimmed_1 $out_trimmed_2 $out_untrimmed_1 $out_untrimmed_2
        
        trimmed_files.append(out_trimmed_final)
        
    else:
        # If no adapters, just add the file to the new list so we can continue
        # with alignment using the new list
        trimmed_files.append(collapsed_files[i])
```

    Running cutadapt on "testing_2/data/WSN_collapsed.fastq.gz", outputting to "testing_2/cutadapt/WSN_trimmed_1.fq"
    This is cutadapt 1.18 with Python 3.7.0
    Command line parameters: -q 20 -m 15 --no-indels -a AGGGAGGACGATGCGG...GTGTCAGTCACTTCCAGCGG -e 0.25 -O 10 --untrimmed-o testing_2/cutadapt/WSN_untrimmed_1.fq -o testing_2/cutadapt/WSN_trimmed_1.fq testing_2/data/WSN_collapsed.fastq.gz
    Processing reads on 1 core in single-end mode ...
    Finished in 5.08 s (14 us/read; 4.15 M reads/minute).
    
    === Summary ===
    
    Total reads processed:                 351,203
    Reads with adapters:                   151,243 (43.1%)
    Reads that were too short:               1,208 (0.3%)
    Reads written (passing filters):       150,035 (42.7%)
    
    Total basepairs processed:    28,037,671 bp
    Quality-trimmed:                     266 bp (0.0%)
    Total written (filtered):      6,560,260 bp (23.4%)
    
    === Adapter 2 ===
    
    Sequence: AGGGAGGACGATGCGG...GTGTCAGTCACTTCCAGCGG; Type: linked; Length: 16+20; 5' trimmed: 151243 times; 3' trimmed: 150760 times
    
    No. of allowed errors:
    0-3 bp: 0; 4-7 bp: 1; 8-11 bp: 2; 12-15 bp: 3; 16 bp: 4
    
    No. of allowed errors:
    0-3 bp: 0; 4-7 bp: 1; 8-11 bp: 2; 12-15 bp: 3; 16-19 bp: 4; 20 bp: 5
    
    Overview of removed sequences at 5' end
    length	count	expect	max.err	error counts
    16	151243	0.0	4	138916 11394 841 80 12
    
    
    Overview of removed sequences at 3' end
    length	count	expect	max.err	error counts
    10	8	0.3	2	8
    11	8	0.1	2	7 0 1
    12	7	0.0	3	7
    13	16	0.0	3	16
    14	10	0.0	3	10
    15	21	0.0	3	17 3 1
    16	15	0.0	4	14 1
    17	99	0.0	4	89 6 3 0 1
    18	996	0.0	4	820 143 23 2 8
    19	1019	0.0	4	789 119 7 12 92
    20	148131	0.0	5	138443 8808 761 73 31 15
    21	61	0.0	5	29 6 5 1 15 5
    23	1	0.0	5	1
    24	1	0.0	5	1
    26	1	0.0	5	1
    27	1	0.0	5	1
    29	4	0.0	5	0 0 0 0 1 3
    30	7	0.0	5	3 1 0 2 1
    32	2	0.0	5	2
    33	1	0.0	5	1
    36	1	0.0	5	0 1
    37	2	0.0	5	1 1
    38	11	0.0	5	2 9
    39	182	0.0	5	181 1
    40	64	0.0	5	63 1
    41	1	0.0	5	1
    42	1	0.0	5	0 0 0 1
    45	2	0.0	5	2
    46	1	0.0	5	0 0 0 0 1
    51	2	0.0	5	2
    53	1	0.0	5	1
    54	1	0.0	5	1
    55	1	0.0	5	1
    57	1	0.0	5	1
    58	6	0.0	5	6
    59	1	0.0	5	1
    60	1	0.0	5	1
    64	1	0.0	5	1
    66	1	0.0	5	1
    67	1	0.0	5	1
    68	4	0.0	5	4
    70	1	0.0	5	1
    74	3	0.0	5	3
    75	1	0.0	5	1
    76	2	0.0	5	2
    77	2	0.0	5	2
    78	3	0.0	5	3
    79	4	0.0	5	4
    80	2	0.0	5	1 1
    81	2	0.0	5	2
    82	4	0.0	5	4
    84	1	0.0	5	1
    86	2	0.0	5	2
    87	1	0.0	5	1
    88	3	0.0	5	3
    89	2	0.0	5	2
    90	1	0.0	5	1
    91	2	0.0	5	2
    93	1	0.0	5	0 1
    95	2	0.0	5	2
    96	2	0.0	5	2
    97	1	0.0	5	1
    99	1	0.0	5	1
    101	1	0.0	5	1
    102	1	0.0	5	1
    106	1	0.0	5	1
    107	2	0.0	5	2
    108	3	0.0	5	3
    109	1	0.0	5	1
    110	1	0.0	5	1
    112	2	0.0	5	2
    113	1	0.0	5	1
    115	1	0.0	5	1
    117	1	0.0	5	1
    118	1	0.0	5	1
    119	1	0.0	5	1
    124	2	0.0	5	2
    126	1	0.0	5	1
    130	1	0.0	5	1
    132	1	0.0	5	1
    
    Running cutadapt on "testing_2/cutadapt/WSN_untrimmed_1.fq", outputting to "testing_2/cutadapt/WSN_trimmed_2.fq"
    This is cutadapt 1.18 with Python 3.7.0
    Command line parameters: -q 20 -m 15 --no-indels -a CCGCTGGAAGTGACTGACAC...CCGCATCGTCCTCCCT --untrimmed-o testing_2/cutadapt/WSN_untrimmed_2.fq -o testing_2/cutadapt/WSN_trimmed_2.fq testing_2/cutadapt/WSN_untrimmed_1.fq
    Processing reads on 1 core in single-end mode ...
    Finished in 3.64 s (18 us/read; 3.29 M reads/minute).
    
    === Summary ===
    
    Total reads processed:                 199,960
    Reads with adapters:                   195,278 (97.7%)
    Reads that were too short:               1,234 (0.6%)
    Reads written (passing filters):       194,044 (97.0%)
    
    Total basepairs processed:    16,018,394 bp
    Quality-trimmed:                       0 bp (0.0%)
    Total written (filtered):      8,591,642 bp (53.6%)
    
    === Adapter 2 ===
    
    Sequence: CCGCTGGAAGTGACTGACAC...CCGCATCGTCCTCCCT; Type: linked; Length: 20+16; 5' trimmed: 195278 times; 3' trimmed: 193827 times
    
    No. of allowed errors:
    0-9 bp: 0; 10-19 bp: 1; 20 bp: 2
    
    No. of allowed errors:
    0-9 bp: 0; 10-16 bp: 1
    
    Overview of removed sequences at 5' end
    length	count	expect	max.err	error counts
    20	195278	0.0	2	182093 12224 961
    
    
    Overview of removed sequences at 3' end
    length	count	expect	max.err	error counts
    3	11	3124.4	0	11
    5	4	195.3	0	4
    6	12	48.8	0	12
    7	3	12.2	0	3
    8	10	3.1	0	10
    9	7	0.8	0	7
    10	38	0.2	1	36 2
    11	51	0.0	1	50 1
    12	111	0.0	1	111
    13	213	0.0	1	206 7
    14	295	0.0	1	228 67
    15	1012	0.0	1	754 258
    16	191040	0.0	1	177854 13186
    17	257	0.0	1	36 221
    18	5	0.0	1	0 5
    27	1	0.0	1	1
    28	4	0.0	1	0 4
    29	2	0.0	1	2
    30	5	0.0	1	2 3
    31	10	0.0	1	2 8
    32	252	0.0	1	251 1
    33	6	0.0	1	6
    34	2	0.0	1	1 1
    35	7	0.0	1	7
    36	9	0.0	1	7 2
    37	2	0.0	1	2
    38	4	0.0	1	1 3
    39	7	0.0	1	4 3
    40	4	0.0	1	3 1
    41	10	0.0	1	10
    42	1	0.0	1	1
    43	2	0.0	1	1 1
    44	5	0.0	1	2 3
    45	2	0.0	1	1 1
    46	15	0.0	1	9 6
    47	5	0.0	1	5
    48	67	0.0	1	67
    49	16	0.0	1	16
    50	7	0.0	1	7
    51	7	0.0	1	7
    52	7	0.0	1	7
    54	13	0.0	1	12 1
    55	3	0.0	1	3
    56	8	0.0	1	7 1
    57	4	0.0	1	4
    58	12	0.0	1	8 4
    59	20	0.0	1	20
    60	16	0.0	1	14 2
    61	7	0.0	1	5 2
    62	5	0.0	1	5
    63	2	0.0	1	2
    64	18	0.0	1	18
    65	4	0.0	1	4
    66	4	0.0	1	4
    67	4	0.0	1	4
    68	8	0.0	1	8
    69	7	0.0	1	6 1
    70	3	0.0	1	3
    71	3	0.0	1	3
    72	4	0.0	1	4
    73	6	0.0	1	6
    74	3	0.0	1	3
    75	7	0.0	1	6 1
    76	13	0.0	1	11 2
    77	6	0.0	1	6
    78	4	0.0	1	4
    79	9	0.0	1	8 1
    80	8	0.0	1	8
    81	3	0.0	1	2 1
    82	2	0.0	1	2
    83	4	0.0	1	4
    84	5	0.0	1	4 1
    85	7	0.0	1	7
    86	3	0.0	1	3
    87	4	0.0	1	4
    88	1	0.0	1	1
    89	4	0.0	1	4
    90	2	0.0	1	2
    91	2	0.0	1	2
    92	2	0.0	1	2
    93	3	0.0	1	3
    94	2	0.0	1	2
    95	3	0.0	1	3
    96	3	0.0	1	3
    97	2	0.0	1	2
    98	4	0.0	1	4
    99	4	0.0	1	4
    100	4	0.0	1	4
    101	5	0.0	1	5
    102	1	0.0	1	1
    103	2	0.0	1	2
    104	3	0.0	1	3
    105	4	0.0	1	4
    106	1	0.0	1	1
    107	3	0.0	1	3
    110	1	0.0	1	1
    111	2	0.0	1	2
    113	2	0.0	1	2
    116	1	0.0	1	1
    119	2	0.0	1	2
    121	2	0.0	1	2
    123	1	0.0	1	1
    127	1	0.0	1	1
    140	1	0.0	1	1
    155	1	0.0	1	0 1
    249	1	0.0	1	1
    


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

    [344079]


## Trim Reads by length??? (80cut)

# Align with STAR

I use STAR because it runs significantly faster than other aligners, e.g. novoalign and hisat2, and provides comparable accuracy to novoalign. Since the question will probably come up for why tophat2 or bowtie2 isn't mentioned, that's because [hisat2 has replaced both of them as superior in every way](https://www.nature.com/articles/nprot.2016.095).

If running on a large genome, run STAR on the htc cluster. STAR takes a lot of memory, as in 32G+ of memory if running with multiple threads and on a big genome. It is also able to massively leverage multithreading. If you need to run a large-genome alignment locally, you might want to consider hisat2 instead.

## Get STAR Genome SA Index Parameter

I don't know how STAR indexing works, but this parameter relies on the size of the genome you're mapping to. For human genomes and anything relatively large, it should be 14. If it starts getting too small, like < 16kbp (roughly 2^14), it will need to change. Though I don't anticipate working with genomes that small, I put the formula here to calculate it just in case (did I have to do this for flu?).


```python
genome_size = 0
for line in open(chrom_file, 'r'):
    chrom, length = line.strip().split('\t')
    genome_size += int(length)
log2_genome_size = int(math.log(genome_size, 2))
genomeSAindexNbases = min(14, log2_genome_size)
```

## Run Alignment


```python
bam_files = []

!mkdir -p $output_dir"/star"

for i in range(0, len(trimmed_files)):
    
    file_to_align = trimmed_files[i]
    prefix = output_dir + '/star/' + file_prefixes[i] + '_'
    
    params = ('--outFileNamePrefix ' + prefix +
              ' --runThreadN ' + str(num_threads) +
              # Folder containing previously-built STAR indices
              ' --genomeDir ' + star_idx_folder +
              # Produce a sorted BAM file
              ' --outSAMtype BAM SortedByCoordinate' +
              # ???
              ' --alignIntronMax 1' +
              # ???
              ' --alignMatesGapMax ' + str(trim_size) +
              ' --readFilesIn ' + file_to_align +
              # See previous cell
              ' --genomeSAindexNbases ' + str(genomeSAindexNbases) +
              # This is to limit the amount of RAM used by STAR 
              ' --limitBAMsortRAM 1838117759' +
              ' --outSAMunmapped Within' +
              # Need to ensure we have the MD tag in particular so
              # finding mutations goes smoothly
              ' --outSAMattributes NH HI NM MD')
              
    !STAR $params

    bam_file = prefix + 'Aligned.sortedByCoord.out.bam'
    bam_files.append(bam_file)
```

    Oct 17 10:50:24 ..... started STAR run
    Oct 17 10:50:24 ..... loading genome
    Oct 17 10:50:24 ..... started mapping
    Oct 17 10:50:30 ..... started sorting BAM
    Oct 17 10:50:30 ..... finished successfully


## Convert to BED


```python
bed_files = []
mut_files = []

!mkdir -p $output_dir"/bed"

for i in range(0, len(bam_files)):
    bam_file = bam_files[i]
    mut_file = output_dir + '/star/' + file_prefixes[i] + '_mutation.txt'
    bed_file = output_dir + '/bed/' + file_prefixes[i] + '.bed'
    script1 = ('/usr/local/bin/ctk-master/parseAlignment.pl ' +
               '-v' +
               ' --mutation-file ' + mut_file +
               ' --map-qual 255' +
               ' --min-len 18' +
               # This allows input from stdin (pipe)
               ' - ' +
               bed_file)
    
    !samtools view $bam_file | perl $script1
    
    bed_files.append(bed_file)
    mut_files.append(mut_file)
```

    0 ...
    50000 ...
    100000 ...
    150000 ...
    200000 ...
    250000 ...
    300000 ...


### Get Aligned Read Counts


```python
total_aligned_reads = []
for i in range(0, len(bed_files)):
    in_file = bed_files[i]
    aligned_read_count = 0
    aligned_read_count = !wc -l $in_file
    aligned_read_count = aligned_read_count[0].split()[0]

    total_aligned_reads.append(int(aligned_read_count))
        
print(total_aligned_reads)
```

    [215503]


# Post-Alignment Analysis

## Make bedgraphs

### Sort Bed Files


```python
for i in range(0, len(bed_files)):
    bed_file = bed_files[i]
    tmp_bed_file = bed_file + '_tmp'
    # Create a sorted bed file, then replace the unsorted bed file with it 
    !sort -k1,1 -k2,2n $bed_file > $tmp_bed_file
    !mv $tmp_bed_file $bed_file
```

### Run genomeCoverageBed


```python
!mkdir -p $output_dir"/bedgraph"

bedgraph_files = []

for i in range(0, len(bed_files)):
    bed_file = bed_files[i]
    bedgraph_base_name = file_prefixes[i] + '_' + str(trim_size) + '_%s.bedgraph'
    bedgraph_pos_base = bedgraph_base_name % 'pos'
    bedgraph_neg_base = bedgraph_base_name % 'neg'
    bedgraph_pos = output_dir + '/bedgraph/' + bedgraph_pos_base
    bedgraph_neg = output_dir + '/bedgraph/' + bedgraph_neg_base
    
    # Include track lines for viewing in UCSC
    track_line = ('track type=bedGraph windowingFunction=maximum ' +
                  'visibility=full name=%s color=0,0,255 altColor=255,0,0')
    pos_track_line = track_line % bedgraph_pos_base
    neg_track_line = track_line % bedgraph_neg_base

    !echo $pos_track_line > $bedgraph_pos
    !echo $neg_track_line > $bedgraph_neg
    !genomeCoverageBed -bg -strand + -i $bed_file -g $chrom_file >> $bedgraph_pos
    !genomeCoverageBed -bg -strand - -i $bed_file -g $chrom_file >> $bedgraph_neg
    
    bedgraph_files.append(bedgraph_pos)
```

## CIMS

### Generate CIMS bedgraphs


```python
def makeCimsBedgraph(cim_folder, cim_file_base, mut_symbol, mut_name, mut_file):
    awk_str = (('\'{if($9==\"%s\"){print $1\"\\t\"$2"' +
               '\\t"$3"\\t"$4"\\t"$5"\\t"$6}}\'') % mut_symbol)
    cim_file_pos = cim_file_base % (mut_name, 'pos')
    cim_file_neg = cim_file_base % (mut_name, 'neg')
    
    track_line = ('track type=bedGraph windowingFunction=maximum ' +
                  'visibility=full name=%s color=0,0,255 altColor=255,0,0')
    pos_track_line = track_line % cim_file_pos
    neg_track_line = track_line % cim_file_neg
    
    cim_file_pos = cim_folder + cim_file_pos
    cim_file_neg = cim_folder + cim_file_neg
    
    !echo $pos_track_line > $cim_file_pos
    !echo $neg_track_line > $cim_file_neg
    
    script = ('awk ' + awk_str + ' ' + mut_file + ' | sort -k1,1 -k2,2n | ' +
              'genomeCoverageBed -bg -i - -g ' + chrom_file)
    
    !$script -strand + >> $cim_file_pos
    !$script -strand - >> $cim_file_neg
```


```python
for i in range(0, len(mut_files)):
    mut_file = mut_files[i]
    cim_folder = output_dir + '/bedgraph/'
    cim_file_base =  file_prefixes[i] + '_%s_%s.bedgraph'
    makeCimsBedgraph(cim_folder, cim_file_base, '-', 'del', mut_file)
    makeCimsBedgraph(cim_folder, cim_file_base, '+', 'ins', mut_file)
    makeCimsBedgraph(cim_folder, cim_file_base, '>', 'sub', mut_file)
```

## Call Peaks


```python
import pysam

def getCoverage(bam_file, genome_size):
    total_num_nts = 0
    bam_file = pysam.AlignmentFile(bam_file, "rb")
    x = 0
    for alignedSegment in bam_file:
        if not (alignedSegment.is_unmapped):
            total_num_nts += alignedSegment.infer_read_length()
    return int(total_num_nts / genome_size)

    '''
    awk_str = ('\'{if($5=="255") {total = total + length($10)}}END' +
               '{printf "%.0f\\n", int(total/genomeSize)==(total' +
               '/genomeSize)?(total/genomeSize):int(total/genomeSize)+1}\'')
    script = ('samtools view ' + bam_file + ' | awk -v genomeSize=' +
              str(genome_size) + ' ' + awk_str)
    print(script)
    coverage = !script
    print(coverage)
    '''
```


```python
print(getCoverage(bam_files[0], genome_size))
```

    711



```python
peak_counts = []
peak_beds = []

for i in range(0, len(bed_files)):
    cov = getCoverage(bam_files[0], genome_size)
    val_depth = 0.5
    script_path = '/usr/local/bin/ctk-master/tag2peak.pl'
    
    bed_file = bed_files[i]
    out_bed_base = (output_dir + '/bed/' + file_prefixes[i] + '_' +
                    str(val_depth) + '_' + str(cov) + '_%s.bed')
    halfph_bed = out_bed_base % 'halfPH'
    peaks_bed = out_bed_base % 'peaks'
    
    params = ('-v' +
              ' -ss' +
              ' --valley-seeking' +
              ' --valley-depth ' + str(val_depth) +
              ' --minPH ' + str(cov) +
              ' --out-half-PH ' + halfph_bed +
              ' ' + bed_file +
              ' ' + peaks_bed)
    
    !perl $script_path $params
    
    peak_beds.append(halfph_bed)
    
    num_peaks = !wc -l $halfph_bed
    num_peaks = num_peaks[0].split()[0]
    peak_counts.append(num_peaks)
```

    CMD=tag2peak.pl -v -ss --valley-seeking --valley-depth 0.5 --minPH 711 --out-half-PH testing_2/bed/WSN_0.5_711_halfPH.bed testing_2/bed/WSN.bed testing_2/bed/WSN_0.5_711_peaks.bed
    generate exact tag coverage profile ...
    CMD=tag2profile.pl -v -exact -ss -of bed testing_2/bed/WSN.bed ./tag2peak.pl_1539787857_0.936259798313589/tag.count.exact.bed
    reading tags from testing_2/bed/WSN.bed ...
    0 ...
    100000 ...
    200000 ...
    get tag count broken down into chromosomes ...
    HA : 32189 tags
    M : 16253 tags
    NA : 22394 tags
    NP : 31554 tags
    NS : 15304 tags
    PA : 33339 tags
    PB1 : 31292 tags
    PB2 : 33178 tags
    #processing strand + ...
    processing chromsome HA ...
    32189 tags loaded on chromosome HA
    processing chromsome M ...
    16253 tags loaded on chromosome M
    processing chromsome NA ...
    22394 tags loaded on chromosome NA
    processing chromsome NP ...
    31554 tags loaded on chromosome NP
    processing chromsome NS ...
    15304 tags loaded on chromosome NS
    processing chromsome PA ...
    33339 tags loaded on chromosome PA
    processing chromsome PB1 ...
    31292 tags loaded on chromosome PB1
    processing chromsome PB2 ...
    33178 tags loaded on chromosome PB2
    #processing strand - ...
    processing chromsome HA ...
    32189 tags loaded on chromosome HA
    processing chromsome M ...
    16253 tags loaded on chromosome M
    processing chromsome NA ...
    22394 tags loaded on chromosome NA
    processing chromsome NP ...
    31554 tags loaded on chromosome NP
    processing chromsome NS ...
    15304 tags loaded on chromosome NS
    processing chromsome PA ...
    33339 tags loaded on chromosome PA
    processing chromsome PB1 ...
    31292 tags loaded on chromosome PB1
    processing chromsome PB2 ...
    33178 tags loaded on chromosome PB2
    extract candidate peaks ...
    0 ...
    get positions with non-zero tag count broken down into chromosomes ...
    HA : 2813 tags
    M : 1657 tags
    NA : 2369 tags
    NP : 2663 tags
    NS : 1508 tags
    PA : 3666 tags
    PB1 : 3748 tags
    PB2 : 3895 tags
    #processing strand + ...
    processing chrom HA ...
    2813 tags loaded on chromosome HA
    processing chrom M ...
    1657 tags loaded on chromosome M
    processing chrom NA ...
    2369 tags loaded on chromosome NA
    processing chrom NP ...
    2663 tags loaded on chromosome NP
    processing chrom NS ...
    1508 tags loaded on chromosome NS
    processing chrom PA ...
    3666 tags loaded on chromosome PA
    processing chrom PB1 ...
    3748 tags loaded on chromosome PB1
    processing chrom PB2 ...
    3895 tags loaded on chromosome PB2
    #processing strand - ...
    processing chrom HA ...
    2813 tags loaded on chromosome HA
    processing chrom M ...
    1657 tags loaded on chromosome M
    processing chrom NA ...
    2369 tags loaded on chromosome NA
    processing chrom NP ...
    2663 tags loaded on chromosome NP
    processing chrom NS ...
    1508 tags loaded on chromosome NS
    processing chrom PA ...
    3666 tags loaded on chromosome PA
    processing chrom PB1 ...
    3748 tags loaded on chromosome PB1
    processing chrom PB2 ...
    3895 tags loaded on chromosome PB2
    71 candidate peaks with PH>=711 detected
    generating output...


## Run IGV


```python
for i in range(0, len(bed_files)):
    bedgraph_file = bedgraph_files[i]
    peak_file = peak_beds[i]
    
    print(genome_file)
    print(bedgraph_file)
    print(peak_file)
    
    !echo $bedgraph_file
    !igv -g $genome_file $bedgraph_file,$peak_file
```

    /Volumes/Lee_Lab/Jack_Files/Genomes/Influenza/WSN/WSN.fasta
    testing_2/bedgraph/WSN_65_pos.bedgraph
    testing_2/bed/WSN_0.5_711_halfPH.bed
    testing_2/bedgraph/WSN_65_pos.bedgraph
    INFO [2018-10-17T11:13:25,903]  [Main.java:473] [main]  Loading: testing_2/bedgraph/WSN_65_pos.bedgraph,testing_2/bed/WSN_0.5_711_halfPH.bed
    WARN [2018-10-17T11:13:26,083]  [Globals.java:138] [AWT-EventQueue-0]  Development mode is enabled
    INFO [2018-10-17T11:13:26,097]  [DirectoryManager.java:179] [AWT-EventQueue-0]  IGV Directory: /Users/jpk90/igv
    INFO [2018-10-17T11:13:26,110]  [Main.java:155] [AWT-EventQueue-0]  Startup  IGV Version 2.4.13 07/11/2018 02:16 PM
    INFO [2018-10-17T11:13:26,110]  [Main.java:156] [AWT-EventQueue-0]  Java 1.8.0_172
    INFO [2018-10-17T11:13:26,110]  [DirectoryManager.java:84] [AWT-EventQueue-0]  Fetching user directory... 
    INFO [2018-10-17T11:13:26,198]  [Main.java:157] [AWT-EventQueue-0]  Default User Directory: /Users/jpk90
    INFO [2018-10-17T11:13:26,199]  [Main.java:158] [AWT-EventQueue-0]  OS: Mac OS X
    
    
    INFO [2018-10-17T11:13:33,839]  [GenomeManager.java:186] [pool-3-thread-1]  Loading genome: /Volumes/Lee_Lab/Jack_Files/Genomes/Influenza/WSN/WSN.fasta
    INFO [2018-10-17T11:13:33,865]  [GenomeComboBox.java:79] [AWT-EventQueue-0]  Enter genome combo box
    INFO [2018-10-17T11:13:34,151]  [GenomeComboBox.java:79] [AWT-EventQueue-0]  Enter genome combo box
    INFO [2018-10-17T11:13:34,157]  [GenomeManager.java:277] [pool-3-thread-1]  Genome loaded.  id= /Volumes/Lee_Lab/Jack_Files/Genomes/Influenza/WSN/WSN.fasta
    INFO [2018-10-17T11:13:34,317]  [IGV.java:1386] [pool-3-thread-2]  Loading 2 resources.
    INFO [2018-10-17T11:13:34,318]  [TrackLoader.java:132] [pool-3-thread-2]  Loading resource, path testing_2/bedgraph/WSN_65_pos.bedgraph
    INFO [2018-10-17T11:13:34,323]  [CommandListener.java:120] [Thread-8]  Listening on port 60151
    INFO [2018-10-17T11:13:34,422]  [TrackLoader.java:132] [pool-3-thread-2]  Loading resource, path testing_2/bed/WSN_0.5_711_halfPH.bed
    INFO [2018-10-17T11:21:45,080]  [ShutdownThread.java:47] [Thread-3]  Shutting down


## General Stats:

| Sample Name  | Raw | Paired | Preprocessed | With Adapters | Aligned | Peaks |
| :----------  | --- | ------ | ------------ | ------------- | ------- | ----- |
| WSN  | {{total_raw_reads[0]}} | {{total_peared_reads[0]}} | {{total_preprocessed_reads[0]}} | {{total_trimmed_reads[0]}} | {{total_aligned_reads[0]}} | {{peak_counts[0]}}
