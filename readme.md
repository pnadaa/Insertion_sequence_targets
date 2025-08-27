## Insertion Sequences Target Prober:

A blast-based pipeline to determine insertion sequence targets.

This is a python program to blast insertion sequence flanking regions to identify their targets. Takes a fasta file input containing insertion sequences of interest and a local blast database for searching.

## Dependencies

Insertion Sequences Target Prober is tested to work under Python 3.12.

The required dependencies are: blast, biopython, pandas

A suitable locally stored blast database is also required to search against. Both the non-redundant nucleotide and locally built databases using makeblastdb have been successfully used.

## Usage

Below is an example of how to use the Insertion Sequence Targets program.

``` bash
python identify_targets.py \
    --query examples/IS15.fasta \
    --database /path/to/blastdb \
    --output examples/IS15 \
    --threads 1
```

## Insertion sequence target identification

To identify insertion sequence targets, run the script using `python insertion_sequence_targets_prober.py` with the following required options:

-   `-q` `--query` str, The path to a fasta file containing insertion sequences
-   `-db` `--database` str, The path to the blast database

### Optional options:

-   `-o` `--output` str, Output name or directory. Default: results/query
-   `-t` `--threads` int, Number of threads for the blast application to use. Default: 1
-   `-f` `--flank_length` int, Length of each flanking region. Default: 200
-   `-m` `--minimal` bool, Removes any unnecessary steps for target identification to produce only target sequences. Reduces compute time by about 40% but excludes target alignments in a human-readable format. Disabled by default\

#### Blast search parameters and filtering options:

-   `-r` `--reward` int, Reward value. Default: 2
-   `-p` `--penalty` int, Penalty value. Default: -3
-   `-go` `--gapopen` int, Gap open cost. Default: 0
-   `-ge` `--gapextend` int, Gap extension cost. Default: 4
-   `-a` `--alignments` int, Max number of alignments per target sequence (max hsps). Not recommended to change. Default: 1
-   `-e` `--evalue` float, Threshold evalue for filtering insertion sequence blast searches. Default: 0
-   `-i` `--identify` float, Threshold percent identity for filtering insertion sequence blast searches. Default: 0.95
-   `-bi` `--bitscore_insertion_sequence` float, Threshold bitscore multiplier for filtering insertion sequence blast searches. Lower is more lenient. Default: 1.5
-   `-bt` `--bitscore_target` float, Threshold bitscore multiplier for filtering target sequence blasts. Lower is more lenient but would not recommend changing. Default: 1.1
-   `--other_insertion_sequence` str, Type any additional parameters for the insertion sequence blast in full.
-   `--other_target` str, Type any additional parameters for the target blast in full

## Outputs

-   `all_targets.fasta` A fasta of all identified targets including replicates
-   `unique_targets.fasta` A fasta of unique identified targets - replicates are removed.
-   `blasted_flanks.out` Blast output file containing the alignments between insertion sequence flanks and their identified targets. Good for finding potential insertions or deletions which occur during the transposition event.