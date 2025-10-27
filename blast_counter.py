#!/usr/bin/env python3
"""
Extract sequences from BLAST results and count copies of each sequence.
Processes BLAST tabular format output and uses blastdbcmd to retrieve sequences.
"""


import csv
import subprocess
from collections import Counter
from typing import Dict
import sys
import argparse
import os


def extract_sequence_from_blast_db(db_path: str, sacc: str, start: int, end: int, strand: str = "plus") -> str:
    """
    Extract a sequence from a BLAST database using blastdbcmd.
    
    Args:
        db_path: Full path to the BLAST database
        sacc: Subject sequence ID
        start: Start position
        end: End position
        strand: Strand orientation ("plus" or "minus")
    
    Returns:
        Extracted sequence as a string (header removed)
    """
    try:
        # Ensure start <= end for blastdbcmd
        if start > end:
            start, end = end, start
            strand = "plus" if strand == "plus" else "minus"
        
        cmd = [
            "blastdbcmd",
            "-db", db_path,
            "-entry", sacc,
            "-range", f"{start}-{end}",
            "-strand", strand
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse FASTA output and extract sequence (skip header line)
        data = result.stdout
        lines = result.stdout.strip().split('\n')
        sequence = ''.join(lines[1:])  # Skip the first line (FASTA header)
        return sequence, data
        
    except subprocess.CalledProcessError as e:
        print(f"Error extracting sequence for {sacc}: {e.stderr}", file=sys.stderr)
        return ""
    except Exception as e:
        print(f"Unexpected error for {sacc}: {str(e)}", file=sys.stderr)
        return ""


def process_blast_results(csv_file: str, db_path: str, min_pident: float = 0.0, max_evalue: float = float('inf'), min_pct_length: float = 0.0) -> Dict[str, int]:
    """
    Process BLAST results CSV and extract sequences.
    
    Args:
        csv_file: Path to the BLAST results CSV file
        db_path: Full path to the BLAST database
        min_pident: Minimum percent identity threshold (default: 0.0)
        max_evalue: Maximum e-value threshold (default: infinity)
        min_pct_length: Minimum percent of query length that must be aligned (default: 0.0)
    
    Returns:
        Dictionary mapping sequences to their counts
    """
    sequences = []
    filtered_count = 0
    datas = []
    
    with open(csv_file, 'r') as f:
        # Create a regular csv.reader first
        reader = csv.reader(f, delimiter='\t')  # BLAST format 6 uses tabs
        
        # Check if first line looks like a header (contains column names)
        first_line = next(reader)
        
        # If the first line contains text headers (not numeric), skip it
        # Otherwise, process it as the first data row
        if first_line and not first_line[0].replace('.', '').replace('_', '').replace('-', '').isalnum():
            # Looks like a header, skip it
            pass
        else:
            # First line is data, process it
            try:
                pident = float(first_line[11])
                qlength = int(first_line[1])
                # length = int(first_line[4])
                sstart = int(first_line[9])
                send = int(first_line[10])
                evalue = float(first_line[14])
                
                # Calculate query coverage percentage
                length = abs(sstart - send) + 1
                pct_length = (length / qlength) * 100
                
                # Apply filters
                if pident >= min_pident and evalue <= max_evalue and pct_length >= min_pct_length:
                    sacc = first_line[7]
                    sstart = int(first_line[9])
                    send = int(first_line[10])
                    
                    strand = "plus" if sstart <= send else "minus"
                    print(f"Processing hit 1: {sacc} ({sstart}-{send}), pident={pident}, evalue={evalue}, pct_length={pct_length:.1f}%", file=sys.stderr)
                    
                    sequence, data = extract_sequence_from_blast_db(db_path, sacc, sstart, send, strand)
                    if sequence:
                        sequences.append(sequence)
                        datas.append(data)
                else:
                    filtered_count += 1
                    print(f"Filtering hit 1: pident={pident} (min={min_pident}), evalue={evalue} (max={max_evalue}), pct_length={pct_length:.1f}% (min={min_pct_length})", file=sys.stderr)
            except (ValueError, IndexError) as e:
                print(f"Warning: Skipping first line due to parsing error: {e}", file=sys.stderr)
        
        # Process remaining rows
        for i, row in enumerate(reader, 2):
            if not row or len(row) < 11:
                continue
                
            try:
                pident = float(row[11])
                qlength = int(row[1])
                # length = int(row[4])
                sstart = int(row[9])
                send = int(row[10])
                evalue = float(row[14])
                
                # Calculate query coverage percentage
                length = abs(sstart - send) + 1
                pct_length = (length / qlength) * 100
                
                # Apply filters
                if pident >= min_pident and evalue <= max_evalue and pct_length >= min_pct_length:
                    sacc = row[7]
                    sstart = int(row[9])
                    send = int(row[10])
                    
                    # Determine strand based on start/end positions
                    strand = "plus" if sstart <= send else "minus"
                    
                    print(f"Processing hit {i}: {sacc} ({sstart}-{send}), pident={pident}, evalue={evalue}, pct_length={pct_length:.1f}%", file=sys.stderr)
                    
                    sequence, data = extract_sequence_from_blast_db(db_path, sacc, sstart, send, strand)
                    
                    if sequence:
                        sequences.append(sequence)
                        datas.append(data)
                else:
                    filtered_count += 1
            except (ValueError, IndexError) as e:
                print(f"Warning: Skipping row {i} due to error: {e}", file=sys.stderr)
                continue
    
    if filtered_count > 0:
        print(f"\nFiltered out {filtered_count} hits based on thresholds", file=sys.stderr)
    
    # Count unique sequences
    sequence_counts = Counter(sequences)
    
    return sequence_counts, datas


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Extract sequences from BLAST results and count unique sequences.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s -i blast_results.csv -db /path/to/database/name
  %(prog)s --input results.csv --database ~/blastdb/nt -o counts.csv
  %(prog)s -i results.csv -db ~/db/mydb --min-pident 95 --max-evalue 1e-10
  %(prog)s -i results.csv -db ~/db/mydb --min-pct-length 80
        """
    )
    
    # Required arguments
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input CSV file containing BLAST results (tabular format -outfmt 6)'
    )
    
    parser.add_argument(
        '-db', '--database',
        required=True,
        help='Path to BLAST database (full path including database name, e.g., /path/to/db/mydb)'
    )
    
    # Optional arguments
    parser.add_argument(
        '-o', '--output',
        default='sequence_counts.csv',
        help='Output CSV file for sequence counts (default: sequence_counts.csv)'
    )
    
    parser.add_argument(
        '--min-pident', "-mip",
        type=float,
        default=0.0,
        help='Minimum percent identity threshold (default: 0.0, accepts all)'
    )
    
    parser.add_argument(
        '--max-evalue', "-me",
        type=float,
        default=float('inf'),
        help='Maximum e-value threshold (default: infinity, accepts all)'
    )
    
    parser.add_argument(
        '--min-pct-length', "-mipl",
        type=float,
        default=0.0,
        help='Minimum percentage of query length that must be aligned (default: 0.0, accepts all)'
    )
    
    return parser.parse_args()


def main():
    """Main function to run the sequence extraction and counting."""
    
    # Parse command-line arguments
    args = parse_arguments()
    
    csv_file = args.input
    db_path = args.database
    output_file = args.output
    min_pident = args.min_pident
    max_evalue = args.max_evalue
    min_pct_length = args.min_pct_length
    
    # Validate input file exists
    if not os.path.exists(csv_file):
        print(f"Error: Input file '{csv_file}' not found.", file=sys.stderr)
        sys.exit(1)
    
    # Validate database exists (check for .nhr, .nin, .nsq files for nucleotide DB)
    db_files_exist = (
        os.path.exists(f"{db_path}.nhr") or 
        os.path.exists(f"{db_path}.phr") or
        os.path.exists(f"{db_path}.00.nhr") or
        os.path.exists(f"{db_path}.00.phr")
    )
    
    if not db_files_exist:
        print(f"Warning: Could not find BLAST database files for '{db_path}'.", file=sys.stderr)
        print("Make sure the path includes the database name (e.g., /path/to/db/mydb)", file=sys.stderr)
    
    print(f"Processing BLAST results from: {csv_file}", file=sys.stderr)
    print(f"Using BLAST database: {db_path}", file=sys.stderr)
    print(f"Filters: min_pident={min_pident}, max_evalue={max_evalue}, min_pct_length={min_pct_length}%", file=sys.stderr)
    
    # Process results and count sequences
    sequence_counts, data = process_blast_results(csv_file, db_path, min_pident, max_evalue, min_pct_length)
    
    # Output results
    print(f"\nTotal unique sequences: {len(sequence_counts)}", file=sys.stderr)
    print(f"Total sequences extracted: {sum(sequence_counts.values())}", file=sys.stderr)
    
    # Write results to CSV
    with open(f"{output_file}.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['sequence', 'count'])
        
        for sequence, count in sequence_counts.most_common():
            writer.writerow([sequence, count])
    # Write results to fasta
    with open(f"{output_file}.fasta", "w", newline = "\n") as f:
        for data in data:
            f.write(data)
    
    print(f"\nResults written to: {output_file}", file=sys.stderr)
    
    # Print summary of top sequences
    print("\nTop 10 most common sequences:", file=sys.stderr)
    for i, (seq, count) in enumerate(sequence_counts.most_common(10), 1):
        preview = seq[:50] + "..." if len(seq) > 50 else seq
        print(f"{i}. Count: {count}, Sequence: {preview}", file=sys.stderr)


if __name__ == "__main__":
    main()
