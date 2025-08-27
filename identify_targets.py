# Dependencies: python 3.12, Biopython, pandas, blast

from typing import Iterator
from pathlib import Path
import argparse
import subprocess
import re

from Bio import Blast, SeqIO
import Bio.Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd


def parse_args() -> object:
    """
    Parses arguments, sets defaults
    """

    parser = argparse.ArgumentParser(
                        prog='blastn IS element flanking regions',
                        description='This program blasts the flanking regions of an insertion sequence, using a high gapopen and gapextend penalty',
                        epilog='Thanks for using! -Chris')
    parser.add_argument("-q", "--query", required = True, help="Path of the query sequence fasta file including extensions")
    parser.add_argument("-db", "--database", required = True, help = "Path of the database name excluding extensions")
    parser.add_argument("-o", "--output", required = False, help = "Output file name. Default is query.out")
    parser.add_argument("-f", "--flank_length", required = False, type = int, default = "200", help = "Length of how long each flanking region is. Default = 200")
    parser.add_argument("-t", "--threads", required = False, default = "1", type = int, help = "Number of threads to use. Default = 1")
    parser.add_argument("-m", "--minimal", action = "store_true", help = "Removes any unnecessary steps for target identification. Disabled by default.")
    parser.add_argument("-r", "--reward", required = False, default = "2", type = int, help = "Reward value. Default = 2. Only change if you need to alter blastn parameters")
    parser.add_argument("-p", "--penalty", required = False, default = "-3", type = int, help = "Penalty value. Default = -3. Only change if you need to alter blastn parameters")
    parser.add_argument("-go", "--gapopen", required = False, default = "0", type = int, help = "Gap open penalty. Default = 0. Only change if you need to alter blastn parameters")
    parser.add_argument("-ge", "--gapextend", required = False, default = "4", type = int, help = "Gap extend penalty. Default =4 . Only change if you need to alter blastn parameters")
    parser.add_argument("-a", "--alignments", required = False, default = "1", type = int, help = "Maximum number of alignments per target sequence. Default = 1. Only change if you need to alter blastn parameters")
    parser.add_argument("-e", "--evalue", required = False, default = "0", type = float, help = "Threshold evalue for filtering insertion sequence blast searches. Default = 0")
    parser.add_argument("-i", "--identity", required = False, default = "0.95", type = float, help = "Threshold percent identity for filtering insertion sequence blast searches. Default = 0.95")
    parser.add_argument("-bi", "--bitscore_insertion_sequence", required = False, type = float, default = "1.5", help = "Threshold bitscore multiplier for filtering insertion sequence blast searches. Lower is more lenient. Default = 1.5")
    parser.add_argument("-bt", "--bitscore_target", required = False, default = "1.1", type = float, help = "Threshold bitscore multiplier for filtering target sequence blasts. Lower is more lenient but would not recommend changing. Default = 1.1")

    parser.add_argument("--other_insertion_sequence", default = "", required = False, type = str, help = "Type any additional parameters for the insertion sequence blast in full")
    parser.add_argument("--other_target", default = "", required = False, type = str, help = "Type any additional parameters for the target blast in full")

    args = parser.parse_args()

    return args


def process_args_output(args: object) -> object:
    """
    Processes the args.output object to direct the exported files to results/output_name \n
    If a directory is provided as args.output, output to that directory instead. \n
    Returns the args object with arguments other than args.output unchanged.
    """

    output_name: str = args.output
    path_string: str = str(args.query)
    path_object: object = Path(path_string)
    script_directory: str = str(Path(__file__).parent.resolve())

    # If no output was provided, get args.query and set the output name to what is provided by args.query.
    # If a directory was provided as the output, add args.query to the end of the output path for final filenames.
    if args.output == "":
        out_path: str = str(path_object.stem)
    elif Path(args.output).is_dir():
        output_name = str(output_name) + "/" + str(path_object.stem)
    else:
        out_path = output_name

    # If the output path doesn't contain a directory, create a results folder in the same directory as the script and save all files in the results folder.
    if str(Path(output_name).parent) == ".":
        results_path: object = Path(script_directory + "/results")
        results_path.mkdir(parents = True, exist_ok = True)
        
        output_name = script_directory + f"/results/{out_path}"
    # Ensure the output directory exists and is avaiable to write in. Create if it does not already exist.
    Path(output_name).mkdir(parents = True, exist_ok = True)
    args.output = str(Path(output_name + f"/{str(path_object.stem)}"))

    return args


def expand_multifasta(args: object) -> Iterator[str]:
    """
    Takes an input fasta file. If there is more than one sequence present in the fasta file, write the individual fastas onto the disk. \n
    Then yield the path to each of the fasta files on disk as a string.
    """

    scriptpath: object = Path(__file__).parent.resolve()
    temp_fastas_path: str = scriptpath / ".temp_fastas"
    temp_fastas_path.mkdir(parents=True, exist_ok=True)

    for record in SeqIO.parse(args.query, "fasta"):
        fasta_name: str = re.sub(r"[/\\?%*:|\"<>\x7F\x00-\x1F ]", "_", record.id)
        output_path: str = temp_fastas_path / f"{fasta_name}.fasta"
        SeqIO.write(record, output_path, "fasta")

        yield str(output_path)


def blast_insertion_sequence(args: object) -> str:
    """
    Takes an input of a fasta file containing an insertion sequence and uses blastn to blast it against a database. \n
    Provides a csv file for easier glancing of the data, and a xml for downstream processing by this program. \n
    
    """

    # I tried to simplify this double search into one by outputting into ASN then parsing through blast_formatter,
    # however, I found that the results seemed to be inconsistent to running the xml separately.
    # The csv blast run is excluded if --minimal is parsed.

    blastn_cmd_csv: str = f'blastn -db {args.database} \
    -query {args.query} \
    -out {args.output}_insertion_seq_blast.csv \
    -outfmt "6 qseqid qlen length qstart qend sseqid stitle sacc slen sstart send pident mismatch gapopen evalue bitscore" \
    -num_threads {args.threads} \
    {args.other_insertion_sequence}'

    blastn_cmd_xml: str = f'blastn -db {args.database} \
    -query {args.query} \
    -out {args.output}_insertion_seq_blast.xml \
    -outfmt 5 \
    -num_threads {args.threads} \
    {args.other_insertion_sequence}'

    if not args.minimal:
        subprocess.run(blastn_cmd_csv, shell = True)
    subprocess.run(blastn_cmd_xml, shell = True)
    print("Insertion sequence blasted")

    # Returns the location of the blast xml output file.
    blast_output_path: str = f"{args.output}_insertion_seq_blast.xml"

    return blast_output_path


def read_blastn(blastn_output: str) -> Blast.Record:
    """
    Reads the blastn output xml file and returns the output data as an Bio.Blast.Record object
    """
    result_stream = open(f"{blastn_output}", "rb")
    blastn_results: Blast.Record = Blast.read(result_stream)
    print("Blast results read")

    return blastn_results
    

def process_and_write_hits(blastn_record: object, args: object) -> list[object]:
    """
    Takes a Bio.Blast.Record object and makes a dictionary containing: \n
    The target id as the key and an array containing the aligned sequence and the coordinates of the aligned sequence as the value. \n
    It then writes this dictionary as a json and fasta format and returns it.
    """
    # Create an array of the alignments sequentially for parsing into the extract_alignment_target_seq function
    i = 0
    alignments: list[object] = []
    for hit in blastn_record:
        j = 0
        for alignment in blastn_record[i]:
            alignments.append(blastn_record[i][j])
            j += 1
        i += 1
    
    filtered_alignments: list[object] = filter_alignments(alignments, args)
    print(f"Filtering {len(alignments)} alignments: resulting in {len(filtered_alignments)} satisfactory insertion sequence alignments.")

    return filtered_alignments

    
def filter_alignments(alignments: object, args: object) -> list[object]:
    """
    Takes a list of Biopython alignment objects and filters reads according to their evalue, bit score, and percent identity to the query. \n
    Returns a list object of Biopython alignment objects.
    """
    i = 0
    satisfactory_alignments: list[object] = []
    for alignment in alignments:
        if ((alignments[i].annotations["evalue"] <= args.evalue) 
        and (alignments[i].annotations["bit score"] > args.bitscore_insertion_sequence * len(alignments[i].query.seq)) 
        and (alignments[i].annotations["identity"]/len(alignments[i].query.seq)) >= args.identity):
            satisfactory_alignments.append(alignments[i])
        i += 1
    
    return satisfactory_alignments


def extract_filtered_flanks(filtered_blast_file: str, args: object) -> str:
    """
    Takes a list of Biopython alignment objects and stores the coordinates of the upstream and downstream flanking regions. \n
    Writes to file a list of flank coordinates for input into blastdbcmd. \n
    Then runs the blastdbcmd program to extract the sequence data from each range of coordinates. Returns the location of the blastdbcmd file containing flanks sequence data.
    """

    filtered_alignments: list = []
    records: Iterator[Blast.Record] = Blast.parse(filtered_blast_file)
    for record in records:
        try:
            alignment: object = record
        except IndexError:
            continue

        if alignment.annotations["bit score"] > args.bitscore_target * len(alignment.query.seq):
            filtered_alignments.append(alignment)

    # Extract coordinates of flanking regions
    seen: set = set()
    extracted_results: list = []
    for alignment in filtered_alignments:
        for result in extract_alignment_target_seq([alignment], "name", True):
            target_name, target_seq, low, high, _ = result
            coords: str = f"{low}-{high}"
            uniq_id: str = f"{target_name} {coords}"
            if uniq_id not in seen:
                extracted_results.append((target_name, coords))
                seen.add(uniq_id)

    # DataFrame of target ranges
    alignment_df: object = pd.DataFrame([coord for _, coord in extracted_results], columns=["target_range"])
    flanking_df: object = get_flanking_coords(alignment_df, args.flank_length)

    # Write coords to file for blastdbcmd
    output_path: object = Path(f"{args.output}_filtered_target_coords.fasta")
    with output_path.open("w") as f_out:
        for (target_name, _), (lower, upper) in zip(extracted_results, flanking_df.values.tolist()):
            f_out.write(f"{target_name} {lower}\n")
            f_out.write(f"{target_name} {upper}\n")

    extracted_output_filepath: str = f"{args.output}_flanking_regions.fasta"
    cmd: str = f'blastdbcmd -db {args.database} -entry_batch {output_path} -outfmt "%f" -out {extracted_output_filepath}'
    subprocess.run(cmd, shell=True) 
    print("Alignment flanking regions extracted")

    return extracted_output_filepath


def get_flanking_coords(alignment_df: object, flank_length: int) -> object:
    """
    Takes a pandas dataframe containing the coordinates of an insertion sequence and grabs the coordinates of the upstream and downstream flanks by a specified number of base pairs. \n
    Returns a pandas datafame containing the lower flank and upper flank values.
    """

    target_ranges: object = alignment_df["target_range"].tolist()
    flanking_coords: list[list] = []

    for range in target_ranges:
        ranges: str = range.split("-")
        i = 0
        for coordinate in ranges:
            if i == 0:
                lower_range: str = f"{str(int(coordinate) - int(flank_length) - 1)}-{str(int(coordinate) - 1)}"
                i = 1
            elif i == 1:
                upper_range: str = f"{str(int(coordinate) + 1)}-{str(int(coordinate) + int(flank_length) + 1)}"
        flanking_coords.append([lower_range, upper_range])
    flanking_coords_df = pd.DataFrame(flanking_coords, columns = ["lower_flank", "upper_flank"])

    return(flanking_coords_df)

        
def combine_flanks(fasta_dir: str, args: object) -> str:
    """
    Takes an input fasta file containing separated fasta flank locations, then concatenates them together to create a single sequence.
    Writes combined sequences to a fasta file.
    """

    flanks_length_limit: int = int(args.flank_length) + 50
    flanking_regions: list[SeqRecord] = []
    not_in_seq: int = 0
    # Open the file containing separate flanking regions, append the downstream flank to the end of the upstream one and save in a list.
    with open(f"{fasta_dir}") as handle:
        records: object = SeqIO.parse(handle, "fasta")
        i = 0
        for record in records:
            if i % 2 == 0:
                flanks_upstream: str = record.seq
                flanks_upstream_coords: str = record.id.split(":")[len(record.id.split(":"))-1]
            elif i % 2 != 0:
                flanks_downstream: str = record.seq
                flanks_flanks_downstream_coords: str = record.id.split(":")[len(record.id.split(":"))-1]
                if len(flanks_upstream) <= flanks_length_limit and len(flanks_downstream) <= flanks_length_limit:
                    flanking_regions.append(SeqRecord(flanks_upstream + flanks_downstream, 
                    id = record.id.split(":")[0] + ":" + flanks_upstream_coords + ":" + flanks_flanks_downstream_coords, 
                    description = record.description.split(" ", 1)[len(record.description.split(" ", 1)) - 1]))
                else:
                    not_in_seq += 1
            i += 1
    # Write the combined flanks into a file
    output_path: str = f"{args.output}_concat_flanks.fasta"
    SeqIO.write(flanking_regions, output_path, "fasta")
    if not_in_seq > 0:
        print(f"Flanks combined: excluded {not_in_seq} sequences due to flanks not present in genome file")
    else:
        print("Flanks combined")
    
    print(f"Blasting {len(flanking_regions)} flanking regions")


    return output_path


def blastn_flanking_regions(args: object, flanking_regions_path: str) -> str:
    """
    Parses and stores argparse variables and sends them to the blastn program. \n 
    It then returns the location of the xml output as a string for further processing of the output alignment files by other functions.
    """

    """
    The reward/penalty and gapopen/gapextend values are optional, but by default are set to a balanced reward/penalty of 2/-3 
    and a gapextend and gapopen penalty of 0/4
    """
    flanking_blast_output: str = f"{args.output}_blasted_flanks"

    # Run the assembled blast command using the provided arguments

    blastn_cmd_txt: str = f"blastn -query {flanking_regions_path} -db {args.database} -out {flanking_blast_output}.out -max_hsps {args.alignments} -num_threads {args.threads} " \
        f"-reward {args.reward} -penalty {args.penalty} -gapopen {args.gapopen} -gapextend {args.gapextend}" \
        f" {args.other_target}"
    
    blastn_cmd_xml: str = f"blastn -query {flanking_regions_path} -db {args.database} -out {flanking_blast_output}.xml -max_hsps {args.alignments} -num_threads {args.threads} -outfmt 5 " \
        f"-reward {args.reward} -penalty {args.penalty} -gapopen {args.gapopen} -gapextend {args.gapextend}" \
        f" {args.other_target}"

    """
    # It was faster to run the blast commands separately on my machine (10 core 16GB M4 Macbook Air) than to blast_formatter the ASN file.

    blastn_cmd_ASN = f"blastn -query {flanking_regions_path} -db {args.database} -out {flanking_blast_output}.asn -max_hsps {args.alignments} -num_threads {args.threads} -outfmt 11 " \
        f"-reward {args.reward} -penalty {args.penalty} -gapopen {args.gapopen} -gapextend {args.gapextend}" \
        f" {args.other_target}"
    
    blast_formatter_cmd_xml = f"blast_formatter -archive {flanking_blast_output}.asn -outfmt '5' -out {flanking_blast_output}.xml"

    blast_formatter_cmd_txt = f"blast_formatter -archive {flanking_blast_output}.asn -outfmt '0' -out {flanking_blast_output}.txt"

    subprocess.run(blastn_cmd_ASN, shell=True)
    subprocess.run(blast_formatter_cmd_xml, shell=True)
    subprocess.run(blast_formatter_cmd_txt, shell=True)
    """

    if not args.minimal:
        subprocess.run(blastn_cmd_txt, shell=True)
    subprocess.run(blastn_cmd_xml, shell=True)

    print("Flanks blasted")

    return f"{str(flanking_blast_output)}.xml"


def process_target_hits(flanking_blast_output_dir: str, args: object) -> None:
    """
    Filters blast hits of the concatenated flanking regions by their bit score, and writes satisfactory targets to a file.
    """

    records: Iterator[Blast.Record] = Blast.parse(flanking_blast_output_dir)

    seen_targets: set = set()
    skipped_targets: int = 0
    total_targets: int = 0

    unique_output_path: object = Path(f"{args.output}_unique_targets.fasta")
    all_output_path: object = Path(f"{args.output}_all_targets.fasta")

    with unique_output_path.open("w") as uniq_out, all_output_path.open("w") as all_out:
        record_indexerror = 0
        for i, record in enumerate(records):
            total_targets += 1
            try:
                alignment: object = record[0][0]
            except IndexError:
                record_indexerror += 1
                continue
            # Filter the alignment with the specified bit score multiplier
            if alignment.annotations["bit score"] <= args.bitscore_target * len(alignment.query.seq):
                skipped_targets += 1
                continue

            for result in extract_alignment_target_seq([alignment], "name", True):
                target_name, target_seq, low, high, _ = result
                coords: str = f"{low}-{high}"
                unique_id: str = f"{target_name}:{coords}"
                all_id: str = f"{unique_id}[{i}]"

                all_out.write(f">{all_id}\n{target_seq}\n")
                if unique_id not in seen_targets:
                    uniq_out.write(f">{unique_id}\n{target_seq}\n")
                    seen_targets.add(unique_id)
    if record_indexerror > 0:
        print(f"Record indexerror events:  {record_indexerror}")
    print(f"Excluded {skipped_targets} of {total_targets} ({(skipped_targets / total_targets * 100):.2f}%) target sequences due to low bit score.")
    print(f"Total of {len(seen_targets)} unique targets written to: {unique_output_path}")
    print(f"All {total_targets - skipped_targets} targets written to: {all_output_path}")


def extract_alignment_target_seq(alignments: list, id_type: str, process_reverse_comp: bool) -> Iterator[list]:
    """
    Takes a list of Biopython alignment objects and extracts the target sequence of the alignment. Yields a list of the target id, target sequence, target coordinates, and a reverse comp indicator.
    """

    """
    Thank you Michael Chai @https://github.com/trashmode for optimising this function for me and making it over 95% faster.
    This program would not be viable without your contributions.
    """

    for alignment in alignments:
        reverse_comp: int = 0
        try:
            for ranges in alignment.target.seq._data.defined_ranges:
                target_sequence: str = str(alignment.target.seq[ranges[0]:ranges[1]])
        except AttributeError:
            target_sequence = extract_alignment_sequence_alternate(alignment)

        target_id: str = str(alignment.target.name) if id_type == "name" else str(alignment.target.id).split('|')[1]
        target_lower: int = alignment.coordinates[0][0]
        target_upper: int = alignment.coordinates[0][-1]

        if target_lower > target_upper:
            reverse_comp = 1

        if reverse_comp and process_reverse_comp:
            target_lower, target_upper = target_upper, target_lower
            target_sequence = str(Seq(target_sequence).reverse_complement())
            reverse_comp = 0

        yield [target_id, target_sequence, target_lower, target_upper, reverse_comp]


def extract_alignment_sequence_alternate(alignment: object) -> str:
    """
    An alternative and original method of extracting the target sequence from a Biopython alignment file. \n
    Handles some edge cases where the the alignment.target.seq._data doesn't have the defined_ranges property and the normal method errors out. \n
    Super inefficient because it checks if there is a base present for each coordinate of the target sequence range. 
    """

    coords_sorted = []
    alignment_target_sequence = ""
    alignments_coordinates_length = len(alignment.coordinates[0])
    coords_sorted = sorted(alignment.coordinates[0])

    """
        From the lower coordinate, go basewise to the upper coordinate. If a base exists at the coordinate, append it to the sequence string. If the base does not
        exist, a Bio.Seq.UndefinedSequenceError will be parsed, in which case nothing is appended to the sequence string.
    """

    alignments_coordinates_lower = coords_sorted[0]
    alignments_coordinates_upper = coords_sorted[alignments_coordinates_length - 1]
    j = alignments_coordinates_lower
    while j < alignments_coordinates_upper:
        try:
            alignment_target_sequence += alignment.target.seq[j]
        except Bio.Seq.UndefinedSequenceError:
            pass
        j += 1

    return alignment_target_sequence


def write_file(content: str, filename: str, type: str) -> None:
    """
    Writes or appends any content parsed into this function into a plaintext file. \n
    Has two writing types: \n
    parse the string "w" into type to write(overwrites existing data). \n
    parse the string "a" or "append" into type to append (adds to any existing text in the file).
    """
    mode: str = "w"
    if type == "append" or type == "a":
        mode = "a"
    file = open(filename, f"{mode}")
    file.write(content)
    file.close



def main() -> None:
    args: object = parse_args()
    target_dir = args.output
    for fasta_path in expand_multifasta(args):
        args.output = target_dir
        args.query = fasta_path
        args = process_args_output(args)
        insertion_seq_blastn_output: str = blast_insertion_sequence(args)
        blastn_results: Blast.Record = read_blastn(insertion_seq_blastn_output)
        filtered_alignments: list[object] = process_and_write_hits(blastn_results, args)
        flanks_output_location: str = extract_filtered_flanks(filtered_alignments, args)
        flanking_regions_path: str = combine_flanks(flanks_output_location, args)
        blasted_flanks_path: str = blastn_flanking_regions(args, flanking_regions_path)
        process_target_hits(blasted_flanks_path, args)


if __name__ == "__main__":
    main()
