from tabulate import tabulate
import argparse
from bio_seq_v1.fasta import fasta_parser
from bio_seq_v1.stats import sequence

def print_sequence_lengths_formatted(sequences):
    """
    Print a formatted table of sequence IDs and their lengths.

    Args:
        sequences (list of sequence): List of Sequence objects to process.
    """
    table = [[s.id, s.sequence_length()] for s in sequences]
    print(tabulate(table, headers=["Sequence ID", "Length"], tablefmt="grid"))

def print_gc_content_table(sequences):
    """
    Print a formatted table of sequence IDs and their GC content percentages.

    Args:
        sequences (list of sequence): List of Sequence objects to process.
    """
    table = [[s.id, f"{s.gc_content():.2f}%"] for s in sequences]
    print(tabulate(table, headers=["Sequence", "GC%"], tablefmt="grid"))

def print_revcomp(sequences):
    """
    Print the reverse complement of each sequence in the list.

    Args:
        sequences (list of sequence): List of Sequence objects to process.
    """
    for s in sequences:
        print(f">{s.id} reverse complement")
        print(s.rev_complement())
        print("-" * 30)

def print_base_count(sequences):
    """
    Print a table of the counts of each base for each sequence.

    Args:
        sequences (list of sequence): List of Sequence objects to process.
    """
    all_counts = [s.base_count() for s in sequences]
    bases_present = [b for b in sequence.valid]
    table = []
    for counts, s in zip(all_counts, sequences):
        row = [s.id] + [counts.get(b, 0) for b in bases_present]
        table.append(row)
    headers = ["Sequence"] + bases_present
    print(tabulate(table, headers=headers, tablefmt="grid"))

def print_summary(sequences):
    """
    Print a full summary of sequences including lengths, GC content, and base composition.

    Args:
        sequences (list of sequence): List of Sequence objects to process.
    """
    print("SEQUENCE LENGTHS")
    print_sequence_lengths_formatted(sequences)
    print()

    print("GC CONTENT")
    print_gc_content_table(sequences)
    print()

    print("BASE COMPOSITION")
    print_base_count(sequences)

def main():
    """
    Command-line interface entry point.

    Parses arguments to specify FASTA input file and analysis options,
    and prints the requested sequence information.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", "-f", help="Path to the FASTA file", required=True)
    parser.add_argument("--length", "-l", help ="Compute sequence length per sequence", action="store_true")
    parser.add_argument("--gc", help ="Compute GC content per sequence", action="store_true")
    parser.add_argument("--revcomp", "-rc", help ="Compute reverse complements per sequence", action="store_true")
    parser.add_argument("--basecount", "-b", help ="Compute total count for bases per sequence", action="store_true")
    parser.add_argument("--summary", help="Print summary statistics", action="store_true")
    args = parser.parse_args()
    sequences = fasta_parser(args.file)
    if not any([args.length, args.gc, args.revcomp, args.basecount, args.summary]):
        print_summary(sequences)
        exit()
    if not sequences:
        raise ValueError("No sequences found in FASTA file")
    if args.length:
        print_sequence_lengths_formatted(sequences)
        print()
    if args.gc:
        print_gc_content_table(sequences)
        print()
    if args.revcomp:
        print_revcomp(sequences)
        print()
    if args.basecount:
        print_base_count(sequences)
        print()
    if args.summary:
        print_summary(sequences)
        exit()