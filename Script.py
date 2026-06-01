from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_sequences(input_file: Path, num_sequences: int = 10):
    records = list(SeqIO.parse(str(input_file), "fasta"))
    extracted_sequences = records[:num_sequences]
    remaining_sequences = records[num_sequences:]
    return extracted_sequences, remaining_sequences


def create_bed_from_sequences(sequences):
    lines = []
    for i, seq_record in enumerate(sequences, start=1):
        lines.append(f"{seq_record.id}\t0\t{len(seq_record)}\tSequence_{i}")
    return "\n".join(lines) + ("\n" if lines else "")


def reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement())


def modify_sequences(sequences):
    modified_sequences = []
    for seq_record in sequences:
        original_sequence = str(seq_record.seq)
        reversed_complement_sequence = reverse_complement(original_sequence)
        if reversed_complement_sequence:
            modified_sequence = reversed_complement_sequence[:-1] + "A"
        else:
            modified_sequence = "A"
        modified_sequences.append(
            SeqRecord(
                Seq(modified_sequence),
                id=seq_record.id,
                description=seq_record.description,
            )
        )
    return modified_sequences


def process_fasta(
    input_file: Path,
    output_prefix: str,
    num_sequences: int,
    output_dir: Path,
):
    extracted_sequences, remaining_sequences = extract_sequences(input_file, num_sequences=num_sequences)

    output_dir.mkdir(parents=True, exist_ok=True)
    extracted_path = output_dir / f"{output_prefix}_extracted.fasta"
    SeqIO.write(extracted_sequences, extracted_path, "fasta")

    bed_content = create_bed_from_sequences(extracted_sequences)
    bed_path = output_dir / f"{output_prefix}.bed"
    bed_path.write_text(bed_content)

    modified_sequences = modify_sequences(extracted_sequences)
    all_sequences = remaining_sequences + modified_sequences
    modified_output_file_path = output_dir / f"{output_prefix}_modified.fasta"
    SeqIO.write(all_sequences, modified_output_file_path, "fasta")

    return {
        "extracted_path": extracted_path,
        "bed_path": bed_path,
        "modified_path": modified_output_file_path,
        "bed_content": bed_content,
        "extracted_sequences": extracted_sequences,
        "modified_sequences": modified_sequences,
        "all_sequences": all_sequences,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract sequences from a FASTA file, generate BED, and write modified FASTA output.",
    )
    default_input = Path(__file__).with_name("influenza.fna")
    parser.add_argument(
        "--input",
        type=Path,
        default=default_input,
        help="Path to the input FASTA file.",
    )
    parser.add_argument(
        "--num-sequences",
        type=int,
        default=10,
        help="Number of sequences to extract from the input file.",
    )
    parser.add_argument(
        "--output-prefix",
        default="influenza",
        help="Prefix for output files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path.cwd(),
        help="Directory to write output files.",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    if not args.input.exists():
        parser.error(f"Input file not found: {args.input}")

    process_fasta(
        input_file=args.input,
        output_prefix=args.output_prefix,
        num_sequences=args.num_sequences,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
