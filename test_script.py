import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import SeqIO

from Script import process_fasta


class TestProcessFasta(unittest.TestCase):
    def test_process_fasta_outputs(self):
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            input_path = tmp_path / "input.fasta"
            input_path.write_text(">seq1\nATGC\n>seq2\nTTAA\n")

            output_dir = tmp_path / "out"
            result = process_fasta(input_path, "test", 1, output_dir)

            self.assertTrue(result["extracted_path"].exists())
            self.assertTrue(result["bed_path"].exists())
            self.assertTrue(result["modified_path"].exists())

            bed_content = result["bed_path"].read_text()
            self.assertEqual(bed_content, "seq1\t0\t4\tSequence_1\n")

            modified_records = list(SeqIO.parse(result["modified_path"], "fasta"))
            ids = {record.id for record in modified_records}
            self.assertEqual(ids, {"seq1", "seq2"})

            seq1_record = next(record for record in modified_records if record.id == "seq1")
            self.assertEqual(str(seq1_record.seq), "GCAA")


if __name__ == "__main__":
    unittest.main()
