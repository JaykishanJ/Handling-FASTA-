🧬 Handling FASTA — Generate BED from Selected Sequences

This repository contains a Python script that performs a series of operations on a **FASTA file** containing *Influenza DNA sequences*.
It leverages the **Biopython** library for sequence extraction, manipulation, and conversion between formats.

⚙️ Overview

The script performs the following tasks:

1. **Extracts specific sequences** from a FASTA file
2. **Generates a BED file** from the extracted sequences
3. **Performs sequence modifications** (reverse complement and nucleotide changes)
4. **Combines modified and remaining sequences** into a final FASTA file


## 🧩 Functions

### 1. `extract_sequences(input_file, num_sequences=10)`

* Extracts a user-defined number of sequences (default: 10) from the input FASTA file.
* Returns both the **extracted** and **remaining** sequences.

### 2. `create_bed_from_sequences(sequences)`

* Converts a list of FASTA sequences into **BED format**.
* BED format is commonly used to represent genomic features and their chromosomal coordinates.

### 3. `reverse_complement(sequence)`

* Computes the **reverse complement** of a given DNA sequence.
* Useful for orientation correction or strand-based analysis.

---

## 🧠 Workflow

1. **Input and Output Setup**

   * Input file: `Influenza.fna`
   * Output files:

     * `influenza_extracted.fasta`
     * `influenza.bed`
     * `influenza_modified.fasta`

2. **Extract Sequences**

   * Extracts the first 10 sequences from the input FASTA.

3. **Save Extracted Sequences**

   * Writes the extracted sequences to `influenza_extracted.fasta`.

4. **Generate BED File**

   * Converts the extracted sequences into BED format and saves to `influenza.bed`.

5. **Modify Sequences**

   * Computes the **reverse complement** of each extracted sequence.
   * Changes the **last nucleotide** of each sequence to `'A'`.

6. **Reinsert Modified Sequences**

   * Combines modified sequences with the remaining ones to create an updated FASTA dataset.

7. **Save Final Output**

   * Writes all modified sequences to `influenza_modified.fasta`.

---

## 📦 Example Usage

```bash
python handle_fasta.py
```

The script will:

* Extract 10 sequences from `Influenza.fna`
* Create `influenza_extracted.fasta` and `influenza.bed`
* Modify sequences (reverse complement + nucleotide change)
* Produce a final combined FASTA file `influenza_modified.fasta`

---

## 🧰 Requirements

* **Python ≥ 3.8**
* **Biopython**

Install dependencies:

```bash
pip install biopython
```

---

---

Would you like me to also draft a **short repository description** (the one-line tagline that appears under your repo title on GitHub)?
For example:

> “Python script for extracting, modifying, and converting Influenza FASTA sequences into BED format using Biopython.”
