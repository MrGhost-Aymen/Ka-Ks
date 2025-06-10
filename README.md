Ka/Ks Batch Processing Script
Overview
This Python script automates the calculation of Ka/Ks
Prerequisites

Python 3.6+
Required Libraries:
biopython - For parsing FASTA alignment files
pandas - For data manipulation and CSV output
seaborn - For visualization
matplotlib - For plotting Ka/Ks results
tqdm - For progress bars during processing


KaKs_Calculator3.0 - External tool for Ka/Ks computation (place KaKs.exe in a specified path)
Install dependencies via:pip install biopython pandas seaborn matplotlib tqdm



Installation

Clone this repository:git clone https://github.com/yourusername/kaks-batch-processor.git


Place the KaKs_Calculator3.0 executable in a directory (e.g., bin/ within the project).
Update the KAKS_CALCULATOR_PATH variable in the script to point to the executable.

Usage
Run the script from the command line with optional arguments to customize paths, methods, and settings:
python kaks_batch.py --input_dir "path/to/aligned/fasta" --output_dir "path/to/results" --method YN --codon_table 11

Command-Line Arguments

--kaks_path: Path to KaKs_Calculator executable (default: C:\Users\moi\Desktop\Bioinfo_tool\KaKs_Calculator3.0\bin\KaKs.exe)
--input_dir: Directory containing aligned FASTA files (.fasta or .fa)
--output_dir: Directory for output files (CSV, JSON, plots)
--method: Ka/Ks calculation method (e.g., NG, YN, MYN, default: YN)
--codon_table: Codon table ID (e.g., 11 for Bacterial, 2 for Vertebrate Mitochondrial, default: 11)
--min_sequences: Minimum sequences per alignment (default: 2)
--max_sequences: Maximum sequences per alignment (default: 1000)
--save_intermediates: Save intermediate AXT/KaKs files (default: True)
--debug: Enable detailed debug logging (default: True)

Outputs

Full Summary CSV (kaks_summary_full_YYYYMMDD_HHMMSS.csv): All results, including invalid/NaN values.
Clean Summary CSV (kaks_summary_clean_YYYYMMDD_HHMMSS.csv): Filtered results with valid Ka and Ks values.
Gene Statistics CSV (gene_kaks_stats_YYYYMMDD_HHMMSS.csv): Mean, median, std, and count of Ka/Ks per gene (top 20).
Bar Plot (ka_ks_barplot_YYYYMMDD_HHMMSS.png): Visualization of average Ka/Ks values for top 20 genes.
Debug Log (debug_log.txt): Timestamped logs, enhanced if --debug is True.

Features

Batch Processing: Handles multiple aligned FASTA files in parallel using multiprocessing.
Pairwise Analysis: Computes Ka/Ks for all sequence pairs within each alignment.
Error Handling: Robust validation, logging, and timeout for KaKs_Calculator execution.
Flexibility: Configurable codon tables, methods, and sequence limits.
Visualization: Bar plot of average Ka/Ks values per gene using Seaborn.

Notes

Ensure aligned FASTA files are properly formatted and sequences are of equal length.
Mitochondrial data may require specific codon tables (e.g., 2 for vertebrates).
Intermediate files (AXT, KaKs output) are saved if --save_intermediates is True, otherwise deleted.
Adjust --max_sequences to avoid combinatorial explosion for large alignments.

Example
python kaks_batch.py --input_dir "C:\Users\xxx\Desktop\KA_KS\plastome\shared_genes_fasta\aligned" --output_dir "C:\Users\xxx\Desktop\KA_KS\plastome\results" --method YN --codon_table 11

Troubleshooting

KaKs_Calculator not found: Verify the executable path and update KAKS_CALCULATOR_PATH.
No FASTA files found: Check the --input_dir path for .fasta or .fa files.
Invalid codon table: Use a value from the supported list (e.g., 1 for Standard, 2 for Vertebrate Mitochondrial).
Check debug_log.txt for detailed error messages if --debug is enabled.

License
MIT License - feel free to use, modify, and distribute this script.
Contributing
Pull requests and issues are welcome! Please test changes and update this README if needed.
Contact
For questions, open an issue on GitHub or contact ouamoa@gmail.com.
