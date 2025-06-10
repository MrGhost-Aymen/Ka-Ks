import os
import glob
import subprocess
import csv
import json
from Bio import AlignIO
from itertools import combinations
from datetime import datetime
import sys
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# === USER CONFIGURATION ===
KAKS_CALCULATOR_PATH = r"C:\Users\xxx\Desktop\Bioinfo_tool\KaKs_Calculator3.0\bin\KaKs.exe"
INPUT_DIR = r"C:\Users\xxx\Desktop\KA_KS\plastome\shared_genes_fasta\aligned"
OUTPUT_DIR = r"C:\Users\xxx\Desktop\KA_KS\plastome\results"
METHOD = "YN"  # Options: NG, YN, MYN, MLWL, MLPB, etc.
CODON_TABLE = "11"  # Mitochondrial codon table
MIN_SEQUENCES = 2  # Minimum number of sequences required in alignment
MAX_SEQUENCES = 1000  # Safety limit to prevent combinatorial explosion
SAVE_INTERMEDIATE_FILES = True  # Set to False to delete AXT files after processing
DEBUG_MODE = True  # Enable detailed debugging output

# Dictionary of available codon tables
CODON_TABLE_OPTIONS = {
    "0": "User-defined",
    "1": "Standard",
    "2": "Vertebrate Mitochondrial",
    "3": "Yeast Mitochondrial",
    "4": "Mold/Protozoan Mitochondrial",
    "5": "Invertebrate Mitochondrial",
    "6": "Ciliate Nuclear",
    "9": "Echinoderm Mitochondrial",
    "10": "Euplotid Nuclear",
    "11": "Bacterial",
    "12": "Alternative Yeast Nuclear",
    "13": "Ascidian Mitochondrial",
    "14": "Flatworm Mitochondrial",
    "15": "Blepharisma Macronuclear",
    "16": "Chlorophycean Mitochondrial",
    "21": "Trematode Mitochondrial",
    "22": "Scenedesmus obliquus Mitochondrial",
    "23": "Thraustochytrium Mitochondrial"
}

# Parse command-line arguments
import argparse
parser = argparse.ArgumentParser(description="Batch process aligned FASTA files for Ka/Ks analysis")
parser.add_argument("--kaks_path", default=KAKS_CALCULATOR_PATH, help="Path to KaKs_Calculator executable")
parser.add_argument("--input_dir", default=INPUT_DIR, help="Directory containing aligned FASTA files")
parser.add_argument("--output_dir", default=OUTPUT_DIR, help="Directory for output files")
parser.add_argument("--method", default=METHOD, help="KaKs calculation method (e.g., NG, YN, MYN)")
parser.add_argument("--codon_table", default=CODON_TABLE, help="Codon table ID (e.g., 2 for Vertebrate Mitochondrial)")
parser.add_argument("--min_sequences", type=int, default=MIN_SEQUENCES, help="Minimum sequences per alignment")
parser.add_argument("--max_sequences", type=int, default=MAX_SEQUENCES, help="Maximum sequences per alignment")
parser.add_argument("--save_intermediates", type=bool, default=SAVE_INTERMEDIATE_FILES, help="Save intermediate AXT/KaKs files")
parser.add_argument("--debug", type=bool, default=DEBUG_MODE, help="Enable debug logging")
args = parser.parse_args()

# Override defaults with CLI arguments
KAKS_CALCULATOR_PATH = args.kaks_path
INPUT_DIR = args.input_dir
OUTPUT_DIR = args.output_dir
METHOD = args.method
CODON_TABLE = args.codon_table
MIN_SEQUENCES = args.min_sequences
MAX_SEQUENCES = args.max_sequences
SAVE_INTERMEDIATE_FILES = args.save_intermediates
DEBUG_MODE = args.debug


def log_message(message):
    """Print timestamped messages with debug support."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")
    if DEBUG_MODE:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        with open(os.path.join(OUTPUT_DIR, "debug_log.txt"), 'a') as f:
            f.write(f"[{timestamp}] {message}\n")


def validate_environment():
    """Check if all required components are available."""
    try:
        assert os.path.exists(KAKS_CALCULATOR_PATH), f"KaKs_Calculator not found at {KAKS_CALCULATOR_PATH}"
        assert os.path.exists(INPUT_DIR), f"Input directory not found: {INPUT_DIR}"
        os.makedirs(OUTPUT_DIR, exist_ok=True)

        if CODON_TABLE not in CODON_TABLE_OPTIONS:
            raise ValueError(f"Invalid codon table: {CODON_TABLE}. Valid options: {list(CODON_TABLE_OPTIONS.keys())}")

        if "mito" in INPUT_DIR.lower() and CODON_TABLE not in ["2", "5", "9", "13", "14", "16", "21", "22", "23"]:
            log_message(f"Warning: Using codon table {CODON_TABLE} ({CODON_TABLE_OPTIONS[CODON_TABLE]}) for mitochondrial data")

        return True
    except Exception as e:
        log_message(f"Environment validation failed: {str(e)}")
        return False


def sanitize_filename(name):
    """Remove or replace characters invalid for file names."""
    invalid_chars = '<>:"/\\|?*'
    for char in invalid_chars:
        name = name.replace(char, '_')
    return name.strip()


def fasta_pair_to_axt(alignment, id1, id2, output_path):
    """Extract two sequences by ID and write AXT format with validation."""
    try:
        seq1 = next((rec.seq for rec in alignment if rec.id == id1), None)
        seq2 = next((rec.seq for rec in alignment if rec.id == id2), None)
        if not seq1 or not seq2:
            log_message(f"Missing sequences for pair {id1} vs {id2}")
            return False
        if len(seq1) != len(seq2):
            log_message(f"Length mismatch: {id1} ({len(seq1)}) vs {id2} ({len(seq2)})")
            return False

        safe_id1, safe_id2 = sanitize_filename(id1), sanitize_filename(id2)
        with open(output_path, 'w') as out:
            out.write(f"{safe_id1} {safe_id2}\n{seq1}\n{seq2}\n")
        log_message(f"Created AXT file: {output_path}")
        return True
    except Exception as e:
        log_message(f"Error creating AXT file {output_path}: {str(e)}")
        return False


def run_kaks(axt_path, output_path):
    """Run KaKs_Calculator on a single AXT file with robust error handling."""
    try:
        if not os.path.exists(axt_path):
            log_message(f"AXT file not found: {axt_path}")
            return False

        cmd = [KAKS_CALCULATOR_PATH, '-i', axt_path, '-o', output_path, '-m', METHOD, '-c', CODON_TABLE]
        log_message(f"Executing: {' '.join(cmd)}")

        creationflags = subprocess.CREATE_NO_WINDOW if sys.platform == "win32" else 0
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120, creationflags=creationflags)
        if result.returncode != 0:
            log_message(f"KaKs_Calculator failed for {axt_path}")
            log_message(f"Error output: {result.stderr[:500]}")
            return False

        log_message(f"Successfully processed: {output_path}")
        return True
    except subprocess.TimeoutExpired:
        log_message(f"Timeout running KaKs_Calculator on {axt_path}")
        return False
    except Exception as e:
        log_message(f"Unexpected error running KaKs: {str(e)}")
        return False


def parse_kaks_result(file_path):
    """Read and validate KaKs result file using headers for field mapping."""
    try:
        if not os.path.exists(file_path):
            log_message(f"Result file not found: {file_path}")
            return None

        with open(file_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        if len(lines) < 2:
            log_message(f"Incomplete results in {file_path}")
            return None

        header = lines[0].split('\t')
        data = lines[1].split('\t')

        result_map = dict(zip(header, data))

        def safe_float(key):
            val = result_map.get(key, 'NA')
            try:
                return float(val) if val not in ('NA', '') else float('nan')
            except ValueError:
                return float('nan')

        sequence_pair = result_map.get("Sequence", "")
        if " + " in sequence_pair:
            seq1, seq2 = sequence_pair.split(" + ")
        elif " " in sequence_pair:
            parts = sequence_pair.split()
            if len(parts) >= 2:
                seq1, seq2 = parts[0], parts[1]
            else:
                log_message(f"Invalid sequence pair format: {sequence_pair}")
                return None
        else:
            log_message(f"Invalid sequence pair format: {sequence_pair}")
            return None

        ka = safe_float("Ka")
        ks = safe_float("Ks")
        kaks = safe_float("Ka/Ks")
        p_gt = safe_float("P-Value(Fisher)")  # Assuming one-sided Fisher test
        p_lt = float('nan')  # Not present

        return [seq1, seq2, ka, ks, kaks, p_gt, p_lt]

    except Exception as e:
        log_message(f"Error parsing {file_path}: {str(e)}")
        return None


def process_fasta(fasta_path, out_dir):
    """Process a single FASTA file, generating all pairwise comparisons."""
    try:
        aln = AlignIO.read(fasta_path, "fasta")
        gene = os.path.splitext(os.path.basename(fasta_path))[0]
        ids = [rec.id for rec in aln]
        results = []

        if len(ids) < MIN_SEQUENCES:
            log_message(f"Skipping {gene}: Needs at least {MIN_SEQUENCES} sequences")
            return results
        if len(ids) > MAX_SEQUENCES:
            log_message(f"Skipping {gene}: Too many sequences ({len(ids)})")
            return results

        pair_count = len(list(combinations(ids, 2)))
        log_message(f"Processing {gene} with {len(ids)} sequences ({pair_count} pairs)")

        for i, (id1, id2) in enumerate(combinations(ids, 2), 1):
            pair_name = f"{gene}_{sanitize_filename(id1)}_vs_{sanitize_filename(id2)}"
            axt_path = os.path.join(out_dir, f"{pair_name}.axt")
            kaks_path = os.path.join(out_dir, f"{pair_name}.kaks")

            if os.path.exists(kaks_path) and os.path.getsize(kaks_path) > 0:
                log_message(f"Skipping already processed pair: {pair_name}")
                data = parse_kaks_result(kaks_path)
                if data:
                    results.append([gene] + data)
                continue

            log_message(f"Processing pair {i} of {pair_count}: {id1} vs {id2}")
            if fasta_pair_to_axt(aln, id1, id2, axt_path):
                if run_kaks(axt_path, kaks_path):
                    data = parse_kaks_result(kaks_path)
                    results.append([gene] + data)
                if not SAVE_INTERMEDIATE_FILES:
                    for f in [axt_path, kaks_path]:
                        try:
                            os.remove(f)
                            log_message(f"Removed intermediate file: {f}")
                        except Exception as e:
                            log_message(f"Error removing {f}: {str(e)}")
        return results
    except Exception as e:
        log_message(f"Error processing {fasta_path}: {str(e)}")
        return []


def process_fasta_wrapper(args):
    """Wrapper function to allow unpacking of arguments for multiprocessing"""
    fasta_path, out_dir = args
    return process_fasta(fasta_path, out_dir)


def batch_process():
    """Process all FASTA files in the input directory."""
    if not validate_environment():
        log_message("Exiting due to environment validation failure")
        return

    fasta_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.fasta"))) + \
                  sorted(glob.glob(os.path.join(INPUT_DIR, "*.fa")))
    if not fasta_files:
        log_message(f"No FASTA files found in {INPUT_DIR}")
        return

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    full_csv = os.path.join(OUTPUT_DIR, f"kaks_summary_full_{timestamp}.csv")
    clean_csv = os.path.join(OUTPUT_DIR, f"kaks_summary_clean_{timestamp}.csv")
    json_path = os.path.join(OUTPUT_DIR, f"kaks_summary_{timestamp}.json")
    stats_csv = os.path.join(OUTPUT_DIR, f"gene_kaks_stats_{timestamp}.csv")
    plot_file = os.path.join(OUTPUT_DIR, f"ka_ks_barplot_{timestamp}.png")

    log_message(f"Starting Ka/Ks analysis with method: {METHOD}")
    log_message(f"Using codon table: {CODON_TABLE} ({CODON_TABLE_OPTIONS[CODON_TABLE]})")
    log_message(f"Found {len(fasta_files)} FASTA files to process")
    log_message(f"Results will be saved to: {full_csv}, {clean_csv}, {stats_csv}, {plot_file}")

    try:
        num_processes = min(cpu_count(), len(fasta_files))
        log_message(f"Using {num_processes} processes for parallel execution")

        with Pool(processes=num_processes) as pool:
            results = list(tqdm(pool.imap_unordered(
                process_fasta_wrapper,
                [(f, OUTPUT_DIR) for f in fasta_files]
            ), total=len(fasta_files)))

        flat_results = [item for sublist in results for item in sublist]
        log_message(f"Collected {len(flat_results)} total results across all genes")

        # Convert to DataFrame
        df = pd.DataFrame(flat_results, columns=[
            "Gene", "Sequence1", "Sequence2", "Ka", "Ks", "Ka/Ks", "P_value(Ka>Ks)", "P_value(Ka<Ks)"
        ])
        df["Codon_Table"] = CODON_TABLE_OPTIONS[CODON_TABLE]

        # Save full report (including NaN values)
        df.to_csv(full_csv, index=False)
        log_message(f"Saved full summary to: {full_csv}")

        # Clean report (only valid rows)
        clean_df = df[df['Ka'].notna() & df['Ks'].notna()]
        clean_df = clean_df[clean_df['Ks'] > 0]  # Avoid division by zero
        clean_df.to_csv(clean_csv, index=False)
        log_message(f"Saved clean summary to: {clean_csv}")

        # Gene-level statistics
        if not clean_df.empty:
            gene_stats = clean_df.groupby('Gene')['Ka/Ks'].agg(['mean', 'median', 'std', 'count'])
            gene_stats = gene_stats.sort_values(by='mean', ascending=False).head(20)
            gene_stats.to_csv(stats_csv)
            log_message(f"Saved gene statistics to: {stats_csv}")

            # Plotting
            plt.figure(figsize=(12, 6))
            sns.barplot(x=gene_stats.index, y=gene_stats['mean'], palette="viridis", ci=None)
            plt.xticks(rotation=45, ha='right')
            plt.ylabel("Average Ka/Ks")
            plt.title("Average Ka/Ks Values Across Genes")
            plt.tight_layout()
            plt.savefig(plot_file, dpi=300)
            log_message(f"Saved Ka/Ks bar plot to: {plot_file}")

        log_message(f"Success! Processed {len(flat_results)} sequence pairs")
        log_message(f"Summary saved to: {full_csv}, {clean_csv}, {stats_csv}, {plot_file}")

    except Exception as e:
        log_message(f"Fatal error during batch processing: {str(e)}")


if __name__ == "__main__":
    start_time = datetime.now()
    log_message("=== Script started ===")
    log_message(f"Configuration:")
    log_message(f"- Method: {METHOD}")
    log_message(f"- Codon Table: {CODON_TABLE} ({CODON_TABLE_OPTIONS[CODON_TABLE]})")
    log_message(f"- Min sequences: {MIN_SEQUENCES}")
    log_message(f"- Max sequences: {MAX_SEQUENCES}")
    log_message(f"- Save intermediates: {SAVE_INTERMEDIATE_FILES}")
    try:
        batch_process()
    except Exception as e:
        log_message(f"CRITICAL ERROR: {str(e)}")
    finally:
        runtime = datetime.now() - start_time
        log_message(f"Total runtime: {runtime}")
        log_message("=== Script completed ===")