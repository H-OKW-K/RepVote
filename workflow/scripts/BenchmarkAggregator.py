import pandas as pd
import glob
import os

def aggregate_benchmarks(root_dir="benchmarks", output_file="total_benchmark_summary.tsv"):
    # 1. Find all .tsv files recursively
    # recursive=True allows searching in subdirectories of 'benchmarks'
    file_pattern = os.path.join(root_dir, "**", "*.tsv")
    tsv_files = glob.glob(file_pattern, recursive=True)
    
    if not tsv_files:
        print(f"No .tsv files found in directory: {root_dir}")
        return

    print(f"Found {len(tsv_files)} benchmark files. Reading...")

    dataframes = []
    
    for file_path in tsv_files:
        try:
            # Read the TSV file
            df = pd.read_csv(file_path, sep='\t')
            
            # Add a column for the filename so you know which job produced the metric
            df['filename'] = os.path.basename(file_path)
            
            # Optional: extract rule name if your filenames follow a pattern like rule.sample.tsv
            # df['rule_guess'] = os.path.basename(file_path).split('.')[0]
            
            dataframes.append(df)
        except Exception as e:
            print(f"Warning: Could not read {file_path}. Reason: {e}")

    if not dataframes:
        print("No valid dataframes to concatenate.")
        return

    # 2. Concatenate all dataframes into one
    combined_df = pd.concat(dataframes, ignore_index=True)

    # 3. Sort by Memory Usage (max_rss) to show peak load at the top
    # Snakemake usually reports max_rss in MB
    if 'max_rss' in combined_df.columns:
        combined_df = combined_df.sort_values(by='max_rss', ascending=False)

    # 4. Save to a single TSV file
    combined_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSuccess! Aggregated benchmarks saved to: {output_file}")

    # 5. Display Peak Load Summary
    print("\n--- Top 5 Peak Memory (max_rss) Jobs ---")
    if 'max_rss' in combined_df.columns:
        # Select relevant columns to display
        display_cols = ['filename', 'max_rss', 's']
        # Filter for existing columns only
        cols_to_show = [c for c in display_cols if c in combined_df.columns]
        
        print(combined_df[cols_to_show].head(5).to_string(index=False))
    else:
        print("Column 'max_rss' not found in benchmark files.")

if __name__ == "__main__":
    aggregate_benchmarks(root_dir=r"/home/cellularsenescence/Downloads/Okw_vol/RepVote/benchmarks")

