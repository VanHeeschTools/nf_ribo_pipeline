#!/usr/bin/env python

import sys
import polars as pl
import re
import gffutils
from typing import List
from transcript_transformer.processing import csv_to_gtf, filter_CDS_variants

def main():

    # Step 1: Read input arguments
    h5_path = sys.argv[1]
    csv_files = sys.argv[2].split(',')
    min_samples = sys.argv[3]

    # Step 2: Load and concatenate all CSV files into a single Polars DataFrame
    print("Load and concatenate all CSV files")
    merged_df = load_and_concat_csvs(csv_files)
    save_dataframe(merged_df, "RiboTIE_unfiltered_merged.csv")

    # Step 3: Group by ORF_id and retain the row with the highest ribotie_score
    # Also filters to remove ORFs found in less than min_samples amount of samples
    print("Group and filter ORF dataframe")
    top_orfs_df = get_top_orf_scores(merged_df, min_samples)
    save_dataframe(top_orfs_df, "RiboTIE_duplicate_filtered_merged.csv")

    # Step 4: Filter CDS overlap
    filtered_top_orfs_df = filter_CDS_variants(top_orfs_df)

    # Step 5: Convert the result to GTF format using the HDF5 file
    print("Convert df to gtf file")
    convert_to_gtf(h5_path, filtered_top_orfs_df)
    merged_gtf = add_cds_summary_to_orfs("RiboTIE_merged.gtf", filtered_top_orfs_df)

    # Step 6: Save the filtered and merged output to a CSV file
    print("Save filtered df to csv file")
    save_dataframe(merged_gtf, "RiboTIE_merged.csv")

    print("Finished")

def load_and_concat_csvs(csv_files: List[str]) -> pl.DataFrame:
    """
    Loads multiple CSV files and concatenates them into a single Polars DataFrame.

    Parameters:
    csv_files : List[str]
        List of paths to CSV files.

    Returns:
    pl.DataFrame
        A concatenated DataFrame containing rows from all CSV files.
    """
    schema = {
        "orf_id": pl.Utf8,
        "seqname": pl.Utf8,       
        "start": pl.Int64,
        "end": pl.Int64,
        "strand": pl.Utf8,
        "score": pl.Float64,      
        "CDS_coords": pl.Utf8     
    }

    return pl.concat([
        pl.read_csv(f, schema_overrides=schema, infer_schema_length=1000, null_values=["."])
        for f in csv_files
    ])

def get_top_orf_scores(df: pl.DataFrame, min_samples: int = 1) -> pl.DataFrame:
    """
    Groups the DataFrame by 'ORF_id' and retains the row with the highest 'ribotie_score' per group.
    Only includes ORFs found in at least `min_samples` distinct samples.

    Parameters:
    df : pl.DataFrame
        Input DataFrame with columns 'ORF_id', 'ribotie_score', and 'sample_id'.

    min_samples : int
        Minimum number of distinct samples an ORF must appear in to be retained.

    Returns:
    pl.DataFrame
        A new DataFrame with the top-scoring row for each ORF_id that passes the sample count filter.
    """
    # Ensure ribotie_score is treated as float
    df = df.with_columns([
        pl.col("ribotie_score").cast(pl.Float64)
    ])

    # Count occurence of each ORF_id
    sample_counts = (
        df.with_columns(pl.col("ORF_id").cast(pl.Utf8))
        .group_by("ORF_id")
        .agg(pl.len().alias("sample_count"))
    )

    # Find ORFs with at least min_samples
    valid_orfs = sample_counts.filter(pl.col("sample_count") >= int(min_samples))["ORF_id"]

    # Filter out ORFs that don't reach occurence threshold
    df_filtered = df.filter(pl.col("ORF_id").is_in(valid_orfs))

    # Return the highest ribotie_score of each ORF_id
    highest_score = (
        df_filtered.sort("ribotie_score", descending=True)
            .group_by("ORF_id")
            .agg(pl.all().first())
    )

    # Filter out identical ORFs on different transcripts keep the ORF 
    # with the highest ribotie_score
    df_dedup = (
        highest_score.sort(["ribotie_score"], descending=[True])  # highest score first
        .group_by(["gene_id", "protein_seq", "TIS_coord", "TTS_coord"])
        .agg(pl.all().first())
    )

    # Return the highest ribotie_score of each ORF_id
    return (df_dedup)

def save_dataframe(df: pl.DataFrame, output_path: str) -> None:
    """
    Saves a Polars DataFrame to a CSV file.

    Parameters:
    df : pl.DataFrame
        The DataFrame to save.
    
    output_path : str
        Path to the output CSV file.
    
    Returns:
    None
        The function writes a csv file.
    """
    df.write_csv(output_path)
    print(f"Merged output written to {output_path}")

def convert_to_gtf(h5_path, top_orfs_df):
    """
    Converts the filtered ORF dataframe to GTF format using genomic features from an HDF5 file.

    Parameters:
    h5_path : str
        Path to the genomic features HDF5 database.

    top_orfs_df : pl.DataFrame or pd.DataFrame
        DataFrame containing ORFs, filtered to include only the top-scoring entries.

    Returns:
    None
        The function writes a GTF file using csv_to_gtf.
    """

    csv_to_gtf(
        h5_path,          # HDF5 file with genomic features
        top_orfs_df,      # Filtered ORF DataFrame
        "RiboTIE_merged", # Output filename prefix
        "RiboTIE"         # Source name for the GTF file
    )

def add_cds_summary_to_orfs(gtf_path: str, orf_df: pl.DataFrame) -> pl.DataFrame:
    """
    Appends a 'CDS_coords' column to the ORF DataFrame by summarizing
    all CDS regions per ORF_id from a GTF file, using gffutils for parsing.

    Parameters
    gtf_path : str
        Path to the GTF file with CDS lines
    orf_df : pl.DataFrame
        ORF DataFrame with 'ORF_id' column

    Returns
    pl.DataFrame
        Original ORF DataFrame with a new 'CDS_coords' column added
    """
    # Create a temporary gffutils DB (in-memory)
    db = gffutils.create_db(
        gtf_path,
        dbfn=":memory:",
        force=True,
        keep_order=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True
    )

    # Collect CDS coords grouped by ORF_id
    cds_map = {}
    for feature in db.features_of_type("CDS"):
        orf = feature.attributes.get("ORF_id", [None])[0]
        if orf is None:
            continue
        coord = f"{feature.seqid}:{feature.start}-{feature.end}({feature.strand})"
        key = feature.start if feature.strand == "+" else -feature.start
        cds_map.setdefault(orf, []).append((key, coord))

    # Create a summary DataFrame from sorted coordinate strings
    cds_df = pl.DataFrame(
        [(orf, "; ".join([c for _, c in sorted(coords)]))
        for orf, coords in cds_map.items()],
        schema=["ORF_id", "CDS_coords"]
    )

    # Join with original ORF dataframe
    return orf_df.join(cds_df, on="ORF_id", how="left")


if __name__ == "__main__":
    main()
