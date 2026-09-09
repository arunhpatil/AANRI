from pathlib import Path
import pandas as pd

# --------------------------------------------------
# Paths
# --------------------------------------------------

input_dir = Path("/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gwas_SS_fixed")
#input_dir = Path("/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gwas_ss")
output_dir = Path("/dcs04/lieber/hwanglab/Arun/snRNA_aanri/0_annotation/gwas_ss_chrSpecific")

# Create chromosome directories
for chr_no in range(1, 23):
    (output_dir / f"chr{chr_no}").mkdir(
        parents=True,
        exist_ok=True
    )


# --------------------------------------------------
# Process each summary statistics file
# --------------------------------------------------

for input_file in input_dir.iterdir():

    if not input_file.is_file():
        continue
    
    if input_file.name == "SORTED_PTSD_AA7_ALL_study_specific_PCs1.txt":
        print(f"{input_file.name}: This file was ignored, because it doesn't have chromosome coordinates\n")
        continue
    print(f"\nProcessing: {input_file.name}")

    # Keep the same filename in the chromosome-specific folders
    filename = input_file.name

    total_removed = 0
    total_rs_replaced = 0
    # Read file in chunks to avoid loading huge GWAS files
    for chunk in pd.read_csv(
        input_file,
        sep="\t",
        chunksize=500_000
    ):

        # Check that CHR exists
        if "CHR" not in chunk.columns:
            raise ValueError(
                f"'CHR' column not found in {input_file}"
            )

        # ---------------------------------------------------------
        # Remove defective GWAS records
        # ---------------------------------------------------------
        defective = (
            chunk["SNP"].isna()
            | chunk["SNP"].astype(str).eq(".")
        )

        n_defective = defective.sum()
        if n_defective > 0:
            print(
                f"  Removing {n_defective} defective rows "
                f"from chunk"
            )

        chunk = chunk.loc[~defective].copy()

        if chunk.empty:
            continue

        # ---------------------------------------------------------
        # Normalize SNP IDs
        #
        # If SNP contains an rsID anywhere in the field,
        # replace the entire field with that rsID.
        #
        # Example:
        # chr22:16269838:G:A;chr22:16269838;rs1175530675
        #                     ->
        # rs1175530675
        # ---------------------------------------------------------
        snp_before = chunk["SNP"].copy()

        chunk["SNP"] = (
            chunk["SNP"]
            .str.extract(r"(rs\d+)", expand=False)
            .fillna(chunk["SNP"])
        )

        n_rs_replaced = (
            snp_before != chunk["SNP"]
        ).sum()

        if n_rs_replaced > 0:
            print(
                f"  Replaced {n_rs_replaced} complex SNP IDs "
                f"with rsIDs"
            )

        total_rs_replaced += n_rs_replaced

        # Split by chromosome
        for chr_no in range(1, 23):

            chr_name = f"chr{chr_no}"

            chr_chunk = chunk[
                chunk["CHR"] == chr_name
            ]

            if chr_chunk.empty:
                continue

            output_file = (
                output_dir /
                chr_name /
                filename
            )

            # Write header only for the first chunk
            write_header = not output_file.exists()

            chr_chunk.to_csv(
                output_file,
                sep="\t",
                index=False,
                mode="a",
                header=write_header
            )

    print(f"Completed: {input_file.name}")

print("\nAll files processed.")
