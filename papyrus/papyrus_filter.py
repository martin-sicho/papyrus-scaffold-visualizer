import os.path

from papyrus_scripts.reader import read_papyrus
from papyrus_scripts.preprocess import keep_quality, keep_accession
from papyrus_scripts.preprocess import consume_chunks

from papyrus import OUTDIR


def papyrus_filter(acc_key: list, quality: str, drop_duplicates: bool = True, chunk_size : int = 1e5, outfile: str = None):
    """
    Filters the downloaded papyrus dataset for quality and accession key (UniProt) and outputs a .tsv file of all compounds fulfilling these requirements.

    Args:
        outfile: path to the output file
        acc_key: list of UniProt accession keys
        quality: str with minimum quality of dataset to keep
        drop_duplicates: boolean to drop duplicates from the final dataset
        chunk_size: integer of chunks to process one at the time
    Output:
        .tsv file with all compounds fulfilling the requirements
    """
    # read data
    sample_data = read_papyrus(is3d=False, chunksize=chunk_size, source_path=OUTDIR)
    print("Read all data.")

    # data filters
    filter1 = keep_quality(data=sample_data, min_quality=quality)
    filter2 = keep_accession(data=filter1, accession=acc_key)
    print("Initialized filters.")

    # filter data per chunk
    filtered_data = consume_chunks(generator=filter2)
    print(f"Number of compounds:{filtered_data.shape[0]}")

    # filter out duplicate InChiKeys
    if drop_duplicates:
        print("Filtering out duplicate molecules")
        amnt_mols_i = len(filtered_data["InChIKey"])
        filtered_data.drop_duplicates(subset=["InChIKey"], inplace=True, ignore_index=True)
        amnt_mols_f = len(filtered_data["InChIKey"])
        print(f"Filtered out {amnt_mols_i - amnt_mols_f} duplicate molecules")

    # write filtered data to .tsv file
    if outfile:
        filtered_data.to_csv(outfile, sep= "\t", index=False)
        print(f"Wrote data to file {name}.tsv")

    return filtered_data

if __name__ == "__main__":
 
    # variables

    # acc_key = ["P32246", "P51681", "P51685",] # CCR1,5,8 
    acc_key = ["P41597"] # CCR2
    quality = "low" # choose from {"high", "medium", "low"}
    name = "hCCR2_LIGANDS"

    # call function
    papyrus_filter(acc_key, quality, drop_duplicates=False, chunk_size=1000_000, outfile=os.path.join(OUTDIR, f"./{name}.tsv"))



