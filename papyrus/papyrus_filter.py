import pandas as pd
import papyrus_scripts
from papyrus_scripts.reader import read_papyrus, read_protein_set 
from papyrus_scripts.preprocess import (keep_quality, keep_source, keep_type,
                                        keep_organism, keep_accession, keep_protein_class,
                                        keep_match, keep_contains
                                       )
from papyrus_scripts.preprocess import consume_chunks
import math



def papyrus_filter(acc_key: list, chunksize: int, quality: str, name: str):
    '''
    This script should filter the downloaded papyrus dataset for quality and accession key (UniProt) and output a .tsv file of all compounds fullfilling these requirements.
    args:
        acc_key: list of UniProt accession keys
        chunksize: integer of chunks to process one at the time
        quality: str with quality of dataset
        name: string to name the file to save
    Output: .tsv file
    '''
    # read data
    sample_data = read_papyrus(is3d=False, chunksize=1000_000, source_path="/zfsdata/data/yorick/DrugEx_Papyrus/")
    print("Read all data.")

    # data filters
    filter1 = keep_quality(data=sample_data, min_quality=quality)
    filter2 = keep_accession(data=filter1, accession=acc_key)
    print("Initialized filters.")

    # filter data per chunk
    tot_chunks = math.ceil(60_000_000 / chunksize) # amount of chunks
    print(f"Starting to filter data for {tot_chunks} chunks of size {chunksize}.")
    filtered_data = consume_chunks(generator=filter2, progress=True, total=tot_chunks)
    print(f"Number of compounds:{filtered_data.shape[0]}")

    # filter out duplicate InChiKeys
    print("Filtering out duplicate molecules")
    amnt_mols_i = len(filtered_data["InChIKey"])
    filtered_data.drop_duplicates(subset=["InChIKey"], inplace=True, ignore_index=True)
    amnt_mols_f = len(filtered_data["InChIKey"])
    print(f"Filtered out {amnt_mols_i - amnt_mols_f} duplicate molecules")

    # write filtered data to .tsv file
    filtered_data.to_csv(f"/zfsdata/data/yorick/Internship/Files/{name}.tsv", sep= "\t", index=False)
    print(f"Wrote data to file {name}.tsv")




if __name__ == "__main__":
 
    # variables

    # acc_key = ["P32246", "P51681", "P51685",] # CCR1,5,8 
    acc_key = ["P41597"] # CCR2
    chunksize = 1000_000
    quality = "low" # choose from {"high", "medium", "low"}
    name = "hCCR2_LIGANDS"

    # call function
    papyrus_filter(acc_key, chunksize, quality, name)



