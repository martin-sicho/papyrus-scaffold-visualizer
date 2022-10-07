"""
example

Created by: Martin Sicho
On: 06.10.22, 14:17
"""
from data.papyrus.papyrus_class import Papyrus

if __name__ == "__main__":

    # acc_key = ["P32246", "P51681", "P51685",] # CCR1,5,8
    acc_keys = ["P41597"] # CCR2
    quality = "low" # choose from {"high", "medium", "low"}
    # name = "hCCR2_LIGANDS_stereo"
    name = "hCCR2_LIGANDS_nostereo"

    # call function
    papyrus = Papyrus(data_dir="./data", stereo=False)
    data_set = papyrus.getData(acc_keys, quality, name, drop_duplicates=False, use_existing=True)
    df = data_set.asDataFrame()
    print(df.columns)
    print(df.head())
    print(df.shape)
