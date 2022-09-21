import pandas as pd


def datasplitter(df: pd.DataFrame, col_name: str, threshold: float, up = True, down = True) -> pd.DataFrame:
    '''
    Splits data from threshold to above and lower than threshold. The splits are returned as DataFrames.
    
    args:
        df: panda DataFrame
        col_name: column to select for threshold
        threshold: integer or float as threshold
        up: boolean when you want to save the upper subset
        down: boolean when you want to save lower subset

    returns: dataframes representing subsets 'up' if True and 'down' if True.
    
    '''
    if up:
        df_up = df[df[col_name].astype(float) >= threshold] # cast type of col_name from string into float
    if down:
        df_down = df[df[col_name].astype(float) < threshold]
    
    if up and down:
        return df_up, df_down # returns tuple
    elif up:
        return df_up
    else:
        return df_down

     


def main():
    '''
    # 1. read dataset .tsv file
    # 2. filter upper part threshhold and write to new df
    # 3. filter lower part threshold and write to new df
    # 4. save new df as .tsv files
    
    '''

    # Variables
    file = "/zfsdata/data/yorick/Internship/Files/hCCR158_LIGANDS.tsv"
    col_name = 'pchembl_value_Mean'
    # threshold = 6.5 

    # read .tsv file into dataframe
    df = pd.read_csv(file, sep='\t')
    data = df["pchembl_value_Mean"].astype(float) # cast pchembl as float to calculate median to set as threshold
    threshold = data.median()

    # Call function
    up, down = datasplitter(df, col_name, threshold)

    #save new files
    up.to_csv("/zfsdata/data/yorick/Internship/Files/hCCR158_LIGANDS_active_median.tsv", sep='\t', index=False)
    down.to_csv("/zfsdata/data/yorick/Internship/Files/hCCR158_LIGANDS_inactive_median.tsv", sep='\t', index=False)
    
    
if __name__ == "__main__":
    main()