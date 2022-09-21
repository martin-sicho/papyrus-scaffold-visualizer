"""
environment necessary for running the app: drugex 
"""

import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem

from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.Scaffolds import MurckoScaffold

from rdkit.Chem import (
    AllChem,
    PandasTools,
    Draw,
    Descriptors,
    MACCSkeys,
    rdFingerprintGenerator,
    Scaffolds
)
from drugex.training.scorers.predictors import Predictor


import matplotlib.pyplot as plt
from sklearn.datasets import load_digits
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
import plotly.express as px
import molplotly



def loadfile(filename):
    '''
    Loads the .csv file and returns a dataframe
    '''
    return pd.read_csv(filename, sep='\t',)

def preparedata(data):
    """ 
    Prepares data for clustering

    data is a dataframe based on the .tsv file from papyrus

    returns a dataframe 
    """
    # get interesting data
    # descr_df = data[["SMILES", "InChIKey", "doc_id", "all_doc_ids", "Year", "pchembl_value_Mean"]]
    descr_df = data.copy()

    # add scaffolds and molecules to df
    descr_df["SMILES_Scaffold"] = descr_df.apply(lambda row: MurckoScaffold.MurckoScaffoldSmilesFromSmiles(row["SMILES"]), axis=1) # add MurckoScaffold in SMILES

    # convert smiles into rdkit molecules. 
    PandasTools.AddMoleculeColumnToFrame(descr_df, smilesCol="SMILES", molCol="Structure")
    PandasTools.AddMoleculeColumnToFrame(descr_df, smilesCol="SMILES_Scaffold", molCol="Structure_Scaffold")

    # add drugex fingerprint
    MorganVectorFP = Predictor.calc_ecfp(descr_df["Structure"]) # Morgan FP as bitvector
    # descr_df["MorganVectorFP"] = MorganVectorFP.tolist()

    # create df solely containing the fingerprints
    MorganFP_df = pd.DataFrame(MorganVectorFP)

    return descr_df, MorganFP_df

def calc_tsne(MorganFP_df, perplexity=40.0, verbose=False):
    '''
    Calculates the t-SNE embedding for the MorganFP_df
    '''
    tsne = TSNE( n_components=2, init='pca', verbose=verbose, perplexity=perplexity, learning_rate='auto') # what should I fill in here? What should the perplexity be?    
    # tsne = TSNE(n_components=n_components, verbose=1, perplexity=40, n_iter=300)
    tsne_results = tsne.fit_transform(MorganFP_df)
    tsne_df = pd.DataFrame(tsne_results, columns=["tSNE_1", "tSNE_2"])
    return tsne_df

def plot_tsne(tsne_df, descr_df, caption_cols, port):
    '''
    Plots the t-SNE embedding
    '''
    # plot t-SNE
    tsne_descr_morgan = pd.concat([descr_df, tsne_df], axis=1)

    series = pd.value_counts(tsne_descr_morgan.SMILES_Scaffold)
    mask = (series.lt(10))
    tsne_descr_morgan['scaffold_grouped'] = np.where(tsne_descr_morgan['SMILES_Scaffold'].isin(series[mask].index),'Other',tsne_descr_morgan['SMILES_Scaffold'])
    symbols = ['circle', 'square', 'diamond', 'cross', 'x',  'pentagon', 'hexagram', 'star',
           'diamond', 'hourglass', 'bowtie']

    # interactive visualization of clustering of scaffolds

    fig_tsne_morgan = px.scatter(tsne_descr_morgan, x="tSNE_1", y="tSNE_2",
                    color = "scaffold_grouped",symbol='scaffold_grouped',
                        symbol_sequence = symbols,                   
                        color_discrete_map= {'Other': 'lightgrey'},
                        title = 't-SNE on morganfingerprint data',
                                        width=1800,
                                        height=850, render_mode='SVG')
    fig_tsne_morgan.update_layout(plot_bgcolor='White')
    #fig_tsne_morgan.show()
    app_scatter = molplotly.add_molecules(fig=fig_tsne_morgan,
                                        df=tsne_descr_morgan,
                                        smiles_col=['SMILES', 'SMILES_Scaffold'],
                                        title_col='InChIKey',
                                        color_col = 'scaffold_grouped',
                                        caption_cols = caption_cols,
                                        width = 250)
                                                                            
    app_scatter.run_server(mode='inline', port=port, height=850)
    # when working on a server do not forget to forward the port, otherwise you won't see it in browser or inline. 
    #   in vscode: tab next to terminal "PORTS" add port 9401.

def cluster_generator(filename, port, perplexity=40.0, caption_cols = ["all_doc_ids", "Year", "pchembl_value_Mean"]):
    '''
    Generates the clusters for the dataframe
    '''
    data = loadfile(filename)
    descr_df, MorganFP_df = preparedata(data)
    tsne_df = calc_tsne(MorganFP_df, perplexity)
    plot_tsne(tsne_df, descr_df, caption_cols, port=port)





if __name__ == "__main__":
    filename = "/zfsdata/data/yorick/Internship/Files/hCCR2_LIGANDS.tsv"
    port = 9502
    cluster_generator(filename, port)