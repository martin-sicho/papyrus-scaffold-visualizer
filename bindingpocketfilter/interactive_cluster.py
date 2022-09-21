import pandas as pd
import numpy as np
import plotly.express as px
import molplotly


from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import PandasTools
from sklearn.manifold import TSNE
from drugex.training.scorers.predictors import Predictor

class InteractiveCluster:

    def __init__(self, perplexity = 40.0, n_iter = 1000, n_components = 2, init = "pca", verbose = True,):
        """
        Initializes the InteractiveCluster class

        parameters:
            perplexity: perplexity parameter for t-SNE
            n_iter: number of iterations for t-SNE
            n_components: number of components for t-SNE
            init: initialization method for t-SNE
            verbose: verbosity for t-SNE

        """
        
        self.perplexity = perplexity
        self.n_iter = n_iter
        self.n_components = n_components
        self.init = init
        self.verbose = verbose

        # define unset class variables
        self.data = None
        self.cols = None
        self.descr = None
        self.ecfp6 = None
        self.tsne = None
    
    def load_data(self, data, fromFile = True, cols =["InChIKey", "all_doc_ids", "Year", "pchembl_value_Mean", "database"]):
        '''
        Loads the dataframe from a file if fromFile = True, otherwise loads the dataframe from the data df

        parameters:
            data: either a file (fromFile = True) or a dataframe (fromFile = False) containing column "SMILES"
            fromfile: boolean indicating whether the data is from a file or df
            cols: columns to be loaded in the dataframe
        returns:
            dataframe with the columns specified in cols as self.data
        '''
        self.cols = cols
        if fromFile:
            self.data = pd.read_csv(data, sep='\t', keep_default_na=False, usecols=["SMILES"] + cols) 
        else:
            self.data = data.copy()
            
    def load_prepdata(self, filename, dfs = ["descr", "ecfp6", "tsne"]):
        '''
        Loads the preprocessed data in the specified dataframes from a file

        parameters:
            filename: name of the file(s) to be loaded without the _{df}.tsv extension 
        '''
        for df in dfs:
            try:
                self.__dict__[df] = pd.read_csv(filename + "_" + df + ".tsv", sep='\t', keep_default_na=False) # keep_default_na=False is necessary so that empty strings are not converted to NaN values, which are not processed by PandasTools
                if df == "descr":
                    # add mols to the descr df
                    PandasTools.AddMoleculeColumnToFrame(self.descr, smilesCol="SMILES", molCol="Structure")
                    PandasTools.AddMoleculeColumnToFrame(self.descr, smilesCol="SMILES_Scaffold", molCol="Structure_Scaffold")
            except:
                print("{} was not able to load.".format(df))
        
    
    def save(self, filename, dfs = ["descr", "ecfp6", "tsne"]):
        '''
        Saves the specified dataframes to a file

        parameters:
            filename: name of the file(s) to be saved without the _{df}.tsv extension
        '''
        for df in dfs:
            try:
                if df == "descr":
                    cp = self.descr.copy()
                    cp.drop(["Structure", "Structure_Scaffold"], axis=1).to_csv(filename + "_descr.tsv", sep='\t', index=False) # drop the rdkit molecules and save the df without them
                else:
                    self.__dict__[df].to_csv(filename + "_" + df + ".tsv", sep='\t', index=False)
            except:
                print("{} was not able to save.".format(df))

    
    def prep_data(self):
        """ 
        Prepares RAW data in self.data loaded in load_data() for clustering
        Calculates MurckoScaffold, and (drugex) efcp6 fingerprints
        Adds rdkit molecule objects from SMILES and MurckoScaffold       
        """
        # rename data
        descr= self.data

        # add scaffolds and molecules to df
        descr["SMILES_Scaffold"] = descr.apply(lambda row: MurckoScaffold.MurckoScaffoldSmilesFromSmiles(row["SMILES"]), axis=1) # add MurckoScaffold in SMILES

        # convert smiles into rdkit molecules. 
        PandasTools.AddMoleculeColumnToFrame(descr, smilesCol="SMILES", molCol="Structure")
        PandasTools.AddMoleculeColumnToFrame(descr, smilesCol="SMILES_Scaffold", molCol="Structure_Scaffold")

        # calculate drugex fingerprint
        MorganVectorFP = Predictor.calc_ecfp(descr["Structure"]) # Morgan FP as bitvector
        
        # create df solely containing the fingerprints
        ecfp6 = pd.DataFrame(MorganVectorFP)

        # save dataframes
        self.descr = descr
        self.ecfp6 = ecfp6

    
    def calc_tsne(self):
        '''
        Calculates the t-SNE embedding for the ecfp6

        parameters: None, have been initialized in __init__
            #TODO: move parameters for perplexity, n_iter, n_components, init, verbose to here
        '''
        tsne_env = TSNE(n_components=self.n_components, init=self.init, verbose=self.verbose, perplexity=self.perplexity, learning_rate='auto') 
        # tsne = TSNE(n_components=n_components, verbose=1, perplexity=40, n_iter=300)
        tsne_results = tsne_env.fit_transform(self.ecfp6)
        tsne = pd.DataFrame(tsne_results, columns=["tSNE_1", "tSNE_2"])

        # save tsne
        self.tsne = tsne
    
    def group_scaffold(self, df):
        '''
        Groups scaffold, per scaffold
        if a scaffold is > 10 times present in df, then the scaffold is added to scaffold_grouped, otherwise 'other' is added.

        parameters:
            df: dataframe containing column 'SMILES_Scaffold'
        returns:
            df, with extra column 'Scaffold_group'
        '''
        series = pd.value_counts(df.SMILES_Scaffold)
        mask = (series.lt(10))
        df['scaffold_grouped'] = np.where(df['SMILES_Scaffold'].isin(series[mask].index),'Other',df['SMILES_Scaffold']) # if the SMILES_Scaffold is in the mask, then put it in the 'Other' group, otherwise put it in the SMILES_Scaffold group
        return df

    def plot_tsne(self, port, color_by = "Scaffold", title_col = "InChIKey", caption_cols = ["Year", "pchembl_value_Mean", "all_doc_ids"]):
        '''
        Plots the t-SNE embedding

        parameters:
            port (int): the port to plot the t-SNE embedding on
            color_by (str): the column to color the points by, default is Scaffold, any other column you should add yourself.
            title_col (str): the column to use as the title, default is InChIKey
            caption_cols (list): the columns to use in the caption, default is ["Year", "pchembl_value_Mean", "all_doc_ids"]

        TODO: use dash app to plot these points and be able to interactively change colormapping
        '''
        # assert caption_cols in self.descr.columns, f"caption_cols {caption_cols} must be in the columns of the descr"

        symbols = ['circle', 'square', 'diamond', 'cross', 'x',  'pentagon', 'hexagram', 'star', 'diamond', 'hourglass', 'bowtie']
        tsne_descr = pd.concat([self.descr, self.tsne], axis=1)

        # create scatter plot by either scaffold coloring or color by column. 
        # TODO:
        #   Generalize fig_tsne generation by using *kwargs.

        if color_by == "Scaffold":
            tsne_descr = self.group_scaffold(tsne_descr)
            color_by = "scaffold_grouped"
            color_discrete_map = {'Other': 'lightgrey'}
            fig_tsne = px.scatter(tsne_descr, x="tSNE_1", y="tSNE_2",
                                                    color = color_by, symbol=color_by,
                                                    symbol_sequence = symbols,                   
                                                    color_discrete_map = color_discrete_map,
                                                    title = 't-SNE on ECFP6 data',
                                                    width=1800,
                                                    height=850, render_mode='SVG'
                                                )
            fig_tsne.update_layout(plot_bgcolor='White')
        else: # color by column color_by
            fig_tsne = px.scatter(tsne_descr, x="tSNE_1", y="tSNE_2",
                                                    color = color_by,symbol=color_by,
                                                    symbol_sequence = symbols,
                                                    title = 't-SNE on ECFP6 data',
                                                    width=1800,
                                                    height=850, render_mode='SVG'
                                                )
            fig_tsne.update_layout(plot_bgcolor='White')
        # fig_tsne.show()
        
        # interactive plot:
        app_scatter = molplotly.add_molecules(fig=fig_tsne,
                                                            df=tsne_descr,
                                                            smiles_col=['SMILES', 'SMILES_Scaffold'],
                                                            title_col= title_col,
                                                            color_col = color_by,
                                                            caption_cols = caption_cols,
                                                            width = 250
                                                            )
                                                                                
        app_scatter.run_server(mode='inline', port=port, height=850)
        # when working on a server do not forget to forward the port, otherwise you won't see it in browser or inline. 
        #   in vscode: tab next to terminal "PORTS" add port.

    def gen_tsne(self, pfile, unprep=True, load_cols = ["InChIKey", "all_doc_ids", "Year", "pchembl_value_Mean", "database"], 
                    caption_col = ["Year", "pchembl_value_Mean", "all_doc_ids"], port=9401, title_col = "InChIKey", color_by = "Scaffold"):
        '''
        Goes through the workflow to generate the t-SNE plots from data.

        If unprep is True, then the file provided is loaded and preprocessed, and the t-SNE is calculated.
            the generated dfs are saved in the same directory as the file provided and the tsne is plotted.
            The file provided must contain the columns: "SMILES", "InChIKey",
        If unprep is False, then the files in pfile_{descr/ecfp6/tsne}.tsv are loaded and the tsne is plotted.

        When runnning default paramters, the loaded data is expected to have columns: "SMILES", "InChIKey", "all_doc_ids", "Year", "pchembl_value_Mean", "database"

        parameters:
            pfile (str): the path to the file(s) to load
            unprep (bool): whether file(s) are preprocessed or not
            load_cols (list): the columns to load from the file(s)
            caption_col (list): the columns to use in the caption
            port (int): the port to plot the t-SNE embedding on
            title_col (str): the column to use as the title, default is InChIKey
            color_by (str): the column to color the points by, default is Scaffold, can be any column that has been loaded.
        output:
            interactive plot of the t-SNE embedding at the specified port.
        '''
        
        # unpreprocessed data:
        if unprep:
            self.load_data(data = pfile, cols = load_cols)
            self.prep_data()
            self.calc_tsne()        
            self.save(pfile.removesuffix('.tsv'))
        # preprocessed data:
        else:
            self.load_prepdata(pfile)

        self.plot_tsne(port=port, color_by=color_by, title_col=title_col, caption_cols=caption_col)
        