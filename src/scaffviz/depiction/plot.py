"""
plot

Created by: Martin Sicho
On: 05.10.22, 16:37
"""
import copy
import threading
import time

import molplotly
import pandas as pd
import plotly.express as px
from qsprpred.data.utils.descriptorcalculator import CustomDescriptorsCalculator
from qsprpred.data.utils.descriptorsets import DataFrameDescriptorSet

from scaffviz.clustering.manifold import Manifold
from qsprpred.data.data import MoleculeTable, QSPRDataset
from scaffviz.data.manifold_table import ManifoldTable
from qsprpred.models.models import QSPRModel
from qsprpred.plotting.interfaces import ModelPlot
from typing import List, Literal
from qsprpred.models.tasks import ModelTasks


class Plot:

    def __init__(self, manifold : Manifold, save_manifold : bool = True):
        """
        Initialize a plotting object for the given `Manifold`.

        Args:
            manifold: the `Manifold` class to use to project molecules to 2D
            save_manifold: if `True` the calculated 2D coordinates are saved to the `MoleculeTable` object
        """

        self.symbols = ['circle', 'square', 'diamond', 'cross', 'x',  'pentagon', 'hexagram', 'star', 'diamond', 'hourglass', 'bowtie']
        self.open_apps = dict()
        self.save_manifold = save_manifold
        self.manifold = manifold

    def getOpenApps(self):
        return self.open_apps

    def plot(self, table : MoleculeTable, x : str = None, y : str = None, color_by : str = None, card_data = tuple(), title_data : str = 'SMILES', port=9292, recalculate=False, mols_per_scaffold_group : int = 10, interactive = True, viewport_height = "100%",  **kwargs):
        """
        Plot the dataset using the manifold or custom `DataSet` fields. The plot is interactive and runs as a web app on the specified port.

        Args:
            table: the `MoleculeTable` object to plot molecules and data from
            x: the name of the variable in the data set to use for the x-axis, if not specified the first dimension of the manifold is used
            y: the name of the variable in the data set to use for the y-axis, if not specified the second dimension of the manifold is used
            color_by: the data to color the points by, by default the first scaffold found in the `DataSet` will be used
            card_data: `list` of data names from the `DataSet` to show on the cards displayed when hovering over a molecule in the interactive plot, ignored if `interactive` is `False`
            title_data: the data to get from the `DataSet` as the card title, ignored if `interactive` is `False`
            port: port to run the interactive web app on, ignored if `interactive` is `False`
            recalculate: whether to recalculate the manifold or use the existing data in the dataset
            mols_per_scaffold_group: how many molecules to include in one scaffold group, only applicable if `color_by` is not specified, the scaffolds with the number of molecules lower than this value will be shown in grey in the plot
            interactive: whether to run the plot as an interactive web app or just return the figure object
            viewport_height: height of the viewport in the browser (use this ie. to make the iframe containing the plot bigger), applies only to interactive plots
            **kwargs: various arguments to pass to the plotting function (see `plotly.express.scatter`)

        Returns:
            `None`
        """
        table = ManifoldTable.fromMolTable(table, name=f"{table.name}_manifold")
        manifold_cols = table.addManifoldData(self.manifold, recalculate=recalculate) if self.manifold else (x, y)
        if not manifold_cols[0] and not manifold_cols[1]:
            raise ValueError("Neither manifold nor x and y were specified.")

        kwargs['height'] = 800 if 'height' not in kwargs else kwargs['height']
        kwargs['width'] = 2*kwargs['height'] if 'width' not in kwargs else kwargs['width']
        kwargs['render_mode'] = 'svg' if 'render_mode' not in kwargs else kwargs['render_mode']
        x = manifold_cols[0] if not x else x
        y = manifold_cols[1] if not y else y
        if not color_by and table.hasScaffolds:
            scaffold = table.getScaffoldNames()[0] # FIXME: we should expose this and give a choice of what scaffold to use
            table.createScaffoldGroups(mols_per_group=mols_per_scaffold_group)
            color_by = table.getScaffoldGroups(f"{scaffold}", mols_per_scaffold_group).name
            color_discrete_map = {'Other': 'lightgrey'}
            df = table.getDF()
            fig = px.scatter(df, x=x, y=y,
                color = color_by, symbol=color_by,
                symbol_sequence = self.symbols,
                color_discrete_map = color_discrete_map,
                **kwargs
            )
        elif color_by:
            df = table.getDF()
            fig = px.scatter(df, x=x, y=y,
                color=color_by,
                **kwargs
            )
        else:
            df = table.getDF()
            fig = px.scatter(df, x=x, y=y,
                **kwargs
            )
        fig.update_layout(plot_bgcolor='White')

        if not interactive:
            return fig

        # interactive plot:
        excluded = df.columns[df.columns.str.contains('RDMol')].tolist() + list(table.getDescriptorNames()) + manifold_cols + df.columns[~df.columns.isin(card_data)].tolist()
        included = [title_data] + [col for col in df.columns if col not in excluded]
        smiles_col = [table.smilesCol] + table.getScaffoldNames() if table.hasScaffolds else [table.smilesCol]
        app_scatter = molplotly.add_molecules(fig=fig,
          df=df,
          smiles_col=smiles_col,
          title_col= title_data,
          color_col = color_by,
          caption_cols = included,
          # width = kwargs['width']
        )

        self.open_apps[port] = app_scatter
        app_scatter.run_server(
            mode='inline',
            port=port,
            height=viewport_height,
        )

class ModelPerformancePlot(ModelPlot):

    def __init__(self, manifold : Manifold, models: List[QSPRModel], ports, datasets : List[QSPRDataset] = None, card_props = None, plot_type = Literal["errors", "splits", "predictions", "labels"], async_execution=True):
        super().__init__(models)
        self.manifold = manifold
        self.plotType = plot_type
        self.ports = ports
        self.runningApps = dict()
        self.asyncExecution = async_execution
        self.cardProps = card_props if card_props else []

        # check if we have data for all models
        self.datasets = datasets if datasets else dict()
        if not self.datasets:
            for model in self.models:
                self.datasets[model] = model.data

        # check if ports unique
        if len(self.ports) != len(set(self.ports)):
            raise ValueError("Ports must be unique.")
        
        if len(self.datasets) != len(self.models):
            raise ValueError("Number of models and datasets does not match.")

        if len(self.ports) != len(self.models):
            raise ValueError("Number of models and ports does not match.")

        if not all(self.datasets):
            raise ValueError("Some models have no associated data. Specify the data used to train each model with the 'datasets' argument or provide a model with a 'QSPRDataset' attached.")

    def getSupportedTasks(self):
        """Return a list of tasks supported by this plotter."""
        return [
            ModelTasks.SINGLECLASS, 
            ModelTasks.MULTICLASS, 
            ModelTasks.REGRESSION,
        ]

    def getPerfCols(self, model, target_prop):
        """
        Get the relevant performance columns for a given model and target property.

        Args:
            model: `QSPRModel`
            target_prop: `TargetProperty`

        Returns:
            col_label: column name for the original label/target
            col_pred: column name for the prediction
            cols_probas: column names for the class probabilities if the target property is a classification task, empty list otherwise
        """

        col_label = f"{target_prop.name}_Label"
        if model.task.isClassification():
            col_pred = f"{target_prop.name}_Prediction"
            cols_probas = []
            for i in range(target_prop.nClasses):
                cols_probas.append(f"{target_prop.name}_ProbabilityClass_{i}")
        elif model.task.isRegression():
            col_pred = f"{target_prop.name}_Prediction"
            cols_probas = []
        else:
            raise NotImplementedError(f"Unsupported task: {model.task}")

        return col_label, col_pred, cols_probas


    def getPerfData(self, path, model, target_prop):
        df = pd.read_table(path, index_col=0)
        col_label, col_pred, cols_probas = self.getPerfCols(model, target_prop)
        col_err = f"{target_prop}_Error"
        df[col_err] = df[col_label] - df[col_pred]
        if model.task.isClassification():
            # convert True/False to string labels
            df[col_label] = [f"Class_{int(x)}" for x in df[col_label]]
            df[col_pred] = [f"Class_{int(x)}" for x in df[col_pred]]
        return df, col_label, col_pred, col_err, cols_probas

    def getCVData(self, model, target_prop):
        cv_path = self.cvPaths[model]
        df, col_label, col_pred, col_err, cols_probas = self.getPerfData(cv_path, model, target_prop)
        df["TestSet"] = [f"Fold_{int(x) + 1}" for x in df["Fold"]]
        del df["Fold"]
        return df, col_label, col_pred, col_err, cols_probas

    def getIndData(self, model, target_prop):
        ind_path = self.indPaths[model]
        df, col_label, col_pred, col_err, cols_probas = self.getPerfData(ind_path, model, target_prop)
        df["TestSet"] = "Independent"
        return df, col_label, col_pred, col_err, cols_probas

    def make(self, show=True, save=False):
        """Make the plot."""

        for model_idx, model in enumerate(self.models):
            port = self.ports[model_idx]
            ds = self.datasets[model]
            df_cv, col_label, col_pred, col_err, cols_probas = self.getCVData(model, model.targetProperties[0])
            df_ind, col_label, col_pred, col_err, cols_probas = self.getIndData(model, model.targetProperties[0])
            df_all = pd.concat([df_cv, df_ind])

            # create a molecule table with the required data
            manifold_cols = ds.getSubset(f"{self.manifold}_")
            if manifold_cols is None:
                manifold_cols = []
            else:
                manifold_cols = manifold_cols.columns.tolist()
            ds_subset = ds.getDF()[[ds.smilesCol] + self.cardProps + ds.indexCols + manifold_cols]
            df_all = ds_subset.merge(df_all, left_index=True, right_index=True)
            mt = MoleculeTable(f"{model.name}_perfplot_{self.plotType}_p{port}", df=df_all, smilesCol=ds.smilesCol, index_cols=ds.indexCols)
            features = ds.getFeatures(concat=True)
            calc = CustomDescriptorsCalculator([DataFrameDescriptorSet(features)])
            mt.addCustomDescriptors(calc)

            # create a server for the current plot
            plt_map = {
                "errors" : col_err,
                "splits" : "TestSet",
                "predictions" : col_pred,
                "labels" : col_label,
            }
            def plot_server(port, runningApps):
                def run():
                    plot = Plot(manifold=self.manifold)
                    runningApps[port]["plot"] = plot
                    plot.plot(
                        mt,
                        title_data=mt.indexCols[0],
                        card_data=mt.indexCols + ["TestSet", col_label, col_pred, col_err] + cols_probas + self.cardProps,
                        color_by=plt_map[self.plotType],
                        port=port,
                        interactive=True,
                        recalculate=False,
                    )

                return run

            server = plot_server(port, self.runningApps)
            self.runningApps[port] = dict()
            self.runningApps[port]["model"] = model
            self.runningApps[port]['server'] = server
            self.runningApps[port]['plot_type'] = self.plotType
            self.runningApps[port]['table'] = mt
            if self.asyncExecution:
                thr = threading.Thread(target=plot_server(port, self.runningApps))
                thr.start()
                # sleep to give the server time to start
                time.sleep(1)
                self.runningApps[port]['thread'] = thr
            else:
                server()

        return self.runningApps