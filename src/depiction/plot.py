"""
plot

Created by: Martin Sicho
On: 05.10.22, 16:37
"""
import molplotly
import numpy as np
import pandas as pd
import plotly.express as px

from src.clustering.manifold import Manifold
from src.data.dataset import DataSet


class Plot:

    def __init__(self, dataset : DataSet, manifold : Manifold):
        self.dataset = dataset
        self.manifold = manifold
        self.symbols = ['circle', 'square', 'diamond', 'cross', 'x',  'pentagon', 'hexagram', 'star', 'diamond', 'hourglass', 'bowtie']

    def plot(self, color_by : str = None, card_data = tuple(), port=9292, recalculate=True, mols_per_scaffold_group : int = 10,  **kwargs):
        """
        Plot the dataset using the manifold. The plot is interactive and runs as a web app on the specified port.

        Args:
            color_by: the data to color the points by, by default the first scaffold found will be used
            card_data: list of data to show on the cards displayed when hovering over a molecule
            port: port to run the web app on
            recalculate: whether to recalculate the manifold or use the existing data in the dataset
            mols_per_scaffold_group: how many molecules to include in one scaffold group
            **kwargs: various arguments to pass to the plotting function (see `plotly.express.scatter`)

        Returns:
            `None`
        """

        manifold_data = self.dataset.getSubset(str(self.manifold))
        manifold_cols = []
        if manifold_data is not None:
            manifold_cols = manifold_data.columns.tolist()
        if recalculate or manifold_data is None:
            X = pd.DataFrame(self.manifold.fit_transform(self.dataset))
            manifold_cols = []
            for i, dim in enumerate(X.columns):
                col_name = f"{self.manifold}_{i+1}"
                manifold_cols.append(col_name)
                self.dataset.addData(col_name, X[dim])

        kwargs['height'] = 800 if 'height' not in kwargs else kwargs['height']
        kwargs['width'] = 2*kwargs['height'] if 'width' not in kwargs else kwargs['width']
        kwargs['render_mode'] = 'svg' if 'render_mode' not in kwargs else kwargs['render_mode']
        x = f"{self.manifold}_1"
        y = f"{self.manifold}_2"
        scaffold = None
        if not color_by and self.dataset.hasScaffolds:
            scaffold = self.dataset.getScaffoldNames()[0]
            self.dataset.createScaffoldGroups(mols_per_group=mols_per_scaffold_group)
            color_by = self.dataset.getScaffoldGroups(scaffold).name
            color_discrete_map = {'Other': 'lightgrey'}
            df = self.dataset.asDataFrame(smiles_col='SMILES')
            fig = px.scatter(df, x=x, y=y,
                color = color_by, symbol=color_by,
                symbol_sequence = self.symbols,
                color_discrete_map = color_discrete_map,
                **kwargs
            )
        elif color_by:
            df = self.dataset.asDataFrame(smiles_col='SMILES')
            fig = px.scatter(df, x=x, y=y,
                color = color_by,symbol=color_by,
                symbol_sequence = self.symbols,
                **kwargs
            )
        else:
            df = self.dataset.asDataFrame(smiles_col='SMILES')
            fig = px.scatter(df, x=x, y=y,
                **kwargs
            )
        fig.update_layout(plot_bgcolor='White')

        # interactive plot:
        excluded = df.columns[df.columns.str.contains('RDMol')].tolist() + list(self.dataset.getDescriptorsNames()) + manifold_cols + df.columns[~df.columns.isin(card_data)].tolist()
        included = ['SMILES'] + [col for col in df.columns if col not in excluded]
        app_scatter = molplotly.add_molecules(fig=fig,
          df=df,
          smiles_col=['SMILES'] + self.dataset.getScaffoldNames(),
          title_col= 'SMILES',
          color_col = color_by,
          caption_cols = included,
          # width = kwargs['width']
        )

        app_scatter.run_server(
            mode='inline',
            port=port,
            height="100%",
        )
