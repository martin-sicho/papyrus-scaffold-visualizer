"""
plot

Created by: Martin Sicho
On: 05.10.22, 16:37
"""
import molplotly
import numpy as np
import pandas as pd
import plotly.express as px

from clustering.manifold import Manifold
from data.dataset import DataSet


class Plot:

    def __init__(self, dataset : DataSet, manifold : Manifold):
        self.dataset = dataset
        self.manifold = manifold
        self.symbols = ['circle', 'square', 'diamond', 'cross', 'x',  'pentagon', 'hexagram', 'star', 'diamond', 'hourglass', 'bowtie']

    def plot(self, color_by : str = None, exclude_data = tuple(), port=9292,  **kwargs):
        """
        Plot the dataset using the manifold. The plot is interactive and runs as a web app on the specified port.

        Args:
            color_by: the data to color the points by, by default the first scaffold found will be used
            exclude_data: list of data to exclude from the plot
            port: port to run the web app on
            **kwargs: various arguments to pass to the plotting function (see `plotly.express.scatter`)

        Returns:
            `None`
        """

        X = pd.DataFrame(self.manifold.fit_transform(self.dataset.getDescriptors()))
        manifold_cols = []
        for i, dim in enumerate(np.transpose(X)):
            col_name = f"{self.manifold}_{i+1}"
            manifold_cols.append(col_name)
            self.dataset.addData(col_name, dim)

        kwargs['width'] = 800 if 'width' not in kwargs else kwargs['width']
        kwargs['height'] = 600 if 'height' not in kwargs else kwargs['height']
        kwargs['render_mode'] = 'svg' if 'render_mode' not in kwargs else kwargs['render_mode']
        df = self.dataset.asDataFrame(smiles_col='SMILES')
        x = f"{self.manifold}_1"
        y = f"{self.manifold}_2"
        scaffold = None
        if not color_by and self.dataset.hasScaffolds():
            scaffold = self.dataset.getScaffoldNames()[0]
            color_by = self.dataset.getScaffoldGroups(scaffold)
            color_discrete_map = {'Other': 'lightgrey'}
            fig = px.scatter(df, x=x, y=y,
                color = color_by, symbol=color_by,
                symbol_sequence = self.symbols,
                color_discrete_map = color_discrete_map,
                **kwargs
            )
        elif color_by:
            fig = px.scatter(df, x=x, y=y,
                color = color_by,symbol=color_by,
                symbol_sequence = self.symbols,
                **kwargs
            )
        else:
            fig = px.scatter(df, x=x, y=y,
                **kwargs
            )
        fig.update_layout(plot_bgcolor='White')

        # interactive plot:
        excluded = ['SMILES'] + list(self.dataset.getDescriptorsNames()) + manifold_cols + list(exclude_data)
        app_scatter = molplotly.add_molecules(fig=fig,
          df=df,
          smiles_col=['SMILES'] + self.dataset.getScaffoldNames(),
          title_col= 'SMILES',
          color_col = color_by,
          caption_cols = [col for col in df.columns if col not in excluded],
          width = kwargs['width']
        )

        app_scatter.run_server(mode='inline', port=port, height=kwargs['height'])
