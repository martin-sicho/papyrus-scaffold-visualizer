"""
plot

Created by: Martin Sicho
On: 05.10.22, 16:37
"""
import molplotly
import plotly.express as px

from scaffviz.clustering.manifold import Manifold
from scaffviz.data.dataset import DataSet


class Plot:

    def __init__(self, dataset : DataSet):
        self.dataset = dataset
        self.symbols = ['circle', 'square', 'diamond', 'cross', 'x',  'pentagon', 'hexagram', 'star', 'diamond', 'hourglass', 'bowtie']
        self.open_apps = dict()

    def getOpenApps(self):
        return self.open_apps

    def plot(self, manifold : Manifold = None, x : str = None, y : str = None, color_by : str = None, card_data = tuple(), title_data : str = 'SMILES', port=9292, recalculate=False, mols_per_scaffold_group : int = 10, interactive = True, viewport_height = "100%",  **kwargs):
        """
        Plot the dataset using the manifold or custom `DataSet` fields. The plot is interactive and runs as a web app on the specified port.

        Args:
            manifold: the manifold to use for embedding the data into the plot, if not specified the 'x' and 'y' parameters must be set
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
        manifold_cols = self.dataset.addManifoldData(manifold, recalculate=recalculate) if manifold else (x, y)
        if not manifold_cols[0] and not manifold_cols[1]:
            raise ValueError("Neither manifold nor x and y were specified.")

        kwargs['height'] = 800 if 'height' not in kwargs else kwargs['height']
        kwargs['width'] = 2*kwargs['height'] if 'width' not in kwargs else kwargs['width']
        kwargs['render_mode'] = 'svg' if 'render_mode' not in kwargs else kwargs['render_mode']
        x = manifold_cols[0] if not x else x
        y = manifold_cols[1] if not y else y
        if not color_by and self.dataset.hasScaffolds:
            scaffold = self.dataset.getScaffoldNames()[0]
            self.dataset.createScaffoldGroups(mols_per_group=mols_per_scaffold_group)
            color_by = self.dataset.getScaffoldGroups(f"{scaffold}", mols_per_scaffold_group).name
            color_discrete_map = {'Other': 'lightgrey'}
            df = self.dataset.asDataFrame(smiles_col='SMILES', mol_col='RDMol')
            fig = px.scatter(df, x=x, y=y,
                color = color_by, symbol=color_by,
                symbol_sequence = self.symbols,
                color_discrete_map = color_discrete_map,
                **kwargs
            )
        elif color_by:
            df = self.dataset.asDataFrame(smiles_col='SMILES', mol_col='RDMol')
            fig = px.scatter(df, x=x, y=y,
                color=color_by,
                **kwargs
            )
        else:
            df = self.dataset.asDataFrame(smiles_col='SMILES', mol_col='RDMol')
            fig = px.scatter(df, x=x, y=y,
                **kwargs
            )
        fig.update_layout(plot_bgcolor='White')

        if not interactive:
            return fig

        # interactive plot:
        excluded = df.columns[df.columns.str.contains('RDMol')].tolist() + list(self.dataset.getDescriptorsNames()) + manifold_cols + df.columns[~df.columns.isin(card_data)].tolist()
        included = ['SMILES'] + [col for col in df.columns if col not in excluded]
        app_scatter = molplotly.add_molecules(fig=fig,
          df=df,
          smiles_col=['SMILES'] + self.dataset.getScaffoldNames(),
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
