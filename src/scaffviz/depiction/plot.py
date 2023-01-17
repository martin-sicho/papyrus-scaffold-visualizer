"""
plot

Created by: Martin Sicho
On: 05.10.22, 16:37
"""
import molplotly
import plotly.express as px

from scaffviz.clustering.manifold import Manifold
from qsprpred.data.data import MoleculeTable
from scaffviz.data.manifold_table import ManifoldTable


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
        table = ManifoldTable.fromMolTable(table)
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
        included = ['SMILES'] + [col for col in df.columns if col not in excluded]
        smiles_col = ['SMILES'] + table.getScaffoldNames() if table.hasScaffolds else ['SMILES']
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
