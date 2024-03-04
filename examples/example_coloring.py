"""
example_coloring

Created by: Martin Sicho
On: 02.11.22, 14:09
"""
from qsprpred.data.descriptors.fingerprints import MorganFP

from examples.utils import fetch_example_dataset
from scaffviz.clustering.manifold import TSNE
from scaffviz.depiction.plot import Plot

if __name__ == '__main__':
    # load data
    dataset = fetch_example_dataset()
    dataset.addDescriptors([MorganFP(3, 2048)], recalculate=False)
    # plot
    plt = Plot(TSNE())
    plt.plot(
        dataset,
        color_by='pchembl_value_Median',
        recalculate=False,
        card_data=[
            'pchembl_value_Median'
        ],
        mols_per_scaffold_group=5,
        port=9292
    )
