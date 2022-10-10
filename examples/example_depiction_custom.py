"""
Show depiction from a custom data set.

Created by: Martin Sicho
On: 10.10.22, 12:01
"""
from scaffviz.clustering.descriptors import MorganFP
from scaffviz.clustering.manifold import TSNE
from scaffviz.clustering.scaffolds import Murcko
from scaffviz.data.dataset import DataSetSDF
from scaffviz.depiction.plot import Plot

if __name__ == "__main__":

    dataset = DataSetSDF("./data/A2A_chembl_export.sdf", smiles_prop="CANONICAL_SMILES")

    dataset.addScaffolds([Murcko()])
    dataset.addDescriptors([MorganFP(radius=2, nBits=1024)], recalculate=False)

    plt = Plot(dataset, TSNE())
    plt.plot(recalculate=False, mols_per_scaffold_group=10, card_data=['INCHI_KEY'], title_data='INCHI_KEY')