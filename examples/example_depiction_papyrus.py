"""
example_depiction

Created by: Martin Sicho
On: 07.10.22, 11:16
"""
from scaffviz.clustering.descriptors import MorganFP
from scaffviz.clustering.manifold import TSNE
from scaffviz.clustering.scaffolds import Murcko
from scaffviz.data.dataset import DataSetTSV
from scaffviz.depiction.plot import Plot

if __name__ == "__main__":

    dataset = DataSetTSV("./data/hCCR2_LIGANDS_nostereo.tsv")
    dataset.addScaffolds([Murcko()])
    dataset.addDescriptors([MorganFP(radius=2, nBits=1024)], recalculate=False)

    plt = Plot(dataset, TSNE())
    plt.plot(recalculate=False, mols_per_scaffold_group=5)
