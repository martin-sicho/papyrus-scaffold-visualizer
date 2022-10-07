"""
example_depiction

Created by: Martin Sicho
On: 07.10.22, 11:16
"""
from clustering.descriptors import MorganFP
from clustering.manifold import TSNE
from clustering.scaffolds import Murcko
from data.dataset import DataSetTSV
from depiction.plot import Plot

if __name__ == "__main__":

    dataset = DataSetTSV("./data/hCCR2_LIGANDS_nostereo.tsv")
    dataset.addScaffolds([Murcko()])
    dataset.addDescriptors([MorganFP(radius=2, nBits=1024)], recalculate=False)

    plt = Plot(dataset, TSNE())
    plt.plot(recalculate=False, mols_per_scaffold_group=5)
