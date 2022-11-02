"""
example_coloring

Created by: Martin Sicho
On: 02.11.22, 14:09
"""
from scaffviz.clustering.descriptors import MorganFP
from scaffviz.clustering.manifold import TSNE
from scaffviz.data.dataset import DataSetTSV
from scaffviz.depiction.plot import Plot

if __name__ == '__main__':
    dataset = DataSetTSV("./data/P51681_LIGANDS_nostereo.tsv")
    dataset.addDescriptors([MorganFP(radius=2, nBits=1024)], recalculate=False)

    plt = Plot(dataset)
    plt.plot(TSNE(), color_by='pchembl_value_Median', recalculate=False, mols_per_scaffold_group=5, port=9292)
