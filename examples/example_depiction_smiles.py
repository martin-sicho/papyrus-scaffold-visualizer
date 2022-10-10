"""
example_depiction_smiles

Created by: Martin Sicho
On: 10.10.22, 16:05
"""
from clustering.descriptors import MorganFP
from clustering.manifold import TSNE
from clustering.scaffolds import Murcko
from data.dataset import DataSetSMILES
from depiction.plot import Plot

if __name__ == "__main__":

    dataset = DataSetSMILES('./data/smiles.tsv', ['CN1C2CCC1C(C(C2)OC(=O)C3=CC=CC=C3)C(=O)OC', 'O=C(OCCN(CC)CC)c1ccc(N)cc1', 'CCO'])
    dataset.addData("Name", ["cocaine", "procaine", "ethanol"])

    dataset.addDescriptors([MorganFP(radius=2, nBits=1024)], recalculate=True)

    plt = Plot(dataset, TSNE(perplexity=1))
    plt.plot(recalculate=True, mols_per_scaffold_group=10, card_data=['Name'], title_data='Name', color_by='Name')
