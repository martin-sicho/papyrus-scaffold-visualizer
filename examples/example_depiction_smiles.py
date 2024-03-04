"""
example_depiction_smiles

Created by: Martin Sicho
On: 10.10.22, 16:05
"""
from qsprpred.data import MoleculeTable
from qsprpred.data.chem.scaffolds import Murcko
from qsprpred.data.descriptors.fingerprints import MorganFP

from scaffviz.clustering.manifold import TSNE
from scaffviz.depiction.plot import Plot

if __name__ == "__main__":
    # create a dataset
    dataset = MoleculeTable.fromSMILES('smiles', ['CN1C2CCC1C(C(C2)OC(=O)C3=CC=CC=C3)C(=O)OC', 'O=C(OCCN(CC)CC)c1ccc(N)cc1', 'CCO'])
    dataset.addProperty("Name", ["cocaine", "procaine", "ethanol"])
    # add descriptors and scaffolds
    dataset.addDescriptors([MorganFP(radius=3, nBits=2048)], recalculate=True)
    dataset.addScaffolds([Murcko()])
    # plot the dataset
    plt = Plot(TSNE(perplexity=1))
    plt.plot(dataset, recalculate=True, mols_per_scaffold_group=1, title_data="Name",card_data=["Name"], port=9292)
