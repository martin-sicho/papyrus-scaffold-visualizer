"""
example_depiction_smiles

Created by: Martin Sicho
On: 10.10.22, 16:05
"""
from qsprpred.data.utils.descriptorcalculator import DescriptorsCalculator
from qsprpred.data.utils.descriptorsets import MorganFP
from qsprpred.data.data import MoleculeTable
from qsprpred.data.utils.scaffolds import BemisMurcko

from scaffviz.clustering.manifold import TSNE
from scaffviz.depiction.plot import Plot

if __name__ == "__main__":

    dataset = MoleculeTable.fromSMILES('smiles', ['CN1C2CCC1C(C(C2)OC(=O)C3=CC=CC=C3)C(=O)OC', 'O=C(OCCN(CC)CC)c1ccc(N)cc1', 'CCO'], store_dir='data')
    dataset.addProperty("Name", ["cocaine", "procaine", "ethanol"])

    desc_calculator = DescriptorsCalculator(descsets=[MorganFP(radius=2, nBits=2048)])
    dataset.addDescriptors(desc_calculator, recalculate=True)
    dataset.addScaffolds([BemisMurcko(convert_hetero=False)])

    plt = Plot(TSNE(perplexity=1))
    plt.plot(dataset, recalculate=True, mols_per_scaffold_group=1, card_data=['Name'], title_data='Name', port=9292)
