"""
Show depiction from a custom data set.

Created by: Martin Sicho
On: 10.10.22, 12:01
"""
from qsprpred.data.data import MoleculeTable
from qsprpred.data.utils.descriptorcalculator import MoleculeDescriptorsCalculator
from qsprpred.data.utils.descriptorsets import FingerprintSet
from scaffviz.clustering.manifold import TSNE
from qsprpred.data.utils.scaffolds import Murcko
from scaffviz.depiction.plot import Plot

if __name__ == "__main__":

    dataset = MoleculeTable.fromSDF("A2A", "./data/A2A_chembl_export.sdf", smiles_prop="CANONICAL_SMILES", store_dir="./data")

    dataset.addScaffolds([Murcko()])
    desc_calculator = MoleculeDescriptorsCalculator(descsets=[FingerprintSet(fingerprint_type="MorganFP", radius=3, nBits=2048)])
    dataset.addDescriptors(desc_calculator, recalculate=False)

    plt = Plot(TSNE())
    plt.plot(dataset, recalculate=False, mols_per_scaffold_group=10, card_data=['INCHI_KEY', 'CANONICAL_SMILES'], title_data='CANONICAL_SMILES', port=9292)

