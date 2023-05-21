"""
example_coloring

Created by: Martin Sicho
On: 02.11.22, 14:09
"""
from qsprpred.data.utils.descriptorcalculator import MoleculeDescriptorsCalculator
from qsprpred.data.utils.descriptorsets import FingerprintSet
from scaffviz.clustering.manifold import TSNE
from qsprpred.data.data import MoleculeTable
from scaffviz.depiction.plot import Plot

if __name__ == '__main__':
    # load data
    dataset = MoleculeTable.fromTableFile("P51681", "./data/P51681_LIGANDS_nostereo.tsv", sep="\t", store_dir="./data")
    desc_calculator = MoleculeDescriptorsCalculator(descsets=[FingerprintSet(fingerprint_type="MorganFP", radius=3, nBits=2048)])
    dataset.addDescriptors(desc_calculator, recalculate=False)

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
