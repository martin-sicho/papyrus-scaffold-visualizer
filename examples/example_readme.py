from qsprpred.data.utils.descriptorcalculator import MoleculeDescriptorsCalculator
from qsprpred.data.utils.descriptorsets import FingerprintSet
from scaffviz.clustering.manifold import TSNE
from qsprpred.data.utils.scaffolds import BemisMurcko
from qsprpred.data.sources.papyrus import Papyrus
from scaffviz.depiction.plot import Plot

acc_keys = ["P51681"] # replace with your own accession key
name = "P51681_LIGANDS_nostereo" # replace with your own name for the output data set file
quality = "low" # choose minimum quality from {"high", "medium", "low"}

# fetches the latest version of Papyrus if not already available and filters out the relevant data
papyrus = Papyrus(data_dir="./data", stereo=False)
dataset = papyrus.getData(
    acc_keys,
    quality,
    name=name,
    use_existing=True # use existing data set if it was already compiled before
)

# add Murcko scaffolds to the data set -> will be used to group compounds inside the plot
dataset.addScaffolds([BemisMurcko(convert_hetero=False, force_single_bonds=False)])

# add Morgan fingerprints to the data set -> these will be used to calculate the t-SNE embedding in 2D
desc_calculator = MoleculeDescriptorsCalculator(descsets=[FingerprintSet(fingerprint_type="MorganFP", radius=3, nBits=2048)])
dataset.addDescriptors(desc_calculator, recalculate=False)

# make an interactive plot that will use t-SNE to embed the data set in 2D (all available descriptors in the data set will be used)
plt = Plot(TSNE(perplexity=150))
plt.plot(dataset, recalculate=True, mols_per_scaffold_group=5, card_data=["all_doc_ids"], title_data='InChIKey')