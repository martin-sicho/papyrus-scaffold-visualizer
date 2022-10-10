"""
The example code from the readme file.

Created by: Martin Sicho
On: 10.10.22, 18:06
"""
from scaffviz.clustering.descriptors import MorganFP
from scaffviz.clustering.manifold import TSNE
from scaffviz.clustering.scaffolds import Murcko
from scaffviz.data.papyrus.papyrus_class import Papyrus
from scaffviz.depiction.plot import Plot

if __name__ == "__main__":

    acc_keys = ["P41597"] # replace with your own accession key
    name = "P41597_LIGANDS_nostereo" # replace with your own name for the output data set file
    quality = "low" # choose minimum quality from {"high", "medium", "low"}

    # fetches latest version of Papyrus if not already available and filters out the relevant data
    papyrus = Papyrus(data_dir="./data", stereo=False)
    dataset = papyrus.getData(
        acc_keys,
        quality,
        name,
        use_existing=True # use existing data set if it was already compiled before
    )

    # use the filtered data to create the visualization
    dataset.addScaffolds([Murcko()]) # add Murcko scaffolds to the data set -> will be used to group compounds
    dataset.addDescriptors([MorganFP(radius=2, nBits=1024)], recalculate=False) # add Morgan fingerprints to the data set -> will be used to calculate the t-SNE embedding

    plt = Plot(dataset, TSNE()) # make an interactive plot that will use t-SNE to embed the data set in 2D (all descriptors in the data set will be used)
    plt.plot(recalculate=False, mols_per_scaffold_group=10, card_data=["doc_id"]) # start the server, you can open the plot in the browser


