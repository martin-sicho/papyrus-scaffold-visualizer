# Papyrus Scaffold Visualizer

This package is for everyone who wants to quickly explore publicly available data before conducting further analysis or any kind of machine learning. Extraction and interactive visualization of any public data set with the help of [Papyrus](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-022-00672-x). This package allows you to quickly filter out relevant data for your target of interest and will make an interactive visualization of the compounds in the data set. By default, the compounds are grouped and colored by scaffold, but you can also customize the grouping and coloring by other means. Thanks to the nature of Papyrus, the assembled data set can also be used for subsequent machine learning tasks with very little effort.

## Getting started

You can install the package from this repository using pip:

```bash
   pip install git+https://github.com/martin-sicho/papyrus-scaffold-visualizer.git@main
```

Here is an example script that extracts a data set and makes the visualization ([examples/example_depiction_papyrus.py](./examples/example_depiction_papyrus.py)):

```python
from qsprpred.data.chem.scaffolds import BemisMurcko
from qsprpred.data.descriptors.fingerprints import MorganFP
from scaffviz.clustering.manifold import TSNE
from qsprpred.data.sources.papyrus import Papyrus
from scaffviz.depiction.plot import Plot

acc_keys = ["P51681"]  # accession keys of the proteins to fetch data for
name = "P51681_LIGANDS_nostereo"  # name of the data set
quality = "low"  # choose minimum quality from {"high", "medium", "low"}

# fetch data from Papyrus for the given targets (may take a while)
papyrus = Papyrus(
    data_dir="./data",  # where to store the data
    stereo=False,  # do not consider stereochemistry
    plus_only=True,  # only use the high quality data set
    descriptors=None,  # do not download descriptors (we will calculate them later)
)
dataset = papyrus.getData(
    name,
    acc_keys,
    quality,
    use_existing=True # use existing data set if it was already compiled before
)

# add generic scaffolds to the data set
# these will be used to group the molecules
dataset.addScaffolds([BemisMurcko()])

# add Morgan fingerprints to the data set
# these will be used to calculate the t-SNE embedding in 2D
dataset.addDescriptors([MorganFP(radius=3, nBits=2048)], recalculate=False)

# make an interactive plot that will use t-SNE to embed the data set in 2D
# (all available descriptors in the data set will be used, not just selected features)
plt = Plot(TSNE(perplexity=150))
plt.plot(
    dataset,
    recalculate=False,  # do not recalculate the t-SNE embedding each time this is run
    mols_per_scaffold_group=5,  # smaller groups will be merged into 'Other' group
    card_data=["all_doc_ids"],  # what to show on the molecule cards
    title_data='InChIKey'   # Data to show in the title of the molecule cards
)
```

You can find more example scripts under [examples](./examples).

## License
[MIT License](./LICENSE.md).

## Credits

This work is based on the code developed by Nadieh van Marrewijk and Yorick van Aalst during their internships at the Computational Drug Discovery Group, LACDR, Leiden University. We would also like to thank the authors of the [molplotly](https://github.com/wjm41/molplotly) package for its continued development and support.
