"""
utils

Created by: Martin Sicho
On: 22.05.23, 10:12
"""
import os

from qsprpred.data.data import QSPRDataset
from qsprpred.data.sources.papyrus import Papyrus
from qsprpred.data.utils.datasplitters import ScaffoldSplit
from qsprpred.data.utils.descriptorcalculator import MoleculeDescriptorsCalculator
from qsprpred.data.utils.descriptorsets import FingerprintSet
from qsprpred.data.utils.featurefilters import LowVarianceFilter, HighCorrelationFilter
from qsprpred.data.utils.scaffolds import Murcko
from qsprpred.models.models import QSPRsklearn


def fetch_example_dataset():
    """
    Use Papyrus to fetch the example set of ligands for the P51681 target.

    Returns:
        MoleculeTable: the example data set
    """

    acc_keys = ["P51681"]
    name = "P51681_LIGANDS_nostereo"
    quality = "low"
    papyrus = Papyrus(
        data_dir="./data",
        stereo=False,
        descriptors=None,
        plus_only=False
    )
    return papyrus.getData(
        acc_keys,
        quality,
        name=name,
        use_existing=True  # use existing data set if it was already compiled before
    )


def prepare_example_dataset(mol_table, target_props, force_build=False):
    """
    Converts a molecule table to a QSPRDataset and prepares it for training.

    Args:
        mol_table: molecule table to convert
        target_props: target properties to use as specified in the QSPRPred package
        force_build: if True, the dataset will be prepared even if it already exists

    Returns:
        QSPRDataset: the prepared data set

    """

    dataset = QSPRDataset.fromMolTable(mol_table, target_props=target_props)
    if not force_build and not dataset.hasDescriptors:
        feature_calculator = MoleculeDescriptorsCalculator(
            desc_sets=[
                FingerprintSet(fingerprint_type="MorganFP", radius=3, nBits=2048)])
        split = ScaffoldSplit(
            scaffold=Murcko(),
            test_fraction=0.2
        )
        lv = LowVarianceFilter(0.05)
        hc = HighCorrelationFilter(0.8)
        dataset.prepareDataset(
            split=split,
            feature_calculators=[feature_calculator],
            feature_filters=[lv, hc]
        )
        dataset.save()
    else:
        print("Data set already prepared. Preparation skipped.")
    print(f"Number of samples train set: {len(dataset.y)}")
    print(
        f"Number of samples test set: {len(dataset.y_ind)}, {len(dataset.y_ind) / len(dataset.df) * 100}%")

    return dataset


def fetch_example_models(models, target_props, force_build=False):
    """
    Use the example data set to build example models if they do not exist. Reload old models otherwise.

    Args:
        models: classes of scikit-learn models to use
        target_props: target properties to use as specified in the QSPRPred package
        force_build: if True, the models will be built even if they already exist

    Returns:
        list of fitted and evaluated models
    """

    dataset = fetch_example_dataset()

    # train the models
    fitted_models = []
    for model, prop in zip(models, target_props):
        dataset = prepare_example_dataset(dataset, prop, force_build=force_build)
        model = QSPRsklearn(
            base_dir='data',
            data=dataset,
            alg=model,
            name=model.__name__
        )

        # only train if required
        if force_build or not os.path.exists(model.metaFile):
            model.evaluate()
            model.fit()

        fitted_models.append(model)

    return fitted_models
