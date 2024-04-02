"""
utils

Created by: Martin Sicho
On: 22.05.23, 10:12
"""
import logging
import os

from qsprpred.data import QSPRDataset, MoleculeTable
from qsprpred.data.sources.papyrus import Papyrus
from qsprpred.data import ScaffoldSplit
from qsprpred.data.descriptors.fingerprints import MorganFP
from qsprpred.data.processing.feature_filters import LowVarianceFilter, HighCorrelationFilter
from qsprpred.data.chem.scaffolds import BemisMurcko
from qsprpred.models.scikit_learn import SklearnModel


def fetch_example_dataset():
    """
    Use Papyrus to fetch the example set of ligands for the P51681 target.

    Returns:
        MoleculeTable: the example data set
    """
    # try to load from file
    name = "P51681_LIGANDS_nostereo"
    try:
        return MoleculeTable.fromFile(f"./data/{name}/{name}_meta.json")
    except FileNotFoundError:
        logging.warning("Data set not found. Fetching from Papyrus.")
    # otherwise fetch from Papyrus
    acc_keys = ["P51681"]
    quality = "high"
    papyrus = Papyrus(
        data_dir="./data",
        stereo=False,
        descriptors=None,
        plus_only=True
    )
    return papyrus.getData(
        name,
        acc_keys,
        quality,
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
    dataset.featureNames = mol_table.getDescriptorNames()
    dataset.loadDescriptorsToSplits()
    if not force_build and not dataset.hasDescriptors():
        split = ScaffoldSplit(
            scaffold=BemisMurcko(),
            test_fraction=0.2
        )
        lv = LowVarianceFilter(0.05)
        hc = HighCorrelationFilter(0.8)
        dataset.prepareDataset(
            split=split,
            feature_calculators=[MorganFP(radius=3, nBits=2048)],
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
    used_datasets = []
    for model, prop in zip(models, target_props):
        dataset = prepare_example_dataset(dataset, prop, force_build=force_build)
        model = SklearnModel(
            base_dir='data',
            alg=model,
            name=model.__name__
        )

        # only train if required
        if force_build or not os.path.exists(model.metaFile):
            from qsprpred.models.assessment.methods import CrossValAssessor, \
                TestSetAssessor
            metric = "roc_auc" if dataset.targetProperties[0].task.isClassification() else "r2"
            CrossValAssessor(scoring=metric)(model, dataset)
            TestSetAssessor(scoring=metric)(model, dataset)
            model.fitDataset(dataset)

        fitted_models.append(model)
        used_datasets.append(dataset)

    return fitted_models, used_datasets
