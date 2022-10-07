"""
dataset

Created by: Martin Sicho
On: 05.10.22, 11:47
"""
import os
from abc import ABC, abstractmethod
from typing import List

import numpy as np
import pandas as pd
from pandas import DataFrame
from rdkit.Chem import PandasTools

from src.clustering.descriptors import Descriptor
from src.clustering.scaffolds import Scaffold


class DataSet(ABC):

    @abstractmethod
    def asDataFrame(self, smiles_col = 'SMILES', mol_col = 'RDMol') -> DataFrame:
        pass

    @abstractmethod
    def addDescriptors(self, descriptors : List[Descriptor]):
        """
        Add descriptors to the dataset.

        Args:
            descriptors: a list of descriptors to calculate and add to the dataset

        Returns:
            `None`
        """
        pass

    @abstractmethod
    def addScaffolds(self, scaffolds : List[Scaffold]):
        """
        Add scaffolds to the dataset.

        Args:
            scaffolds: scaffolds to calculate and add to the dataset

        Returns:
            `None`
        """
        pass

    @abstractmethod
    def getDescriptors(self):
        """
        Get the descriptors that are currently in the dataset.

        Returns:
            a `DataFrame` with the descriptors
        """

        pass

    @abstractmethod
    def getScaffolds(self):
        """
        Get the scaffolds that are currently in the dataset.

        Returns:
            a `DataFrame` with the scaffolds
        """
        pass

    @abstractmethod
    def getDescriptorsNames(self):
        """
        Get the names of the descriptors that are currently in the dataset.

        Returns:
            a `list` of descriptor names
        """
        pass

    @abstractmethod
    def getScaffoldNames(self):
        """
        Get the names of the scaffolds that are currently in the dataset.

        Returns:
            a `list` of scaffold names
        """
        pass

    @abstractmethod
    def getMols(self, asSmiles = True):
        """
        Get the molecules that are currently in the dataset.

        Args:
            asSmiles: if `True` return the molecules as SMILES strings, otherwise return RDKit molecules

        Returns:
            a `list` of rdkit molecules or SMILES strings
        """
        pass

    @property
    @abstractmethod
    def hasDescriptors(self):
        pass

    @property
    @abstractmethod
    def hasScaffolds(self):
        pass

    @property
    @abstractmethod
    def hasScaffoldGroups(self):
        pass

    @abstractmethod
    def getSubset(self, prefix : str):
        pass

    @abstractmethod
    def getNames(self):
        pass

    @abstractmethod
    def addData(self, name, data):
        pass

    @abstractmethod
    def removeData(self, name):
        pass

    @abstractmethod
    def createScaffoldGroups(self, mols_per_group : int = 10):
        pass

    @abstractmethod
    def getScaffoldGroups(self, scaffold_name : str):
        pass


class DataSetTSV(DataSet):

    def __init__(self, path, mol_col = 'SMILES', data : DataFrame = None):
        self.smilesCol = mol_col
        self.molCol = "RDMol"
        self.path = path
        if data is not None:
            self._data = data
        else:
            self._data = pd.read_csv(self.path, sep='\t') if os.path.exists(self.path) else None
        PandasTools.AddMoleculeColumnToFrame(self._data, smilesCol=mol_col, molCol=self.molCol)

    def asDataFrame(self, smiles_col = 'SMILES'):
        if smiles_col == self.smilesCol:
            return self._data
        else:
            return self._data.rename(columns={self.smilesCol: smiles_col})

    def addDescriptors(self, descriptors: List[Descriptor], recalculate = False):
        for descriptor in descriptors:
            recalculate = recalculate or len([x for x in self.getDescriptorsNames() if x.startswith(f"Descriptor_{descriptor}")]) == 0
            if not recalculate:
                continue
            values = self._data.apply(lambda row: descriptor(row[self.smilesCol]), axis=1).to_list()
            values = pd.DataFrame(values, columns=[f"Descriptor_{descriptor}_{idx}" for idx in range(len(values[0]))])
            self._data = pd.concat([self._data, values], axis=1)
            self.save()

    def addScaffolds(self, scaffolds: List[Scaffold]):
        for scaffold in scaffolds:
            if f"Scaffold_{scaffold}" in self._data.columns:
                continue

            self._data[f"Scaffold_{scaffold}"] = self._data.apply(lambda row: scaffold(row[self.smilesCol]), axis=1)
            PandasTools.AddMoleculeColumnToFrame(self._data, smilesCol=f"Scaffold_{scaffold}", molCol=f"Scaffold_{scaffold}_{self.molCol}")
            self.save()

    def getDescriptorsNames(self):
        return [col for col in self._data.columns if col.startswith("Descriptor_")]

    def getScaffoldNames(self):
        return [col for col in self._data.columns if col.startswith("Scaffold_")]

    def getDescriptors(self):
        return self._data[self.getDescriptorsNames()]

    def getScaffolds(self, includeMols = False):
        if includeMols:
            return self._data[[col for col in self._data.columns if col.startswith("Scaffold_")]]
        else:
            return self._data[self.getScaffoldNames()]

    def getMols(self, asSmiles=True):
        return self._data[self.smilesCol] if asSmiles else self._data[self.molCol]

    @property
    def hasDescriptors(self):
        return len(self.getDescriptorsNames()) > 0

    @property
    def hasScaffolds(self):
        return len(self.getScaffoldNames()) > 0

    @property
    def hasScaffoldGroups(self):
        return len([col for col in self._data.columns if col.startswith("ScaffoldGroup_")]) > 0

    def addData(self, name, data):
        self._data[name] = data
        self.save()

    def removeData(self, name):
        self._data.drop(name, axis=1, inplace=True)
        self.save()

    def createScaffoldGroups(self, mols_per_group = 10):
        scaffolds = self.getScaffolds(includeMols=False)
        for scaffold in scaffolds.columns:
            counts = pd.value_counts(self._data[scaffold])
            mask = counts.lt(mols_per_group)
            self._data[f'ScaffoldGroup_{scaffold}'] = np.where(self._data[scaffold].isin(counts[mask].index),'Other', self._data[scaffold])
        self.save()

    def getScaffoldGroups(self, scaffold_name : str):
        return self._data[f'ScaffoldGroup_{scaffold_name}']

    def save(self, path = None):
        self._data.to_csv(path if path else self.path, sep='\t', index=False, header=True)

    def getSubset(self, prefix : str):
        if self._data.columns.str.startswith(prefix).any():
            return self._data[self._data.columns[self._data.columns.str.startswith(prefix)]]

    def getNames(self):
        return self._data.columns

