"""
scaffold

Created by: Martin Sicho
On: 05.10.22, 10:07
"""
from abc import ABC, abstractmethod

import pandas as pd
from rdkit.Chem import PandasTools

class Scaffold(ABC):
    """
    Abstract base class for calculating scaffold.
    """

    @abstractmethod
    def __call__(self, mol):
        """
        Calculate the scaffold for a molecule.

        Args:
            mol: smiles or rdkit molecule

        Returns:
            smiles of the scaffold
        """
        pass

    @abstractmethod
    def __str__(self):
        pass

class Murcko(Scaffold):

    def __call__(self, mol):
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
        mol = Chem.MolFromSmiles(mol) if isinstance(mol, str) else mol
        return MurckoScaffold.MurckoScaffoldSmilesFromSmiles(Chem.MolToSmiles(mol))

    def __str__(self):
        return "Murcko"