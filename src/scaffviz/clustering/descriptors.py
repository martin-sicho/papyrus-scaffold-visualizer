"""
descriptors

Created by: Martin Sicho
On: 05.10.22, 10:13
"""
from abc import abstractmethod, ABC

import numpy as np


class Descriptor(ABC):
    """
    Abstract base class for descriptors.

    """

    @abstractmethod
    def __call__(self, mol):
        """
        Calculate descriptors for a molecule.

        Args:
            mol: smiles or rdkit molecule

        Returns:
            a numpy array of descriptor values
        """
        pass

    @abstractmethod
    def __str__(self):
        pass

class MorganFP(Descriptor):

    def __init__(self, *args, **kwargs):
        """
        Initialize the descriptor with the same arguments as you would pass to `GetMorganFingerprintAsBitVect` function of RDKit.

        Args:
            *args: `GetMorganFingerprintAsBitVect` arguments
            **kwargs: `GetMorganFingerprintAsBitVect` keyword arguments
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
        self._convertMol = Chem.MolFromSmiles
        self._convertFP = DataStructs.ConvertToNumpyArray
        self._morgan = AllChem.GetMorganFingerprintAsBitVect
        self._args = args
        self._kwargs = kwargs
        self._ln = None

    def __call__(self, mol):
        mol = self._convertMol(mol) if isinstance(mol, str) else mol
        ret = np.zeros(self._ln)
        if not mol:
            return ret
        fp = self._morgan(mol, *self._args, **self._kwargs)
        if not self._ln:
            self._ln = len(fp)
        self._convertFP(fp, ret)
        return ret

    def __str__(self):
        return "MorganFP"
