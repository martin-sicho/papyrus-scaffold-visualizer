"""
descriptors

Created by: Martin Sicho
On: 05.10.22, 10:13
"""
from abc import abstractmethod, ABC


class Descriptor(ABC):
    """
    Abstract base class for descriptors.

    """

    @abstractmethod
    def __call__(self, mol):
        """
        Calculate the descriptor for a molecule.

        Args:
            mol: smiles or rdkit molecule
            *args: optional arguments
            **kwargs: optional keyword arguments

        Returns:
            a `list` of descriptor values
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

    def __call__(self, mol):
        mol = self._convertMol(mol) if isinstance(mol, str) else mol
        return self._convertFP(self._morgan(mol, *self._args, **self._kwargs))

    def __str__(self):
        return "MorganFP"
