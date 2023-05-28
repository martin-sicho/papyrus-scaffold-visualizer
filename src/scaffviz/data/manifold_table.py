"""
manifold_table

Created by: Martin Sicho
On: 17.01.23, 17:20
"""
import numpy as np
from qsprpred.data.data import MoleculeTable
from scaffviz.clustering.manifold import Manifold


class ManifoldTable(MoleculeTable):

    @staticmethod
    def fromMolTable(mol_table : MoleculeTable, name=None):
        name = name if name is not None else mol_table.name
        mt = ManifoldTable(name, mol_table.getDF(), smilescol=mol_table.smilescol, store_dir=mol_table.storeDir)
        mt.descriptors = mol_table.descriptors
        mt.descriptorCalculators = mol_table.descriptorCalculators
        return mt

    def getManifoldData(self, manifold: Manifold):
        return self.getSubset(str(manifold))

    def addManifoldData(self, manifold : Manifold, recalculate=True):
        manifold_data = self.getManifoldData(manifold)
        manifold_cols = []
        if manifold_data is not None:
            manifold_cols = manifold_data.columns.tolist()
        if recalculate or manifold_data is None:
            if not self.hasDescriptors:
                raise ValueError("Descriptors must be calculated before adding manifold data.")
            X = manifold.fit_transform(self.getDescriptors())
            manifold_cols = []
            x = np.transpose(X)
            for i, dim in enumerate(x):
                col_name = f"{manifold}_{i + 1}"
                manifold_cols.append(col_name)
                self.addProperty(col_name, dim)

        return manifold_cols
