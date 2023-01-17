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
    def fromMolTable(mol_table : MoleculeTable):
        return ManifoldTable(mol_table.name, mol_table.getDF(), smilescol=mol_table.smilescol, store_dir=mol_table.storeDir)

    def getManifoldData(self, manifold: Manifold):
        return self.getSubset(str(manifold))

    def addManifoldData(self, manifold : Manifold, recalculate=True):
        manifold_data = self.getManifoldData(manifold)
        manifold_cols = []
        if manifold_data is not None:
            manifold_cols = manifold_data.columns.tolist()
        if recalculate or manifold_data is None:
            X = manifold.fit_transform(self.getDescriptors())
            manifold_cols = []
            x = np.transpose(X)
            for i, dim in enumerate(x):
                col_name = f"{manifold}_{i + 1}"
                manifold_cols.append(col_name)
                self.addProperty(col_name, dim)

        return manifold_cols
