"""
manifold

Created by: Martin Sicho
On: 05.10.22, 9:59
"""
from abc import abstractmethod, ABC

from scaffviz.data.dataset import DataSet


class Manifold(ABC):

    @abstractmethod
    def fit(self, X : DataSet, **kwargs):
        pass

    @abstractmethod
    def transform(self, X : DataSet, **kwargs):
        pass

    def fit_transform(self, X : DataSet, **kwargs):
        self.fit(X, **kwargs)
        return self.transform(X, **kwargs)

    @abstractmethod
    def __str__(self):
        pass

class TSNE(Manifold):

    def __init__(self, *args, **kwargs):
        from sklearn.manifold import TSNE as skTSNE
        self._skTSNE = skTSNE(
            *args, **kwargs
        )

    def fit(self, X, **kwargs):
        self._skTSNE.fit(X.getDescriptors(), **kwargs)
        return self

    def transform(self, X, **kwargs):
        return self._skTSNE.fit_transform(X.getDescriptors(), **kwargs)

    def fit_transform(self, X, **kwargs):
        return self._skTSNE.fit_transform(X.getDescriptors(), **kwargs)

    def __str__(self):
        return "TSNE"
