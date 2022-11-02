"""
manifold

Created by: Martin Sicho
On: 05.10.22, 9:59
"""
from abc import abstractmethod, ABC


class Manifold(ABC):

    @abstractmethod
    def fit(self, X):
        """
        Fit the manifold to the data.

        Args:
            X: a matrix of data to fit the manifold to, can be a pandas `DataFrame` or a numpy `ndarray`

        Returns:
            `None`
        """

        pass

    @abstractmethod
    def transform(self, X):
        """
        Transform the data using the fitted manifold.

        Args:
            X: a matrix of data to transform, can be a pandas `DataFrame` or a numpy `ndarray`

        Returns:

        """

        pass

    def fit_transform(self, X):
        """
        Fit the manifold to the data and transform it in one step.

        Args:
            X: a matrix of data to fit the manifold to and transform, can be a pandas `DataFrame` or a numpy `ndarray`

        Returns:

        """

        self.fit(X)
        return self.transform(X)

    @abstractmethod
    def __str__(self):
        """
        Get a string representation of the defined manifold. Used to save the manifold in the `Dataset` table, for example.

        Returns:
            a string representation of the manifold
        """
        pass

class TSNE(Manifold):

    def __init__(self, *args, **kwargs):
        from sklearn.manifold import TSNE as skTSNE
        self._skTSNE = skTSNE(
            *args, **kwargs
        )

    def fit(self, X):
        self._skTSNE.fit(X)
        return self

    def transform(self, X):
        return self._skTSNE.fit_transform(X)

    def fit_transform(self, X):
        return self._skTSNE.fit_transform(X)

    def __str__(self):
        return "TSNE"
