import numpy as np


def diagM(M: np.ndarray):
    """Diagonalize M

    Args:
        M (np.ndarray): Matrix to be diagonalized

    Returns:
        np.ndarray,np.ndarray: Eigenvalue vector and eigenvector matrix (columns are eigenvectors), sorted by increasing absolute eigenvalue
    """

    w, D = np.linalg.eig(M)

    sortIndices = np.argsort(np.abs(w))

    return w[sortIndices], D[:, sortIndices]


def greenlandMetric(M: np.ndarray) -> float:
    """Greenland metric of given array

    Args:
        M (np.ndarray): Array whose metric should be returned

    Returns:
        float: Greatest value of L1 metric amongst columns of M
    """

    return np.max(np.sum(np.abs(M), axis=1))


def decompT(T: np.ndarray, indices: np.ndarray):
    """Decomposes matrix T into 4 blocks by taking out the subspace determined by passed indices

    Args:
        T (np.ndarray): Array to decomposes
        indices (np.ndarray): Indices to decompose with

    Returns:
        np.ndarrays: 4 block matrices in following order: upper diagonal, lower diagonal, upper off-diagonal, lower off-diagonal
    """

    n = T.shape[0]

    mask = np.zeros(n, dtype=bool)
    mask[indices] = True

    compIndices = np.arange(n)[~mask]

    reshapedT = T[:, np.concatenate([indices, compIndices])]
    reshapedT = reshapedT[np.concatenate([indices, compIndices]), :]

    return reshapedT[0:len(indices), 0:len(indices)], \
        reshapedT[len(indices):, len(indices):], \
        reshapedT[0:len(indices), len(indices):], \
        reshapedT[len(indices):, 0:len(indices)]


class CRM:
    """CRM class containing effective rates and companion evolution matrices
    """

    def __init__(self, Meff: np.ndarray, dTqInv: np.ndarray, dTpInv: np.ndarray) -> None:
        """CRM constructor

        Args:
            Meff (np.ndarray): Effective rate matrix
            dTqInv (np.ndarray): $-\Delta T_Q^{-1}$ matrix in Greenland notation
            dTpInv (np.ndarray): $\delta T_P^{-1}$ matrix in Greenland notation
        """
        self.__Meff__ = Meff
        self.__dTqInv__ = dTqInv
        self.__dTpInv__ = dTpInv

    @property
    def Meff(self):
        return self.__Meff__

    @property
    def dTqInv(self):
        return self.__dTqInv__

    @property
    def dTpInv(self):
        return self.__dTpInv__


class CRMAnalyser:
    """Object responsible for generating CRMs from rate matrix
    """

    def __init__(self, M: np.ndarray) -> None:
        """CRMAnalyser constructor

        Args:
            M (np.ndarray): Rate matrix
        """

        self.__M__ = M
        self.__l__, self.__T__ = diagM(M)

    @property
    def eigenVals(self):
        return self.__l__

    @property
    def fullMat(self):
        return self.__M__

    @property
    def eigenVecs(self):
        return self.__T__

    def getPotentialCRMs(self, tol: float) -> 'list[np.ndarray]':
        """Return list of indices for which the eigenvector Greenland condition is satisfied

        Args:
            tol (float): Eigenvector condition tolerance

        Returns:
            list[np.ndarray]: List of indices for which eigenvector condition is satisfied
        """
        n = self.eigenVecs.shape[0]

        Z = np.zeros(n, dtype=int)

        for i in range(n):
            Z[i] = np.argmax(np.abs(self.eigenVecs[i, :]) >= tol)

        # TODO: better warning
        if any(Z == n):
            print("Warning: Singular indicator matrix")

        sortIndices = np.argsort(Z)
        firstPassIndices = []

        # Sort indices in new list according to eigenvalues in case of ties
        for i in range(1, n + 1):
            zIndices = np.where(Z[sortIndices] == i)
            sIndices = np.argsort(abs(self.eigenVals[sortIndices[zIndices]]))
            sortIndices[zIndices] = sortIndices[zIndices][sIndices]

        for i in range(1, n):
            if Z[sortIndices[i]] >= i:
                firstPassIndices.append(sortIndices[0:i])

        return firstPassIndices

    def getPQTimescales(self, indicesP: np.ndarray):
        """Returns P and Q timescales for CRM with given P indices

        Args:
            indicesP (np.ndarray): P space indices

        Returns:
            float,float: Shortest P timescale and longest Q timescales based on eigenvalues of M
        """

        N = len(self.eigenVals)
        mask = np.zeros(N, dtype=bool)
        mask[indicesP] = True

        tauP = np.abs(1 / np.max(self.eigenVals[mask]))
        tauQ = np.abs(1 / np.min(self.eigenVals[~mask]))

        return tauP, tauQ

    def generateCRM(self, indices: np.ndarray, exactMeff=False):
        """Generate a CRM based on P space indices

        Args:
            indices (np.ndarray): P space indices
            exactMeff (bool, optional): If true will calculate the exact Meff based on M. Defaults to False.

        Returns:
            CRM,float: Generated CRM object and the estimated accuracy/validity based on the Greenland norm of $T_Q^{-1}\delta$
        """

        Tp, Tq, bigDelta, smallDelta = decompT(self.eigenVecs, indices)

        TpInv = np.linalg.inv(Tp)
        TqInv = np.linalg.inv(Tq)
        Dp = np.diag(self.eigenVals[indices])

        if exactMeff:
            Mp, Mq, H, V = decompT(self.fullMat, indices)
            Meff = (np.diag(np.ones([len(indices)])) - H @ np.linalg.inv(Mq) @ smallDelta @ TpInv) @ Tp @ Dp @ TpInv
        else:
            Meff = Tp @ Dp @ TpInv

        dTqInv = -bigDelta @ TqInv
        dTpInv = smallDelta @ TpInv

        validityCoeff = greenlandMetric(TqInv @ smallDelta)

        return CRM(Meff, dTqInv, dTpInv), validityCoeff