
import numpy as np
import tomograph



####################################################################################################
# Exercise 1: Gaussian elimination

def gaussian_elimination(A: np.ndarray, b: np.ndarray, use_pivoting: bool = True) -> (np.ndarray, np.ndarray):
    """
    Gaussian Elimination of Ax=b with or without pivoting.

    Arguments:
    A : matrix, representing left side of equation system of size: (m,m)
    b : vector, representing right hand side of size: (m, )
    use_pivoting : flag if pivoting should be used

    Return:
    A : reduced result matrix in row echelon form (type: np.ndarray, size: (m,m))
    b : result vector in row echelon form (type: np.ndarray, size: (m, ))

    Raised Exceptions:
    ValueError: if matrix and vector sizes are incompatible, matrix is not square or pivoting is disabled but necessary

    Side Effects:
    -

    Forbidden:
    - numpy.linalg.*
    """
    # Create copies of input matrix and vector to leave them unmodified
    A = A.copy()
    b = b.copy()

    # TODO: Test if shape of matrix and vector is compatible and raise ValueError if not
    (n, m) = A.shape
    o = b.shape[0]
    if m != n:
        raise ValueError("matrix not quadratic")
    if m != o:
        raise ValueError("matrix and vector are incompatible")

    # TODO: Perform gaussian elimination

   # piv = np.ones(n)
    #for i in range(0, n):
    #    for j in range(0, n):
    #        if A[i, j] != 0:
    #            piv[i] = j
    #            break
    #for i in range(0, n-1):
    #    for j in range(i, n):
    #        if piv[i] > piv[j]:
    #            raise ValueError("nicht ohne pivoting lösbar")



    if use_pivoting == False:
        for i in range(n):
            div = A[i][i]
            if div == 0:                                        #TODO NEW
                raise ValueError("not solvable without pivot")  #TODO NEW

            for j in range(m):
                A[i][j] = A[i][j] / div
            b[i] = b[i] / div

            for k in range(i+1, n):
                mult = A[k][i]
                for j in range(m):
                    A[k][j] = A[k][j] - mult * A[i][j]
                b[k] = b[k] - mult * b[i]

    elif use_pivoting == True:
        for i in range(n):
            gi = i
            for j in range(i, n):
                if abs(A[j][i]) > abs(A[i][i]):
                    A[i, j] = A[j, i]
                    b[i], b[j] = b[j], b[i]
            div = A[i][i]

            for j in range(m):
                A[i][j] = A[i][j] / div
            b[i] = b[i] / div

            for k in range(i+1, n):
                mult = A[k][i]
                for j in range(m):
                    A[k][j] = A[k][j] - mult * A[i][j]
                b[k] = b[k] - mult * b[i]

    return A, b


def back_substitution(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Back substitution for the solution of a linear system in row echelon form.

    Arguments:
    A : matrix in row echelon representing linear system
    b : vector, representing right hand side

    Return:
    x : solution of the linear system

    Raised Exceptions:
    ValueError: if matrix/vector sizes are incompatible or no/infinite solutions exist

    Side Effects:
    -

    Forbidden:
    - numpy.linalg.*
    """

    # TODO: Test if shape of matrix and vector is compatible and raise ValueError if not



    if A.shape[1] != b.shape[0]:
        raise ValueError("matrix and vector are incompatible")

    # TODO: Initialize solution vector with proper size

    x = np.zeros_like(b)

    # TODO: Run backsubstitution and fill solution vector, raise ValueError if no/infinite solutions exist
    if A[b.shape[0]-1][b.shape[0]-1] == 0:
        raise ValueError("no solution")
    for k in range(A.shape[1] - 1, -1, -1):
        x[k] = (b[k] - np.dot(A[k, k + 1:], x[k + 1:])) / A[k, k]
    return x


####################################################################################################
# Exercise 2: Cholesky decomposition

def compute_cholesky(M: np.ndarray) -> np.ndarray:
    """
    Compute Cholesky decomposition of a matrix

    Arguments:
    M : matrix, symmetric and positive (semi-)definite

    Raised Exceptions:
    ValueError: L is not symmetric and psd

    Return:
    L :  Cholesky factor of M

    Forbidden:
    - numpy.linalg.*
    """

    # TODO check for symmetry and raise an exception of type ValueError
    (n, m) = M.shape

    if m != n:
        raise ValueError("not quadratic")
    if not np.allclose(M, M.T):
        raise ValueError("not symmetric")
    #if not np.all(np.dot(M, M.T) == np.dot(M, M)):
    #    raise ValueError("not a PSD")




    # TODO build the factorization and raise a ValueError in case of a non-positive definite input matrix
    # TODO PSD

    L = np.zeros((n, n))
    try:
        for i in range(n):
            for k in range(i+1):
                tmp = sum(L[i][j] * L[k][j] for j in range(k))
                if i == k:
                    L[i][k] = np.sqrt(M[i][i] - tmp)
                else:
                    L[i][k] = (1.0/L[k][k]*(M[i][k] - tmp))
    except RuntimeWarning:
        raise ValueError("not a PSD")

    return L


def solve_cholesky(L: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Solve the system L L^T x = b where L is a lower triangular matrix

    Arguments:
    L : matrix representing the Cholesky factor
    b : right hand side of the linear system

    Raised Exceptions:
    ValueError: sizes of L, b do not match
    ValueError: L is not lower triangular matrix

    Return:
    x : solution of the linear system

    Forbidden:
    - numpy.linalg.*
    """

    # TODO Check the input for validity, raising a ValueError if this is not the case
    (n, m) = L.shape
    if n != m:
        raise ValueError("not a quadratic matrix")
    if (b.shape[0] != m) or (b.shape[0] != n):
        raise ValueError("matrix not compatible")
    if not np.allclose(np.tril(L), L):
        raise ValueError("not a triangular matrix")
    # TODO Solve the system by forward- and backsubstitution
    x = np.zeros(m)
    y = np.zeros(m)

    for i in range(m):
        y[i] = b[i]
        for j in range(i):
            y[i] = y[i] - L[i, j] * y[j]
        y[i] = y[i] / L[i][i]

    x = back_substitution(L.T, y)
    return x


####################################################################################################
# Exercise 3: Tomography

def setup_system_tomograph(n_shots: np.int, n_rays: np.int, n_grid: np.int) -> (np.ndarray, np.ndarray):
    """
    Set up the linear system describing the tomographic reconstruction

    Arguments:
    n_shots  : number of different shot directions
    n_rays   : number of parallel rays per direction
    n_grid   : number of cells of grid in each direction, in total n_grid*n_grid cells

    Return:
    L : system matrix
    g : measured intensities

    Raised Exceptions:
    -

    Side Effects:
    -

    Forbidden:
    -
    """

    # TODO: Initialize system matrix with proper size
    L = np.zeros((n_shots * n_rays, n_grid * n_grid))
    # TODO: Initialize intensity vector
    g = np.zeros(n_rays * n_shots)

    # TODO: Iterate over equispaced angles, take measurements, and update system matrix and sinogram
    theta = 0
    # Take a measurement with the tomograph from direction r_theta.
    # intensities: measured intensities for all <n_rays> rays of the measurement. intensities[n] contains the intensity for the n-th ray
    # ray_indices: indices of rays that intersect a cell
    # isect_indices: indices of intersected cells
    # lengths: lengths of segments in intersected cells
    # The tuple (ray_indices[n], isect_indices[n], lengths[n]) stores which ray has intersected which cell with which length. n runs from 0 to the amount of ray/cell intersections (-1) of this measurement.
    for i in range(n_shots):
        theta += np.pi/n_shots
        intensities, ray_indices, isect_indices, lengths = tomograph.take_measurement(n_grid, n_rays, theta)
        g[i] = intensities
        L[i] = (ray_indices,isect_indices, lengths)
    return [L, g]


def compute_tomograph(n_shots: np.int, n_rays: np.int, n_grid: np.int) -> np.ndarray:
    """
    Compute tomographic image

    Arguments:
    n_shots  : number of different shot directions
    n_rays   : number of parallel rays per direction
    n_grid   : number of cells of grid in each direction, in total n_grid*n_grid cells

    Return:
    tim : tomographic image

    Raised Exceptions:
    -

    Side Effects:
    -

    Forbidden:
    """

    # Setup the system describing the image reconstruction
    [L, g] = setup_system_tomograph(n_shots, n_rays, n_grid)

    # TODO: Solve for tomographic image using your Cholesky solver
    # (alternatively use Numpy's Cholesky implementation)

    # TODO: Convert solution of linear system to 2D image
    tim = np.zeros((n_grid, n_grid))

    return tim


if __name__ == '__main__':
    print("All requested functions for the assignment have to be implemented in this file and uploaded to the "
          "server for the grading.\nTo test your implemented functions you can "
          "implement/run tests in the file tests.py (> python3 -v test.py [Tests.<test_function>]).")
