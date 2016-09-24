#!/usr/bin/env python

"""Auxiliary functions for sympy matrix manipulations."""

from sympy import Matrix, log, exp, factor, expand, SparseMatrix, solve, sympify
from numpy import array, float64

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


def schur_reduction(A, rowIndices, colIndices):
    """Return the Schur complement of the block identified
    by rowIndices and colIndices.
    Works only if this block is square and invertible."""
    complRowIndices = [i for i in range(A.rows) if i not in rowIndices]
    complColIndices = [i for i in range(A.cols) if i not in colIndices]
    A11 = A.extract(complRowIndices, complColIndices)
    A12 = A.extract(complRowIndices, colIndices)
    A21 = A.extract(rowIndices, complColIndices)
    A22 = A.extract(rowIndices, colIndices)
    return (A11 - A12.multiply(A22.inv().multiply(A21))).applyfunc(lambda x: factor(expand(x)))


def negative(matrix):
    """Return a matrix with the negative elements of matrix."""
    return matrix.applyfunc(lambda x: x if x < 0 else 0)


def clog(A):
    """Componentwise logarithm."""
    return A.applyfunc(lambda x: log(x))


def cexp(A):
    """Componentwise exp."""
    return A.applyfunc(lambda x: exp(x))


def sdiag(vector):
    """Return a sparse diagonal matrix with vector as diagonal."""
    n = len(vector)
    return SparseMatrix(n, n, dict(((i, i), vector[i]) for i in range(n)))


def print_matrix(matrix, rowlabels, collabels, numeric = False, precision = 3):
    """Print a matrix with rowlabels as first column,
    and collabels as first row.
    If numeric = True, the matrix is converted to a numpy array,
    and precision can be used to format the elements (default is 3)."""
    if matrix != Matrix():
        # check that rowlabels and collabels are in the correct number
        nr = matrix.rows
        nc = matrix.cols
        if len(rowlabels) != matrix.rows:
            raise ValueError("Wrong number of row labels")
        if len(collabels) != matrix.cols:
            raise ValueError("Wrong number of column labels")
        if numeric:
            matrix = array(matrix).astype(float64)

        def to_str(x):
            if numeric:
                return '{:.{}f}'.format(x, precision)
            return str(x)

        rLength = max(map(len, rowlabels))
        cLengths = [max([len(collabels[i])] + [len(to_str(matrix[j, i])) for j in range(nr)]) for i in range(nc)]
        print(' ' * (rLength + 3) + '  '.join([collabels[i].rjust(cLengths[i]) for i in range(nc)]))
        for i in range(nr):
            print(rowlabels[i].ljust(rLength) + " | " + \
                  '  '.join(to_str(matrix[i,j]).rjust(cLengths[j]) for j in range(nc)) + ' |')


def dependent(matrix, vector, particular = False):
    """Check if vector can be written as a linear combination of
    the columns of matrix.
    Return the coeffients if possible, None otherwise.
    If particular = True, use the option 'particular' of sympy's solve(),
    i.e. try to find a particular solution with as many zeros as possible."""
    coeffs = Matrix(list(map(lambda i: sympify("x" + str(i)), range(len(matrix)))))
    sols = solve(Matrix(matrix).T * coeffs - Matrix(vector), coeffs, dict = True, particular = particular)
    if len(sols) > 0: return [sols[0][s] if s in sols[0] else 0 for s in coeffs]
    return None


def _pos_dependent(matrix, vector):
    """Check if vector can be written as a positive linear combination of
    the columns of matrix.
    Return the coeffients if possible, None otherwise."""
    #A = np.transpose(np.array(matrix).astype(np.float64))
    #b = np.array(vector).astype(np.float64)
    #x, rnorm = nnls(A, b)
    #return rnorm < 0.0000001, x
    if len(matrix) > 0:
        coeffs = Matrix(list(map(lambda i: sympify("x" + str(i)), range(len(matrix)))))
        sols = solve(Matrix(matrix).T * coeffs - Matrix(vector), coeffs, dict = True, particular = True)
        for sol in sols:
            if all([sol[s] >= 0 for s in sol]): return [sol[s] if s in sol else 0 for s in coeffs]
    return None


def _pos_generators(matrix):
    """Returns a matrix with rows the extreme rays of
    the pointed cone `matrix x = 0, x >= 0'."""
    import cdd

    if matrix == Matrix(): return matrix

    S = matrix
    nr = S.cols

    # create the matrix | S^t -S^t I|^t
    # with an additional first column of zeros for the constant term
    A = [[0] + [int(r[i]) for i in range(len(r))] for r in S.col_join(-S).tolist()] + \
        [[0 if j != i else 1 for j in range(nr + 1)] for i in range(nr + 1)]

    # create the polyhedron
    A = cdd.Matrix(A, number_type = "fraction")
    A.rep_type = cdd.RepType.INEQUALITY

    # find extreme rays
    ers = cdd.Polyhedron(A).get_generators()
    ers = [er[1:] for er in ers if any(e != 0 for e in er[1:])]

    return Matrix(ers)


def _cancel_opposite(a, b):
    m = min(a, b)
    return a - m, b - m
