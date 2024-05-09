# -*- coding: utf-8 -*-
from sympy.combinatorics.partitions import IntegerPartition
from sympy.combinatorics.generators import symmetric
from sympy import symbols, Poly
import numpy as np
from .internal import __partition_to_array__

def perms(mu):
    symgroup = symmetric(len(mu))
    allperms = np.array(
        [p(mu) for p in symgroup]    
    )
    return np.unique(allperms, axis=0)


def msf_poly(n, kappa):
    """
    Monomial symmetric polynomial. 

    Parameters
    ----------
    n : int
        Positive integer, the number of variables of the polynomial.
    kappa : IntegerPartition
        An integer partition obtained with `sympy.combinatorics.partitions`.

    Returns
    -------
    Poly
        The monomial symmetric polynomial corresponding to `kappa` in `n` 
        variables `x_1`, ..., `x_n`, with integer coefficients.

    """
    if not (isinstance(n, int) and n >= 1):
        raise ValueError("`n` must be a strictly positive integer.")
    if not isinstance(kappa, IntegerPartition):
        raise ValueError("`kappa` must be a SymPy integer partition.")
    variables = [symbols(f'x_{i}') for i in range(1, n+1)]
    kappa_ = __partition_to_array__(kappa)
    l = len(kappa_)
    if l > n:
        return Poly(0, *variables, domain='ZZ')
    mu = np.zeros(n, dtype=int)
    mu[:l] = kappa_
    perms_mu = perms(mu)
    x = [Poly(v, *variables, domain='ZZ') for v in variables]
    terms = np.array([np.prod(np.array(
        [x[i]**perm[i] for i in range(n)]    
    )) for perm in perms_mu])
    return np.sum(terms)

