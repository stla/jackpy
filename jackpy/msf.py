# -*- coding: utf-8 -*-
from sympy.combinatorics import Permutation
from sympy.combinatorics.generators import symmetric
from sympy import symbols, Poly
import numpy as np
from .internal import (
    partition_to_array,
    constant_poly
)

def perms(mu):
    symgroup = symmetric(len(mu))
    allperms = np.array(
        [p(mu) for p in symgroup]    
    )
    return np.unique(allperms, axis=0)


def msf_poly(m, kappa):
    kappa_ = partition_to_array(kappa)
    n = len(kappa_)
    if n > m:
        return constant_poly(0)
    mu = np.zeros(m, dtype=int)
    mu[:n] = kappa_
    perms_mu = perms(mu)
    x = [symbols(f'x_{i}') for i in range(m)]
    monomials = np.array([np.prod(np.array(
        [Poly(x[i]**perm[i], x[i]) for i in range(m)]    
    )) for perm in perms_mu])
    return np.sum(monomials)

