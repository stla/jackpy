# -*- coding: utf-8 -*-
from sympy.combinatorics.partitions import IntegerPartition
from sympy import symbols, Poly, parse_expr
import numpy as np
from .internal import (
        __partition_to_array__
    ,   __permutations__
    ,   __drop_trailing_zeros__
    ,   __is_decreasing__
    ,   __msp_symbol__
    )

def monomial_symmetric_polynomial(n, kappa):
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
    if l == 0:
        return Poly(1, *variables, domain='ZZ')
    mu = np.zeros(n, dtype=int)
    mu[:l] = kappa_
    perms_mu = __permutations__(mu)
    x = [Poly(v, *variables, domain='ZZ') for v in variables]
    terms = np.array([np.prod(np.array(
        [x[i]**perm[i] for i in range(n)]    
    )) for perm in perms_mu])
    return np.sum(terms)

def msp_combination(poly):
    """
    Symmetric polynomial as a linear combination of some monomial symmetric polynomials. 

    Parameters
    ----------
    p : Poly
        Polynomial. It must be symmetric, otherwise the output of this function makes no sense. 

    Returns
    -------
    dict
        A dictionary representing the linear combination. A key represents
        the integer partition of a monomial symmetric polynomial, and the value
        attached to this key is the coefficient of this monomial symmetric 
        polynomial.

    """
    out = {}
    d = poly.as_dict()
    for exponents in d.keys():
        if __is_decreasing__(exponents):
            kappa = __drop_trailing_zeros__(list(exponents))
            out[tuple(kappa)] = d.get(exponents)
    return out

def msp_combination_expr(poly):
    """
    Symmetric polynomial as a linear combination of some monomial symmetric polynomials. 

    Parameters
    ----------
    p : Poly
        Polynomial. It must be symmetric, otherwise the output of this function makes no sense. 

    Returns
    -------
    expression
        An expression representing the linear combination. The monomial 
        symmetric polynomials are represented by symbols, e.g. the monomial 
        symmetric polynomial of the integer partition `(3,2,1)` is represented 
        by the symbol `M[3;2;1]`.

    """
    combo = msp_combination(poly)
    kappas = combo.keys();
    out = parse_expr("0")
    for kappa in kappas:
        out = out + combo.get(kappa) * __msp_symbol__(kappa)
    return out
