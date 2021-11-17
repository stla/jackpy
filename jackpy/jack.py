# -*- coding: utf-8 -*-
from gmpy2 import mpq, mpz, fac
import numpy as np
from sympy import symbols, Poly
from .internal import (
    partition_to_array,
    hook_lengths_gmp,
    _N,
    _betaratio,
    constant_poly,
    _is_number
)

def SchurPol(n, kappa):
    """
    Schur polynomial of an integer partition.

    Parameters
    ----------
    n : int
        Positive integer, the number of variables of the polynomial.
    kappa : IntegerPartition
        An integer partition obtained with `sympy.combinatorics.partitions`.

    Returns
    -------
    Poly
        The Schur polynomial of `kappa`. It has integer coefficents.
    
    Examples
    --------
    >>> from sympy.combinatorics.partitions import IntegerPartition
    >>> from jackpy.jack import SchurPol
    >>> p = SchurPol(2, IntegerPartition([2,1]))
    >>> print(p)
    Poly(x_1**2*x_0 + x_1*x_0**2, x_1, x_0, domain='ZZ')
    >>> y = p.eval({x_0: 1, x_1: 1})
    >>> print(y)
    2

    """
    if not isinstance(n, int):
        raise ValueError("`n` must be a integer.")
    if n < 1:
        raise ValueError("`n` must be at least one.")
    if not isinstance(kappa, IntegerPartition):
        raise ValueError("`kappa` must be a SymPy integer partition.")
    def sch(m, k, nu):
        if len(nu) == 0 or nu[0] == 0 or m == 0:
            return constant_poly(1)
        if len(nu) > m and nu[m] > 0:
            return constant_poly(0)
        if m == 1:
            return Poly(x[0]**nu[0], x[0])
        s = S[_N(kappa_, nu)-1, m-1]
        if s is not None:
            return s
        s = sch(m-1, 1, nu)
        i = k
        while len(nu) >= i and nu[i-1] > 0:
            if len(nu) == i or nu[i-1] > nu[i]:
                _nu = nu.copy()
                _nu[i-1] = nu[i-1]-1
                if nu[i-1] > 1:
                    s = s + Poly(x[m-1], x[m-1]) * sch(m, i, _nu)
                else:
                    s = s + Poly(x[m-1], x[m-1]) * sch(m-1, 1, _nu)
            i = i + 1
        if k == 1:
            S[_N(kappa_, nu)-1, m-1] = s
        return s
    x = [symbols(f'x_{i}') for i in range(n)]
    kappa_ = partition_to_array(kappa)
    S = np.full((_N(kappa_,kappa_), n), None)
    return sch(n, 1, kappa_)


def JackPol(n, kappa, alpha):
    """
    Jack polynomial of an integer partition.

    Parameters
    ----------
    n : int
        Positive integer, the number of variables of the polynomial.
    kappa : IntegerPartition
        An integer partition obtained with `sympy.combinatorics.partitions`.
    alpha : number
        A positive number, the parameter of the Jack polynomial.

    Returns
    -------
    Poly
        The Jack polynomial of `kappa` with parameter `alpha`. The type of 
        its coefficients depend on the type of `alpha`.
    
    Examples
    --------
    >>> from gmpy2 import mpq
    >>> from sympy.combinatorics.partitions import IntegerPartition
    >>> from jackpy.jack import JackPol
    >>>
    >>> poly = JackPol(3, IntegerPartition([2, 1]), alpha = mpq(3, 2))
    >>> print(poly)
    Poly(7/2*x_0**2*x_1 + 7/2*x_0**2*x_2 + 7/2*x_0*x_1**2 + 6*x_0*x_1*x_2
    + 7/2*x_0*x_2**2 + 7/2*x_1**2*x_2 + 7/2*x_1*x_2**2, x_0, x_1, x_2, domain='QQ')

    """
    if not isinstance(n, int):
        raise ValueError("`n` must be a integer.")
    if n < 1:
        raise ValueError("`n` must be at least one.")
    if not isinstance(kappa, IntegerPartition):
        raise ValueError("`kappa` must be a SymPy integer partition.")
    if not _is_number(alpha) and type(alpha) != type(mpz(0)) and type(alpha) != type(mpq(0)):
        raise ValueError("`alpha` must be a real number.")
    if alpha <= 0:
        raise ValueError("`alpha` must be positive.")
    def jac(m, k, mu, nu, beta):
        if len(nu) == 0 or nu[0] == 0 or m == 0:
            return constant_poly(1)
        if len(nu) > m and nu[m] > 0:
            return constant_poly(0)
        if m == 1:
            coef = np.prod(alpha * np.arange(1, nu[0]) + 1)
            return Poly(coef * x[0]**nu[0], x[0]) 
        s = S[_N(kappa_,nu)-1, m-1]
        if k == 0 and s is not None:
            return s
        i = max(1, k)
        s = (
            jac(m-1, 0, nu, nu, 1)
            * beta
            * Poly(x[m-1]**(np.sum(mu) - np.sum(nu)), x[m-1])
        )
        while len(nu) >= i and nu[i-1] > 0:
            if len(nu) == i or nu[i-1] > nu[i]:
                _nu = nu.copy()
                _nu[i-1] = nu[i-1]-1
                gamma = beta * _betaratio(mu, nu, i-1, alpha)
                if nu[i-1] > 1:
                    s = s + jac(m, i, mu, _nu, gamma)
                else:
                    s = s + jac(m-1, 0, _nu, _nu, 1) * gamma * Poly(
                        x[m-1]**(np.sum(mu) - np.sum(_nu)), x[m-1]
                    )
            i += 1
        if k == 0:
            S[_N(kappa_, nu)-1, m-1] = s
        return s
    x = [symbols(f'x_{i}') for i in range(n)]
    kappa_ = partition_to_array(kappa)
    S = np.full((_N(kappa_,kappa_), n), None)
    return jac(n, 0, kappa_, kappa_, 1)


def ZonalPol(m, kappa):
    """
    Zonal polynomial of an integer partition. Up to a normalization, this is 
    the Jack polynomial of the integer partition with paramater `alpha = 2`.

    Parameters
    ----------
    m : int
        Positive integer, the number of variables of the polynomial.
    kappa : IntegerPartition
        An integer partition obtained with `sympy.combinatorics.partitions`.

    Returns
    -------
    Poly
        The zonal polynomial of `kappa`. It has rational coefficients.
    
    """
    if not isinstance(m, int):
        raise ValueError("`m` must be a integer.")
    if m < 1:
        raise ValueError("`m` must be at least one.")
    if not isinstance(kappa, IntegerPartition):
        raise ValueError("`kappa` must be a SymPy integer partition.")
    alpha = mpq(2)
    jack = JackPol(m, kappa, alpha)
    (hooku, hookl) = hook_lengths_gmp(kappa, alpha)
    jlambda = np.prod(hooku) * np.prod(hookl)
    n = int(np.sum(partition_to_array(kappa)))
    return (alpha**n * fac(n) / jlambda) * jack


def ZonalQPol(m, kappa):
    """
    Quaternionic zonal polynomial of an integer partition. Up to a 
    normalization, this is the Jack polynomial of the integer partition with 
    paramater `alpha = 1/2`.

    Parameters
    ----------
    m : int
        Positive integer, the number of variables of the polynomial.
    kappa : IntegerPartition
        An integer partition obtained with `sympy.combinatorics.partitions`.

    Returns
    -------
    Poly
        The quaternionic zonal polynomial of `kappa`. 
        It has rational coefficients.
    
    """
    if not isinstance(m, int):
        raise ValueError("`m` must be a integer.")
    if m < 1:
        raise ValueError("`m` must be at least one.")
    if not isinstance(kappa, IntegerPartition):
        raise ValueError("`kappa` must be a SymPy integer partition.")
    alpha = mpq(1, 2)
    jack = JackPol(m, kappa, alpha)
    (hooku, hookl) = hook_lengths_gmp(kappa, alpha)
    jlambda = np.prod(hooku) * np.prod(hookl)
    n = int(np.sum(partition_to_array(kappa)))
    return (alpha**n * fac(n) / jlambda) * jack
