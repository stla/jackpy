# -*- coding: utf-8 -*-
from gmpy2 import mpq
import numpy as np
from sympy import symbols, Poly, Symbol
from sympy.combinatorics.partitions import IntegerPartition
from jackpy.internal import (
    __get_domain__,
    __partition_to_array__,
    __N__,
    __betaratio__,
    __Jack_C_coefficient__
)
from numbers import Real, Number 

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
        The Schur polynomial of `kappa` in `n` variables `x_1`, ..., `x_n`, 
        with integer coefficents.
    
    Examples
    --------
    >>> from sympy.combinatorics.partitions import IntegerPartition
    >>> from jackpy.jack import SchurPol
    >>> p = SchurPol(2, IntegerPartition([2,1]))
    >>> print(p)
    Poly(x_1*x_2**2 + x_1**2*x_2, x_1, x_2, domain='ZZ')
    >>> y = p.eval({x_1: 1, x_2: 1})
    >>> print(y)
    2

    """
    if not (isinstance(n, int) and n >= 1):
        raise ValueError("`n` must be a strictly positive integer.")
    if not isinstance(kappa, IntegerPartition):
        raise ValueError("`kappa` must be a SymPy integer partition.")
    variables = [symbols(f'x_{i}') for i in range(1, n+1)]
    x = [Poly(v, *variables, domain='ZZ') for v in variables]
    def sch(m, k, nu):
        if len(nu) == 0 or nu[0] == 0 or m == 0:
            return Poly(1, *variables, domain='ZZ')
        if len(nu) > m and nu[m] > 0:
            return Poly(0, *variables, domain='ZZ')
        if m == 1:
            return x[0]**nu[0]
        s = S[__N__(kappa_, nu)-1, m-1]
        if s is not None:
            return s
        s = sch(m-1, 1, nu)
        i = k
        while len(nu) >= i and nu[i-1] > 0:
            if len(nu) == i or nu[i-1] > nu[i]:
                _nu = nu.copy()
                _nu[i-1] = nu[i-1]-1
                if nu[i-1] > 1:
                    s = s + x[m-1] * sch(m, i, _nu)
                else:
                    s = s + x[m-1] * sch(m-1, 1, _nu)
            i = i + 1
        if k == 1:
            S[__N__(kappa_, nu)-1, m-1] = s
        return s
    kappa_ = __partition_to_array__(kappa)
    S = np.full((__N__(kappa_,kappa_), n), None)
    return sch(n, 1, kappa_)


def JackPol(n, kappa, alpha):
    """
    Jack polynomial of an integer partition, with given Jack parameter.

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
        The Jack polynomial of `kappa` in `n` variables `x_1`, ..., `x_n`, 
        with Jack parameter `alpha`. The type of 
        its coefficients depends on the type of `alpha`.
    
    Examples
    --------
    >>> from gmpy2 import mpq
    >>> from sympy.combinatorics.partitions import IntegerPartition
    >>> from jackpy.jack import JackPol
    >>>
    >>> poly = JackPol(3, IntegerPartition([2, 1]), alpha = mpq(3, 2))
    >>> print(poly)
    Poly(7/2*x_1**2*x_2 + 7/2*x_1**2*x_3 + 7/2*x_1*x_2**2 + 6*x_1*x_2*x_3
    + 7/2*x_1*x_3**2 + 7/2*x_2**2*x_3 + 7/2*x_2*x_3**2, x_1, x_2, x3, domain='QQ')

    """
    if not (isinstance(n, int) and n >= 1):
        raise ValueError("`n` must be a strictly positive integer.")
    if not isinstance(kappa, IntegerPartition):
        raise ValueError("`kappa` must be a SymPy integer partition.")
    if isinstance(alpha, Number):
        if not isinstance(alpha, Real):
            raise ValueError("`alpha` must be a real number.")
        if alpha <= 0:
            raise ValueError("`alpha` must be positive.")
        domain = __get_domain__(alpha)
    elif isinstance(alpha, Symbol):
        domain = 'QQ(alpha)'
    else:
        raise ValueError("`alpha` must be a number.")
    variables = [symbols(f'x_{i}') for i in range(1, n+1)]
    x = [Poly(v, *variables, domain=domain) for v in variables]
    def jac(m, k, mu, nu, beta):
        if len(nu) == 0 or nu[0] == 0 or m == 0:
            return Poly(1, *variables, domain=domain)
        if len(nu) > m and nu[m] > 0:
            return Poly(0, *variables, domain=domain)
        if m == 1:
            coef = np.prod(alpha * np.arange(1, nu[0]) + 1)
            return coef * x[0]**nu[0]
        s = S[__N__(kappa_,nu)-1, m-1]
        if k == 0 and s is not None:
            return s
        i = max(1, k)
        s = (
            jac(m-1, 0, nu, nu, 1)
            * beta
            * x[m-1]**(np.sum(mu) - np.sum(nu))
        )
        while len(nu) >= i and nu[i-1] > 0:
            if len(nu) == i or nu[i-1] > nu[i]:
                _nu = nu.copy()
                _nu[i-1] = nu[i-1]-1
                gamma = beta * __betaratio__(mu, nu, i-1, alpha)
                if nu[i-1] > 1:
                    s = s + jac(m, i, mu, _nu, gamma)
                else:
                    s = (
                        s + jac(m-1, 0, _nu, _nu, 1) * gamma
                        * x[m-1]**(np.sum(mu) - np.sum(_nu))
                    )
            i += 1
        if k == 0:
            S[__N__(kappa_, nu)-1, m-1] = s
        return s
    kappa_ = __partition_to_array__(kappa)
    S = np.full((__N__(kappa_,kappa_), n), None)
    return jac(n, 0, kappa_, kappa_, 1)


def ZonalPol(n, kappa):
    """
    Zonal polynomial of an integer partition. Up to a normalization, this is 
    the Jack polynomial of the integer partition with paramater `alpha = 2`.

    Parameters
    ----------
    n : int
        Positive integer, the number of variables of the polynomial.
    kappa : IntegerPartition
        An integer partition obtained with `sympy.combinatorics.partitions`.

    Returns
    -------
    Poly
        The zonal polynomial of `kappa` in `n` variables
        `x_1`, ..., `x_n`, with rational coefficients.
    
    """
    alpha = mpq(2)
    jack = JackPol(n, kappa, alpha)
    return __Jack_C_coefficient__(kappa, alpha) * jack


def ZonalQPol(n, kappa):
    """
    Quaternionic zonal polynomial of an integer partition. Up to a 
    normalization, this is the Jack polynomial of the integer partition with 
    paramater `alpha = 1/2`.

    Parameters
    ----------
    n : int
        Positive integer, the number of variables of the polynomial.
    kappa : IntegerPartition
        An integer partition obtained with `sympy.combinatorics.partitions`.

    Returns
    -------
    Poly
        The quaternionic zonal polynomial of `kappa` in `n` variables
        `x_1`, ..., `x_n`, with rational coefficients.
    
    """
    alpha = mpq(1, 2)
    jack = JackPol(n, kappa, alpha)
    return __Jack_C_coefficient__(kappa, alpha) * jack

def JackSymbolicPol(n, kappa):
    """
    Jack polynomial of an integer partition, with symbolic Jack parameter.

    Parameters
    ----------
    n : int
        Positive integer, the number of variables of the polynomial.
    kappa : IntegerPartition
        An integer partition obtained with `sympy.combinatorics.partitions`.

    Returns
    -------
    Poly
        The Jack polynomial of `kappa` in `n` variables `x_1`, ..., `x_n`, 
        with symbolic Jack parameter denoted by `alpha`. The domain of 
        this polynomial is `'QQ(alpha)`.
    
    Examples
    --------
    >>> from gmpy2 import mpq
    >>> from sympy.combinatorics.partitions import IntegerPartition
    >>> from jackpy.jack import JackSymbolicPol
    >>>
    >>> poly = JackSymbolicPol(3, IntegerPartition([2, 1]))
    >>> print(poly)
    Poly(7/2*x_1**2*x_2 + 7/2*x_1**2*x_3 + 7/2*x_1*x_2**2 + 6*x_1*x_2*x_3
    + 7/2*x_1*x_3**2 + 7/2*x_2**2*x_3 + 7/2*x_2*x_3**2, x_1, x_2, x3, domain='QQ')

    """
    return JackPol(n, kappa, symbols("alpha"))
