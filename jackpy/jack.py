# -*- coding: utf-8 -*-
from gmpy2 import mpq
import numpy as np
from sympy import symbols, Poly, Symbol
from .internal import (
    __get_domain__,
    __make_partition__,
    __N__,
    __betaratio__,
    __Jack_C_coefficient__,
    __Jack_P_coefficient__,
    __Jack_Q_coefficient__
)
from numbers import Real, Number, Integral 

def SchurPol(n, kappa):
    """
    Schur polynomial of an integer partition.

    Parameters
    ----------
    n : int
        Positive integer, the number of variables of the polynomial.
    kappa : list of integers
        An integer partition given as a list of decreasing integers. Trailing 
        zeros are dropped.

    Returns
    -------
    Poly
        The Schur polynomial of `kappa` in `n` variables `x_1`, ..., `x_n`, 
        with integer coefficents.
    
    Examples
    --------
    >>> from jackpy.jack import SchurPol
    >>> p = SchurPol(2, [2, 1])
    >>> print(p)
    Poly(x_1*x_2**2 + x_1**2*x_2, x_1, x_2, domain='ZZ')
    >>> y = p.eval({x_1: 1, x_2: 1})
    >>> print(y)
    2

    """
    if not (isinstance(n, int) and n >= 1):
        raise ValueError("`n` must be a strictly positive integer.")
    kappa_ = __make_partition__(kappa)
    variables = [symbols(f'x_{i}') for i in range(1, n+1)]
    x = [Poly(v, *variables, domain='ZZ') for v in variables]
    def sch(S, m, k, nu):
        if len(nu) == 0 or nu[0] == 0 or m == 0:
            return Poly(1, *variables, domain='ZZ')
        if len(nu) > m and nu[m] > 0:
            return Poly(0, *variables, domain='ZZ')
        if m == 1:
            return x[0]**nu[0]
        N = __N__(kappa_, nu)
        s = S[N-1, m-1]
        if s is not None:
            return s
        s = sch(S.copy(), m-1, 1, nu)
        i = k
        while len(nu) >= i and nu[i-1] > 0:
            if len(nu) == i or nu[i-1] > nu[i]:
                _nu = nu.copy()
                _nu[i-1] = nu[i-1]-1
                if nu[i-1] > 1:
                    s = s + x[m-1] * sch(S.copy(), m, i, _nu)
                else:
                    s = s + x[m-1] * sch(S.copy(), m-1, 1, _nu)
            i = i + 1
        if k == 1:
            S[__N__(kappa_, nu)-1, m-1] = s
        return s
    S = np.full((__N__(kappa_,kappa_), n), None)
    return sch(S.copy(), n, 1, kappa_)


def JackPol(n, kappa, alpha, which = 'J'):
    """
    Jack polynomial of an integer partition, with given Jack parameter.

    Parameters
    ----------
    n : int
        Positive integer, the number of variables of the polynomial.
    kappa : list of integers
        An integer partition given as a list of decreasing integers. Trailing 
        zeros are dropped.
    alpha : number
        A positive number, the parameter of the Jack polynomial.
    which: character
        Which Jack polynomial, either `'J'`, `'C'`, `'P'` or `'Q'`.

    Returns
    -------
    Poly
        The Jack polynomial of `kappa` in `n` variables `x_1`, ..., `x_n`, 
        with Jack parameter `alpha`. The type of 
        its coefficients depends on the type of `alpha`.
    
    Examples
    --------
    >>> from gmpy2 import mpq
    >>> from jackpy.jack import JackPol
    >>>
    >>> poly = JackPol(3, [2, 1], alpha = mpq(3, 2))
    >>> print(poly)
    Poly(7/2*x_1**2*x_2 + 7/2*x_1**2*x_3 + 7/2*x_1*x_2**2 + 6*x_1*x_2*x_3
    + 7/2*x_1*x_3**2 + 7/2*x_2**2*x_3 + 7/2*x_2*x_3**2, x_1, x_2, x3, domain='QQ')

    """
    if not (isinstance(n, int) and n >= 1):
        raise ValueError("`n` must be a strictly positive integer.")
    kappa_ = __make_partition__(kappa)
    if isinstance(alpha, Number):
        if not isinstance(alpha, Real):
            raise ValueError("`alpha` must be a real number.")
        if alpha <= 0:
            raise ValueError("`alpha` must be positive.")
        if isinstance(alpha, Integral):
            alpha = mpq(alpha)
        domain = __get_domain__(alpha)
    elif isinstance(alpha, Symbol):
        domain = 'QQ(alpha)'
    else:
        raise ValueError("`alpha` must be a number.")
    if not which in ['J', 'C', 'P', 'Q']:
        raise ValueError("`which` must be either 'J', 'C', 'P' or 'Q'.")
    variables = [symbols(f'x_{i}') for i in range(1, n+1)]
    x = [Poly(v, *variables, domain=domain) for v in variables]
    def jac(S, m, k, mu, nu, beta):
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
            jac(S.copy(), m-1, 0, nu, nu, 1)
            * beta
            * x[m-1]**(np.sum(mu) - np.sum(nu))
        )
        while len(nu) >= i and nu[i-1] > 0:
            if len(nu) == i or nu[i-1] > nu[i]:
                _nu = nu.copy()
                _nu[i-1] = nu[i-1]-1
                gamma = beta * __betaratio__(mu, nu, i-1, alpha)
                if nu[i-1] > 1:
                    s = s + jac(S.copy(), m, i, mu, _nu, gamma)
                else:
                    s = (
                        s + jac(S.copy(), m-1, 0, _nu, _nu, 1) * gamma
                        * x[m-1]**(np.sum(mu) - np.sum(_nu))
                    )
            i += 1
        if k == 0:
            S[__N__(kappa_, nu)-1, m-1] = s
        return s
    S = np.full((__N__(kappa_,kappa_), n), None)
    jp = jac(S.copy(), n, 0, kappa_, kappa_, 1)
    if which != 'J':
        if which == 'C':
            jp = __Jack_C_coefficient__(kappa_, alpha) * jp
        elif which == 'P':
            jp = __Jack_P_coefficient__(kappa_, alpha) * jp
        else:
            jp = __Jack_Q_coefficient__(kappa_, alpha) * jp
    return jp


def ZonalPol(n, kappa):
    """
    Zonal polynomial of an integer partition. This is the Jack C-polynomial of 
    this integer partition with parameter `alpha = 2`.

    Parameters
    ----------
    n : int
        Positive integer, the number of variables of the polynomial.
    kappa : list of integers
        An integer partition given as a list of decreasing integers. Trailing 
        zeros are dropped.

    Returns
    -------
    Poly
        The zonal polynomial of `kappa` in `n` variables
        `x_1`, ..., `x_n`, with rational coefficients.
    
    """
    return JackPol(n, kappa, mpq(2), which = 'C')


def ZonalQPol(n, kappa):
    """
    Quaternionic zonal polynomial of an integer partition. This is the Jack 
    C-polynomial of the integer partition with parameter `alpha = 1/2`.

    Parameters
    ----------
    n : int
        Positive integer, the number of variables of the polynomial.
    kappa : list of integers
        An integer partition given as a list of decreasing integers. Trailing 
        zeros are dropped.

    Returns
    -------
    Poly
        The quaternionic zonal polynomial of `kappa` in `n` variables
        `x_1`, ..., `x_n`, with rational coefficients.
    
    """
    return JackPol(n, kappa, mpq(1, 2), which = 'C')

def JackSymbolicPol(n, kappa, which = 'J'):
    """
    Jack polynomial of an integer partition, with symbolic Jack parameter.

    Parameters
    ----------
    n : int
        Positive integer, the number of variables of the polynomial.
    kappa : IntegerPartition
        An integer partition obtained with `sympy.combinatorics.partitions`.
    which: character
        Which Jack polynomial, either `'J'`, `'C'`, `'P'` or `'Q'`.

    Returns
    -------
    Poly
        The Jack polynomial of `kappa` in `n` variables `x_1`, ..., `x_n`, 
        with symbolic Jack parameter denoted by `alpha`. The domain of 
        this polynomial is `'QQ(alpha)'`.
    
    Examples
    --------
    >>> from sympy.combinatorics.partitions import IntegerPartition
    >>> from jackpy.jack import JackSymbolicPol
    >>>
    >>> poly = JackSymbolicPol(2, IntegerPartition([2, 1]))
    >>> print(poly)
    Poly((alpha + 2)*x_1**2*x_2 + (alpha + 2)*x_1*x_2**2, x_1, x_2, domain='QQ(alpha)')

    """
    return JackPol(n, kappa, symbols("alpha"), which)
