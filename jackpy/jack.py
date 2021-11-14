# -*- coding: utf-8 -*-
from gmpy2 import mpq, fac
#from sympy.combinatorics.partitions import IntegerPartition
import numpy as np
from sympy import symbols, Poly
from .internal import (
    partition_to_array,
    hook_lengths_gmp,
    _N,
    _betaratio,
    constant_poly
)

def SchurPol(n, kappa):
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
    #stopifnot(isPositiveInteger(n), alpha >= 0, isPartition(lambda))
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
    alpha = mpq(2)
    jack = JackPolDK(m, kappa, alpha)
    (hooku, hookl) = hook_lengths_gmp(kappa, alpha)
    jlambda = np.prod(hooku) * np.prod(hookl)
    n = int(np.sum(partition_to_array(kappa)))
    return (alpha**n * fac(n) / jlambda) * jack


def ZonalQPol(m, kappa):
    alpha = mpq(1, 2)
    jack = JackPolDK(m, kappa, alpha)
    (hooku, hookl) = hook_lengths_gmp(kappa, alpha)
    jlambda = np.prod(hooku) * np.prod(hookl)
    n = int(np.sum(partition_to_array(kappa)))
    return (alpha**n * fac(n) / jlambda) * jack
