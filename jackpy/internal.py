# -*- coding: utf-8 -*-
import numpy as np
from sympy import symbols, Poly
from numbers import Rational, Real

def __get_domain__(x):
    if isinstance(x, Rational):
        return 'QQ'
    elif isinstance(x, Real):
        return 'RR'
    else:
        return None    


def __partition_to_array__(mu):
    d = mu.as_dict()
    if len(d) == 0:
        return np.asarray([], dtype=int)
    return np.repeat(list(d.keys()), list(d.values()))


def __hook_lengths_gmp__(mu, alpha):
    mu_prime = np.array(mu.conjugate)
    mu_ = __partition_to_array__(mu)
    i = np.repeat(np.arange(len(mu_)), mu_)
    j = np.concatenate([np.arange(n) for n in mu_])
    x = mu_prime[j] - i + alpha*(mu_[i] - j) 
    return (x - 1, x - alpha)


def __betaratio__(kappa, mu, k, alpha):
    k += 1
    t = k - alpha*mu[k-1] 
    s = np.arange(1, k+1)
    u = t + 1 - s + alpha*kappa[s-1]
    s = np.arange(1, k)
    v = t - s + alpha*mu[s-1]
    s = np.arange(1, mu[k-1])
    w = np.array([np.count_nonzero(mu >= i) for i in s]) - t - alpha*s
    return (
        alpha
        * np.prod(u / (u + alpha - 1))
        * np.prod((v + alpha) / v)
        * np.prod((w + alpha) / w)
    )


def __N__(kappa, mu):
    n = len(kappa)
    if n == 0:
        return 0
    kappa = kappa + 1
    M = np.array([np.prod(kappa[(i+1):]) for i in range(n)])
    return np.sum(mu * M)
