# -*- coding: utf-8 -*-
import numpy as np
from sympy import symbols, Poly

def partition_to_array(mu):
    d = mu.as_dict()
    return np.repeat(list(d.keys()), list(d.values()))


def hook_lengths_gmp(mu, alpha):
    mu_prime = np.array(mu.conjugate)
    mu_ = partition_to_array(mu)
    i = np.repeat(np.arange(len(mu_)), mu_)
    j = np.concatenate([np.arange(n) for n in mu_])
    x = mu_prime[j] - i + alpha*(mu_[i] - j) 
    return (x - 1, x - alpha)


def _betaratio(kappa, mu, k, alpha):
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


def _N(kappa, mu):
    n = len(kappa)
    kappa = kappa + 1
    M = np.array([np.prod(kappa[(i+1):]) for i in range(n)])
    return np.sum(mu * M)


def constant_poly(a):
    x = symbols("x")
    p = Poly(a*x**0, x)
    return p.subs(x, 0)
