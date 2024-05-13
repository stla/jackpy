# -*- coding: utf-8 -*-
from gmpy2 import fac
import numpy as np
from numbers import Rational, Real
from sympy.utilities.iterables import multiset_permutations

def __is_decreasing__(x):
    l = len(x)
    i = 0
    out = True
    while i < l-1 and out:
        out = x[i] >= x[i+1]
        i = i + 1
    return out

def __drop_trailing_zeros__(x):
    n = len(x) - 1
    while n >= 0 and x[n] == 0:
        x.pop()
        n = n - 1
    return x

def __permutations__(mu):
    return list(multiset_permutations(mu))

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

def __hook_lengths_lower__(mu, alpha):
    mu_prime = np.array(mu.conjugate)
    mu_ = __partition_to_array__(mu)
    i = np.repeat(np.arange(len(mu_)), mu_)
    j = np.concatenate([np.arange(n) for n in mu_])
    x = mu_prime[j] - i + alpha*(mu_[i] - j) 
    return x - alpha

def __hook_lengths_upper__(mu, alpha):
    mu_prime = np.array(mu.conjugate)
    mu_ = __partition_to_array__(mu)
    i = np.repeat(np.arange(len(mu_)), mu_)
    j = np.concatenate([np.arange(n) for n in mu_])
    x = mu_prime[j] - i + alpha*(mu_[i] - j) 
    return x - 1

def __hook_lengths__(mu, alpha):
    mu_prime = np.array(mu.conjugate)
    mu_ = __partition_to_array__(mu)
    i = np.repeat(np.arange(len(mu_)), mu_)
    j = np.concatenate([np.arange(n) for n in mu_])
    x = mu_prime[j] - i + alpha*(mu_[i] - j) 
    return (x - alpha, x - 1)

def __Jack_C_coefficient__(kappa, alpha):
    (hookl, hooku) = __hook_lengths__(kappa, alpha)
    jlambda = np.prod(hooku) * np.prod(hookl)
    k = int(np.sum(__partition_to_array__(kappa)))
    return alpha**k * fac(k) / jlambda

def __Jack_P_coefficient__(kappa, alpha):
    if len(kappa.as_dict()) == 0:
        return 1
    hookl = __hook_lengths_lower__(kappa, alpha)
    return 1 / np.prod(hookl)
    
def __Jack_Q_coefficient__(kappa, alpha):
    if len(kappa.as_dict()) == 0:
        return 1
    hooku = __hook_lengths_upper__(kappa, alpha)
    return 1 / np.prod(hooku)


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
    nu = kappa + 1
    M = np.array([np.prod(nu[(i+1):]) for i in range(n)])
    return np.sum(mu * M)
