# -*- coding: utf-8 -*-
from gmpy2 import mpq
from sympy.combinatorics.partitions import IntegerPartition
from jackpy.jack import JackPol
from jackpy.msf import msf_poly

def test_jackpol():
    mu = IntegerPartition([3,1])
    alpha = mpq(5, 2)
    m = 4
    expected = (
        (2*alpha**2 + 4*alpha + 2) * msf_poly(m, IntegerPartition([3,1]))
        + (6*alpha + 10) * msf_poly(m, IntegerPartition([2,1,1]))
        + (4*alpha + 4) * msf_poly(m, IntegerPartition([2,2]))
        + 24 * msf_poly(m, IntegerPartition([1,1,1,1]))
    )
    obtained = JackPol(m, mu, alpha)
    assert obtained == expected