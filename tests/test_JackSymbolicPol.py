# -*- coding: utf-8 -*-
from sympy import symbols
from sympy.combinatorics.partitions import IntegerPartition
from jackpy.jack import JackSymbolicPol
from jackpy.monomial_symmetric import monomial_symmetric_polynomial

def test_jacksymbolicpol():
    mu = IntegerPartition([3,1])
    alpha = symbols("alpha")
    m = 4
    expected = (
        (2*alpha**2 + 4*alpha + 2) 
        * monomial_symmetric_polynomial(m, IntegerPartition([3,1]))
        + (6*alpha + 10) 
        * monomial_symmetric_polynomial(m, IntegerPartition([2,1,1]))
        + (4*alpha + 4) 
        * monomial_symmetric_polynomial(m, IntegerPartition([2,2]))
        + 24 * monomial_symmetric_polynomial(m, IntegerPartition([1,1,1,1]))
    )
    obtained = JackSymbolicPol(m, mu)
    assert obtained == expected.set_domain('QQ(alpha)')