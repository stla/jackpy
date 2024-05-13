# -*- coding: utf-8 -*-
from gmpy2 import mpq
from sympy.combinatorics.partitions import IntegerPartition
from sympy import symbols, Poly
from jackpy.jack import JackPol, SchurPol
from jackpy.monomial_symmetric import (
        monomial_symmetric_polynomial
    ,   msp_combination
    )

def check_mscombo(poly):
    combo = msp_combination(poly)
    n = len(poly.gens)
    variables = [symbols(f'x_{i}') for i in range(1, n+1)]
    kappas = combo.keys()
    q = Poly(0, *variables)
    for kappa in kappas:
        part = IntegerPartition(kappa)
        q = q + combo.get(kappa)*monomial_symmetric_polynomial(n, part)
    return q.as_dict()

def test_msp_combination():
    jp = JackPol(6, IntegerPartition([3,2,1]), mpq(2, 7))
    dic = check_mscombo(jp)
    assert dic == jp.as_dict()

def test_jackpol():
    mu = IntegerPartition([3,1])
    alpha = mpq(5, 2)
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
    obtained = JackPol(m, mu, alpha)
    assert obtained == expected

def test_jack_p_pol_is_schurpol():
    mu = IntegerPartition([3,2,1])
    n = 4
    expected = SchurPol(n, mu)
    obtained = JackPol(n, mu, mpq(1), which = 'P')
    assert obtained == expected.set_domain('QQ')