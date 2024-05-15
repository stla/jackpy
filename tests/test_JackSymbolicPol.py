# -*- coding: utf-8 -*-
from sympy import symbols
from jackpy.jack import JackSymbolicPol
from jackpy.monomial_symmetric_polynomials import (
        monomial_symmetric_polynomial
    ,   msp_combination_expr
    )

def test_jacksymbolicpol():
    mu = [3, 1]
    alpha = symbols("alpha")
    m = 4
    expected = (
        (2*alpha**2 + 4*alpha + 2) 
        * monomial_symmetric_polynomial(m, [3,1])
        + (6*alpha + 10) 
        * monomial_symmetric_polynomial(m, [2,1,1])
        + (4*alpha + 4) 
        * monomial_symmetric_polynomial(m, [2,2])
        + 24 
        * monomial_symmetric_polynomial(m, [1,1,1,1])
    )
    obtained = JackSymbolicPol(m, mu)
    assert obtained == expected.set_domain('QQ(alpha)')

def test_symbolic_msp_combination_expr():
    mu = [3, 1]
    alpha = symbols("alpha")
    m = 4
    jp = JackSymbolicPol(m, mu)
    expr = msp_combination_expr(jp)
    expected = (
        (2*alpha**2 + 4*alpha + 2) 
        * symbols("M[3;1]", commutative=False)
        + (6*alpha + 10) 
        * symbols("M[2;1;1]", commutative=False)
        + (4*alpha + 4) 
        * symbols("M[2;2]", commutative=False)
        + 24 
        * symbols("M[1;1;1;1]", commutative=False)
    )
    assert expr == expected
