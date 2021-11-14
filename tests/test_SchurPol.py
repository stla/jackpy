# -*- coding: utf-8 -*-
from sympy.combinatorics.partitions import IntegerPartition
from jackpy.jack import SchurPol
from sympy import symbols, Poly


def test_schurpol():
    n = 4
    x_0, x_1, x_2, x_3 = symbols("x_0, x_1, x_2, x_3")
    expected = (
        Poly(x_0, x_0) + Poly(x_1, x_1) + Poly(x_2, x_2) + Poly(x_3, x_3) 
    )**n
    obtained = (
        SchurPol(n, IntegerPartition([4]))
        + 3 * SchurPol(n, IntegerPartition([3,1]))
        + 2 * SchurPol(n, IntegerPartition([2,2]))
        + 3 * SchurPol(n, IntegerPartition([2,1,1]))
        + SchurPol(n, IntegerPartition([1,1,1,1]))
    )
    assert obtained == expected
