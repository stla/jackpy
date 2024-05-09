# -*- coding: utf-8 -*-
from sympy.combinatorics.partitions import IntegerPartition
from jackpy.jack import SchurPol
from sympy import symbols, Poly


def test_schurpol():
    n = 4
    v1, v2, v3, v4 = symbols("x_1, x_2, x_3, x_4")
    x_1 = Poly(v1, v1, domain='ZZ')
    x_2 = Poly(v2, v2, domain='ZZ')
    x_3 = Poly(v3, v3, domain='ZZ')
    x_4 = Poly(v4, v4, domain='ZZ')
    expected = (x_1 + x_2 + x_3 + x_4)**n
    obtained = (
        SchurPol(n, IntegerPartition([4]))
        + 3 * SchurPol(n, IntegerPartition([3,1]))
        + 2 * SchurPol(n, IntegerPartition([2,2]))
        + 3 * SchurPol(n, IntegerPartition([2,1,1]))
        + SchurPol(n, IntegerPartition([1,1,1,1]))
    )
    assert obtained == expected
