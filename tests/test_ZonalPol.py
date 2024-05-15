# -*- coding: utf-8 -*-
from jackpy.jack import ZonalPol
from sympy import symbols, Poly


def test_zonalpol():
    n = 4
    v1, v2, v3, v4 = symbols("x_1, x_2, x_3, x_4")
    x_1 = Poly(v1, v1, domain='QQ')
    x_2 = Poly(v2, v2, domain='QQ')
    x_3 = Poly(v3, v3, domain='QQ')
    x_4 = Poly(v4, v4, domain='QQ')
    expected = (x_1 + x_2 + x_3 + x_4)**(n - 1)
    obtained = (
        ZonalPol(n, [3])
        + ZonalPol(n, [2,1])
        + ZonalPol(n, [1,1,1])
    )
    assert obtained == expected

