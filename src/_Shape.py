#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys

sys.path.append(os.path.dirname(__file__))

from _Interval import *

def SG(_n:int, _l:int, _sr:Interval=Interval(1, float('Inf'))):
    __nsr =  Interval(1, _n - (_l - 1)) & _sr 
    if 1 <= _l <= _n and _n in __nsr * _l:
        # CASE 1 NULL DAG
        if _l == 1:
            yield (_n, )

        # CASE 2 CHAIN DAG
        elif _l == _n:
            yield tuple([1 for _ in range(_l)])

        # CASE 3
        elif __nsr.lower == __nsr.upper:
            yield tuple([set(__nsr.lower) for _ in range(_l)])

        # CASE * OTHER DAG
        else:
            for __Xp, __Xl in map(lambda __x: (_n - __x, __x), __nsr.iter()):
                # reserved DB ；
                for _ps in SG(__Xp, _l - 1, __nsr):
                    yield _ps + (__Xl, )
