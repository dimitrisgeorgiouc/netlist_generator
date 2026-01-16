#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import math
import time
import networkx as nx

from itertools import chain, combinations, combinations_with_replacement

sys.path.append(os.path.dirname(__file__))

from _Shape import *
from _Interval import *

# Layer-By-Layer DAG Generation Algorithm 
# Exhaustive Search Based on Point Iteration
# Param：
#   - n;
#   - l;
#   - s;

class LBL_Dag_Generator:
    def __init__(self, _n:int, _lr:Interval, _sr:Interval):
        assert _n in _lr * _sr
        # 0 <= lr[0] * sr[0] <= n <= lr[1] * sr[1]
        self.n = _n
        self.lr = _lr
        self.sr = _sr

    def gen(self):
        for __l in self.lr.iter():
            for __s in SG(self.n, __l, self.sr):
                # yield __s
                for __d in self.LBLGen(__s):
                    yield __d

    def LBLGen(self, _s:tuple):
        __n, __l = sum(_s), len(_s)

        if __l == 1:
            __rd = nx.DiGraph()
            __rd.add_nodes_from(range(__n))
            yield __rd

        else:
            __as = set([sum(_s[:__si]) + __ni for __si, __sn in enumerate(_s[:-2]) for __ni in range(__sn)])
            __ps = set([sum(_s[:-2]) + __ni for __ni in range(_s[-2])])
            for __sd in self.LBLGen(_s[:-1]):
                for __png in combinations_with_replacement(self.__deff(__sd, __ps, __as), _s[-1]):
                    __rd = nx.DiGraph(__sd)
                    __rd.add_nodes_from([__ni + __n - _s[-1] for __ni in range(_s[-1])])
                    __rd.add_edges_from([(__pi, __ni + sum(_s[:-1])) for __ni, __ps in enumerate(__png) for __pi in __ps])
                    yield __rd

    def __deff(self, _sdag, _ps, _as):
        for _plist in chain.from_iterable([combinations(_ps, __pn) for __pn in range(1, len(_ps) + 1)]):
            for _alist in nx.antichains(_sdag.subgraph(_as - set.union(*(nx.ancestors(_sdag, _pi) for _pi in _plist)))):
                    yield tuple(_alist) + _plist


if __name__ == "__main__":
    # _ex_t = 'E1'
    for _n in range(3, 50):
        # Exam 1
        _lr, _sr = Interval(1, 4), Interval(1, math.ceil(_n / 2))

        # Exam 2
        # _lr, _sr = Interval(_n - 3, _n), Interval(1, math.ceil(_n / 2))

        # print(f"{_lr}_{_sr}")
        _g = LBL_Dag_Generator(_n, _lr, _sr)
        _dn = 0
        _st = time.time()
        for _d in _g.gen():
            _dn += 1
        _et = time.time()
        print(f"{_n},\t{_dn},\t{_et - _st:.6f}")
               
