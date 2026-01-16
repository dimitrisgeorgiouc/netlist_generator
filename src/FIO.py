#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import time 
import networkx as nx

sys.path.append(os.path.dirname(__file__))

from _Interval import *


class FIO_Dag_Generator:
    def __init__(self, _n:int, _id:int=float('Inf'), _od:int=float('Inf'), _lr:Interval=Interval(1, float('Inf'))):
        assert 0 < _n 
        self.n = _n
        self.id = min(_id, _n - 1)
        self.od = min(_od, _n - 1)
        self.lr = _lr & Interval(1, _n)

    def gen(self):
        yield from self.FIOGen(self.n, self.lr)

    def FIOGen(self, _n:int, _lr:Interval):
        if 1 in _lr:
            __rd = nx.DiGraph()
            __rd.add_nodes_from([(__i, {'d': 1}) for __i in range(_n)])
            yield __rd
        if 2 <= _lr.upper:
            for __sd in self.FIOGen(_n - 1, (_lr - Interval(1, 0)) & Interval(1, _n)):
                __rni = __sd.number_of_nodes()
                for __pns in nx.antichains(__sd):
                    if 0 <= len(__pns) <= self.id and \
                       max([__sd.out_degree(__pnx) for __pnx in __pns], default=0) < self.od:
                        __rd = nx.DiGraph(__sd)
                        __rd.add_nodes_from([(__rni, {'d': max([__rd.nodes[__pnx]['d'] for __pnx in __pns], default=0) + 1})])
                        __rd.add_edges_from([(__pnx, __rni) for __pnx in __pns])
                        if max([__snd['d'] for __sni, __snd in __rd.nodes(data=True)]) in _lr:
                            yield __rd


if __name__ == "__main__":
    for _n in range(1, 100):
        # print('new')
        _st = time.time()
        # Exam 1
        # _id, _od, _lr = 1, 2, Interval(_n - 4, _n)
        # Exam 2
        # _id, _od, _lr = 1, 2, Interval(1, 5)
        # Exam 3
        # _id, _od, _lr = 2, 1, Interval(_n - 4, _n)
        # Exam 4
        _id, _od, _lr = 2, 1, Interval(1, 5)
        # print(f"{_id}_{_od}_{_lr}")
        _g = FIO_Dag_Generator(_n, _id=_id, _od=_od, _lr=_lr)
        _dn = 0
        _st = time.time()

        for _d in _g.gen():
            # print(d.edges())
            _dn += 1

        _et = time.time()
               
        # print(f"n:{_n}\n\tdn:{_dn}\n\time:{_et - _st:.6f}")
        print(f"{_n},\t{_dn},\t{_et - _st:.6f}")