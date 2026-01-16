#!/usr/bin/python3
# -*- coding: utf-8 -*-

import time 
import math
import networkx as nx

from itertools import chain, combinations

class GNM_Dag_Generator:
    def __init__(self, n:int, mr:tuple=(0, float('Inf'))):
        assert n > 0
        self.n = n
        self.mr = (max(mr[0], 0), min(mr[1], math.comb(n, 2)))
        self.es = [(_i, _j) for _i in range(self.n) for _j in range(_i + 1, self.n)]
        self.ec = chain.from_iterable([combinations(self.es, _m) for _m in range(self.mr[0], self.mr[1] + 1)])

    def gen(self):
        for _el in self.ec:
            _ret = nx.DiGraph()
            _ret.add_nodes_from(range(self.n))
            _ret.add_edges_from(_el)
            yield _ret


if __name__ == "__main__":
    # Experment 1
    for _n in range(1, 100):
        # basic
        # _mr = (0, math.comb(_n, 2))

        # Exam 1.1
        # _mr = (0, _n)

        # Exam 1.2
        # _mr = (0, int(_n / 3))

        # Exam 1.3
        _mr = (math.ceil(_n / 3), int(2 * _n / 3))

        # Exam 1.4
        # _mr = (math.ceil(2 * _n / 3), _n)

        # print(f"{_mr}")
        _g = GNM_Dag_Generator(_n, _mr)
        _dn = 0

        _st = time.time()

        for _d in _g.gen():
            _dn += 1

        _et = time.time()

        # print(f"n:{_n}\n\tdn:{_dn}\n\time:{_et - _st:.6f}")
        print(f"{_n},\t{_dn},\t{_et - _st:.12f}")
