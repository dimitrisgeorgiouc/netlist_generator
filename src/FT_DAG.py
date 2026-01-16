#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import time
import math
import pickle
import hashlib
import sqlite3
import numpy as np
import igraph as ig
import networkx as nx



from collections import Counter
from itertools import chain, product, zip_longest, combinations_with_replacement

sys.path.append(os.path.dirname(__file__))

from _Shape import *
from _Interval import *

global GTableHead, DTableHead, TCONN, TCURS

# L : 0 writing；
# L : 1 excutable；
# L : 2 null data；

# Table Name: GTable
GTableHead = {'ts': ('TEXT', 'PRIMARY KEY'), 'L': ('TINYINT', 'NOT NULL'),
              'id': ('TINYINT UNSIGNED', ''), 'od': ('TINYINT UNSIGNED', ''),
              'jld': ('TINYINT UNSIGNED', ''), 'jlh': ('TINYINT UNSIGNED', ''),
              'w_l': ('TINYINT UNSIGNED', ''), 'w_u': ('TINYINT UNSIGNED', ''),
              'm_l': ('INTEGER UNSIGNED', ''), 'm_u': ('INTEGER UNSIGNED', ''),
              'sl': ('TINYINT UNSIGNED', ''), }

# Table Name: ts (c BLOB PRIMARY KEY)
DTableHead = {'c': ('TEXT', 'PRIMARY KEY'), 'dag': ('BLOB', 'NOT NULL'),
              'm': ('INTEGER UNSIGNED', 'NOT NULL'), 'w': ('TINYINT UNSIGNED', 'NOT NULL'),
              'id': ('TINYINT UNSIGNED', 'NOT NULL'), 'od': ('TINYINT UNSIGNED', 'NOT NULL'),
              'jld': ('TINYINT UNSIGNED', 'NOT NULL'), 'jlh': ('TINYINT UNSIGNED', 'NOT NULL'),
              'sl': ('TINYINT UNSIGNED', 'NOT NULL'), }


class FT_Dag_Generator:
    def __init__(self, _s: tuple, _sl: int = 1,
                 _mr: Interval = Interval(0, float('Inf')),
                 _wr: Interval = Interval(1, float('Inf')),
                 _id: int = float('Inf'), _od: int = float('Inf'),
                 _jld: int = float('Inf'), _jlh: int = float('Inf')):

        self.s = _s
        self.mr = _mr
        self.wr = _wr
        self.id = _id
        self.od = _od
        self.sl = _sl
        self.jld = _jld
        self.jlh = _jlh

    def gen(self):
        for _d in self.DG(self.s, self.mr, self.wr, self.id, self.od, self.jld, self.jlh, self.sl, _ms=True, _mc=True):
            yield _d

    #################################
    def __slength_constrin(self, _s: tuple, _sl: int = 1) -> int:
        return max(1, min(_sl, len(_s)))

    def __in_degree_constrin(self, _s: tuple, _id: int = float('Inf')) -> int:
        if len(_s) > 1:
            return max(1, min(_id, sum(_s[:-1]) - len(_s) + 2))
        else:
            return 0

    def __out_degree_constrin(self, _s: tuple, _od: int = float('Inf')) -> int:
        if len(_s) > 1:
            return max(1, min(_od, sum(_s[1:]) - len(_s) + 2))
        else:
            return 0

    def __jump_layer_constrin(self, _s: tuple, _jld: int = float('Inf')) -> int:
        if len(_s) > 1:
            return max(1, min(_jld, len(_s) - 1))
        else:
            return 0

    def __jump_level_constrin(self, _s: tuple, _jlh: int = float('Inf')) -> int:
        if len(_s) > 1:
            return max(1, min(_jlh, len(_s) - 1))
        else:
            return 0

    def __edge_constrin(self, _s: tuple, _mr: Interval = Interval(0, float('Inf'))) -> Interval:
        return Interval(sum(_s[1:]), int(pow(sum(_s), 2) / 4)) & _mr

    def __width_constrin(self, _s: tuple, _wr: Interval = Interval(0, float('Inf'))) -> Interval:
        return Interval(max(_s), sum(_s) - len(_s) + 1) & _wr

    #################################
    def DagWidth(self, _cdag, _ans: set = set()):
        """ 默认使用霍普克洛夫特-卡普算法计算最大匹配 """
        _ons = _cdag.nodes() - _ans
        __dn, __tg, __us, __te = len(_ons), nx.DiGraph(), list(), list()
        for __pi, __o_d in map(lambda __x: (f"p{__x}", _cdag.nodes[__x]['D'] - _ans), _ons):
            if len(__o_d) > 0:
                __us.append(__pi)
                for __si in map(lambda __y: f"s{__y}", __o_d):
                    __te.append((__pi, __si))
        __tg = nx.DiGraph(__te)
        return __dn - int(len(nx.algorithms.bipartite.matching.hopcroft_karp_matching(__tg, __us)) / 2)

    def __DB_Data(self, _hash_s, _s: tuple, _mr: Interval, _wr: Interval, _id: int, _od: int, _jld: int, _jlh: int,
                  _sl: int, _dt: str):
        global DTableHead, TCURS
        __GDGen, __tname = None, f"{_dt}_{_hash_s}"
        if _dt == 'SCS':
            __GDGen = self.SCSGen(_s, _mr, _wr, _id, _od, _jld, _jlh, _sl)
        elif _dt == 'MC':
            __GDGen = self.MCGen(_s, _mr, _wr, _id, _od, _jld, _jlh, _sl)
        elif _dt == 'MS':
            __GDGen = self.MSGen(_s, _mr, _wr, _id, _od, _jld, _jlh, _sl)
        else:
            assert False

        TCURS.execute(f'SELECT * FROM GTable WHERE ts = "{__tname}"')
        __db_data = TCURS.fetchall()

        # C1 Data Existe
        if len(__db_data) > 0:
            __red = dict(__db_data[0])

            # C1.1 Data Cover；
            if _wr in Interval(__red['w_u'], __red['w_l']) and _mr in Interval(__red['m_l'], __red['m_u']) and \
                    __red['id'] >= _id and __red['od'] >= _od and __red['jld'] >= _jld and __red['jlh'] >= _jlh and \
                    __red['sl'] <= _sl and __red['L'] == 1:
                TCURS.execute(f'''SELECT dag FROM {__tname} WHERE 
                                    sl >= {_sl} AND od <= {_od} AND id <= {_id} AND jld <= {_jld} AND jlh <= {_jlh} AND 
                                    {_mr.upper} >= m AND m >= {_mr.lower} AND {_wr.upper} >= w AND w >= {_wr.lower} ''')
                for __DbDag in TCURS.fetchall():
                    yield pickle.loads(__DbDag[0])

            # C1.2 Other，By Algorithm (Database-independent)；
            else:
                for __rd in __GDGen:
                    yield __rd

        # C2 Data Not Existe
        else:
            # New gtable line and data table（locking L=0）；
            TCURS.execute(
                f"INSERT INTO GTable (ts, m_l, m_u, w_l, w_u, id, od, jld, jlh, sl, L) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", \
                (__tname, _mr.lower, _mr.upper, _wr.lower, _wr.upper, _id, _od, _jld, _jlh, _sl, 0))
            __data_buff = []
            for __rd in __GDGen:
                yield __rd
                __data_buff.append({"c": __rd.graph['c'], "dag": pickle.dumps(__rd),
                                    "m": __rd.graph['m'], "w": __rd.graph['w'],
                                    "id": __rd.graph['id'], "od": __rd.graph['od'],
                                    "jld": __rd.graph['jld'], "jlh": __rd.graph['jlh'],
                                    "sl": __rd.graph['sl'], })

            # C2.1 Data Input (release L = 1)；
            if len(__data_buff) > 0:
                TCURS.execute(
                    f'''CREATE TABLE {__tname} ({', '.join([f'{_i} {_t[0]} {_t[1]}' for _i, _t in DTableHead.items()])});''')

                TCURS.executemany(
                    f''' INSERT INTO {__tname} (c, m, w, id, od, jld, jlh,  sl, dag) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) ''',
                    [(__dx['c'], __dx['m'], __dx['w'], __dx['id'], __dx['od'], __dx['jld'], __dx['jlh'], __dx['sl'],
                      __dx['dag'],) for __dx in __data_buff])

                TCURS.execute(f'UPDATE GTable  SET L = 1 WHERE ts = "{__tname}"')

            # C2.2 Null Data (L = 2)；
            else:
                TCURS.execute(f'UPDATE GTable  SET L = 2 WHERE ts = "{__tname}"')
        # for __rd in __GDGen:
        #     yield __rd
        __GDGen.close()

    """ (*) Main """

    def DG(self, _s: tuple, _mr: Interval, _wr: Interval, _id: int, _od: int, _jld: int, _jlh: int, _sl: int,
           _ms: bool = True, _mc: bool = True):
        __n, __l = sum(_s), len(_s)
        __nmr = self.__edge_constrin(_s, _mr)
        __nwr = self.__width_constrin(_s, _wr)
        __nsl = self.__slength_constrin(_s, _sl)
        __nid = self.__in_degree_constrin(_s, _id)
        __nod = self.__out_degree_constrin(_s, _od)
        __njld = self.__jump_layer_constrin(_s, _jld)
        __njlh = self.__jump_level_constrin(_s, _jlh)

        if all([_s[__i - 1] * __nod >= _s[__i] for __i in range(1, len(_s))]):
            __hash_s = f"s{hashlib.md5(pickle.dumps(_s)).hexdigest()}"

            """ (1) SCS DAG """
            if (__n == 1) or (__n > 3 and _s[0] > 1 and __l > 1 and __nid > 1):
                for __rd in self.__DB_Data(__hash_s, _s, __nmr, __nwr, __nid, __nod, __njld, __njlh, __nsl, 'SCS'):
                    yield __rd

            """ (2) MC DAG """
            if _mc and _s[0] > 1:
                for __rd in self.__DB_Data(__hash_s, _s, __nmr, __nwr, __nid, __nod, __njld, __njlh, __nsl, 'MC'):
                    yield __rd

            """ (3) MS DAG """
            if _ms and __l > 1:
                for __rd in self.__DB_Data(__hash_s, _s, __nmr, __nwr, __nid, __nod, __njld, __njlh, __nsl, 'MS'):
                    yield __rd

    """ (1) SCS-DAG """

    def SCSGen(self, _s: tuple, _mr: Interval, _wr: Interval, _id: int, _od: int, _jld: int, _jlh: int, _sl: int):
        global Tim
        __n, __l = sum(_s), len(_s)

        # (1) Trivial DAG
        if _sl == __n == 1 and 0 in _mr and 1 in _wr:
            yield self.NullDAG(__n)

        # (2) Other DAG (Expect T-DAG)；l & s[0] > 1；
        if __n > 3 and __l > 1 and _s[0] > 1 and _sl > 0 and _id > 0 and _od > 0 and _jld > 0 and _jlh > 0:
            __ss = _s[:-1] if _s[-1] == 1 else _s[:-1] + (_s[-1] - 1,)
            __sn, __cbuff = sum(__ss), set()
            __smr = _mr - Interval(_id, 1)
            __swr = _wr - Interval(1, 0) if _s[-1] > 1 else _wr

            # S1. Generate sub-DAGs that satisfies the condition
            for __sd in self.DG(__ss, __smr, __swr, _id, _od, _jld, _jlh, min(_sl, __l - _jld), _ms=False,
                                _mc=True if _id > 1 else False):

                # S2. Constraint Condition Determination and Connected Component Grouping of Sub-DAG:
                if __cng := self.CCD(__sd, __l, _id, _od, _jld, _jlh, _sl):

                    __sd_sink = set([__ni for __ni, __nd in __sd.nodes(data=True) if len(__nd['S']) == 0])

                    __ss_cn = sum([max(1, len(__mns)) for __mns, _ in __cng])

                    __ss_w, __ss_en = __sd.graph['w'], __sd.number_of_edges()

                    # S3. Based on the sub-DAG that meets the conditions, further exhaust all predecessor combinations.
                    # * The number of necessary connection points：
                    for __pns in self.NPlGen(__cng, __l - 1, Interval(__ss_cn, _id) & (_mr - __ss_en)):

                        # S3_1. sink nodes det
                        __pni_set = set(__pnsx[0] for __pnsx in __pns)
                        if __pni_set != __sd_sink:

                            # S3_2. width det
                            __pi_ance = __pni_set | set().union(*[__pid[1]['A'] for __pid in __pns])
                            __scs_w = max(__ss_w, self.DagWidth(__sd, __pi_ance) + 1) if _s[-1] > 1 else __ss_w

                            if __scs_w in _wr:
                                __n_es = tuple([(__epi, __esi) for __epi, __esi in __sd.edges()]) + tuple(
                                    [(__pi, __sn) for __pi in __pni_set])
                                __n_am = np.zeros((__n, __n), dtype=int)
                                __n_am[tuple(zip(*__n_es))] = 1
                                __cx = hashlib.md5(self.certi_wadjm(__n_am).tobytes()).hexdigest()

                                # S5. Data Return
                                if __cx not in __cbuff:
                                    __cbuff.add(__cx)
                                    __scs = nx.DiGraph()
                                    __n_ns = tuple([(__sni, {'d': __snd['d'], 'h': __snd['h'], 'Equ': __snd['Equ'],
                                                             'P': set() | __snd['P'], 'S': set() | __snd['S'],
                                                             'A': set() | __snd['A'], 'D': set() | __snd['D']}
                                                     ) for __sni, __snd in __sd.nodes(data=True)])
                                    __n_ns += tuple([(__sn, {'d': __l, 'h': 1, 'Equ': None,
                                                             'P': set() | __pni_set, 'S': set(),
                                                             'A': set() | __pi_ance, 'D': set()}), ])

                                    for __n_ni, __n_nd in __n_ns:
                                        if __n_ni in __pi_ance:
                                            __n_nd['D'].add(__sn)
                                            if __n_ni in __pni_set:
                                                __n_nd['S'].add(__sn)

                                    __scs.add_nodes_from(__n_ns)
                                    __scs.add_edges_from(__n_es)

                                    self.__height_update(__scs, __sn)

                                    __rd_equ = list()
                                    for __ni in range(__n):
                                        __nel = (frozenset(__scs.nodes[__ni]['P']), frozenset(__scs.nodes[__ni]['S']))
                                        if __nel not in __rd_equ:
                                            __rd_equ.append(__nel)
                                        __scs.nodes[__ni]['Equ'] = __rd_equ.index(__nel)

                                    __rlist = []
                                    for __souri in __scs.nodes():
                                        if __scs.nodes[__souri]['d'] == 1:
                                            for __sinki in __scs.nodes():
                                                if __scs.nodes[__sinki]['h'] == 1:
                                                    if nx.has_path(__scs, source=__souri, target=__sinki):
                                                        __rlist.append(nx.shortest_path_length(__scs, source=__souri,
                                                                                               target=__sinki))

                                    __scs.graph = {'c': __cx, 's': _s, 'w': __scs_w, 'dt': 'SCS', 'm': len(__n_es),
                                                   'id': max([len(__rd['P']) for _, __rd in __scs.nodes(data=True)]),
                                                   'od': max([len(__rd['S']) for _, __rd in __scs.nodes(data=True)]),
                                                   # 'sl': min([__rd['h'] + __rd['d'] for _, __rd in __scs.nodes(data=True)]),
                                                   'sl': min(__rlist) + 1,
                                                   'jld': max([__scs.nodes[__esi]['d'] - __scs.nodes[__epi]['d'] for
                                                               __epi, __esi in __scs.edges()]),
                                                   'jlh': max([__scs.nodes[__epi]['h'] - __scs.nodes[__esi]['h'] for
                                                               __epi, __esi in __scs.edges()])}
                                    # DAG_Det(__scs)
                                    yield __scs

    def CCD(self, _tdag, _ln: int, _id: int, _od: int, _jld: int, _jlh: int, _sl: int):
        if len([_i for _i, _d in _tdag.nodes(data=True) if _d['d'] == _ln - 1 and len(_d['S']) < _od]) > 0:
            __min_pn, __rdg = 0, tuple()
            for __scs_cx in nx.weakly_connected_components(_tdag):
                __cx_d, __cx_ms, __cx_ns, __cx_md = dict(), set(), set(), set()
                for __ni, __nd in map(lambda __x: (__x, _tdag.nodes[__x]), __scs_cx):
                    if len(__nd['S']) > _od and len(__nd['P']) > _id:
                        # __nd['d'] < min(_sl, _ln - _jld)
                        # __nd['h'] < min(_sl, _ln - _jlh)
                        return False
                    elif __nd['d'] <= _ln - 1:
                        __cx_d[__ni] = {'d': __nd['d'], 'h': __nd['h'], 'Equ': __nd['Equ'],
                                        'P': set() | __nd['P'], 'S': set() | __nd['S'],
                                        'A': set() | __nd['A'], 'D': set() | __nd['D']}
                        if len(__nd['S']) == 0 and __nd['d'] < _sl:
                            __cx_ms.add(__ni)
                            __cx_md |= __nd['D'] | __nd['A']

                        elif len(__nd['S']) < _od and _ln - _jld <= __nd['d'] and __nd['h'] - 1 <= _jlh:
                            __cx_ns.add(__ni)

                __min_pn += max(1, len(__cx_ms))
                if len(__cx_d) > 0 and __min_pn <= _id:
                    __rdg += ((tuple([(__mni, __cx_d[__mni]) for __mni in __cx_ms]),
                               tuple([(__nni, __cx_d[__nni]) for __nni in __cx_ns - __cx_md])),)
                else:
                    return False
            return __rdg
        return False

    def NPlGen(self, _pdg, _pln: int, _pnr: Interval, _pn: bool = False):
        if 1 <= _pnr.lower <= _pnr.upper:
            __smns, __snns = _pdg[0]
            __rpn_min = sum([max(1, len(__y)) for __y, _ in _pdg[1:]])
            __rpn_max = sum([len(__y) + len(__z) for __y, __z in _pdg[1:]])
            __tpn_max = _pnr.upper - __rpn_min - len(__smns)
            __tpn_min = max(1, _pnr.lower - __rpn_max - len(__smns))
            __anti_gen = self.NxAnti(sorted(__snns, key=lambda x: x[0]), __tpn_min, __tpn_max)

            if len(__smns) > 0:
                __anti_gen = chain(__anti_gen, [tuple()])
            for __rpnx in map(lambda __x: __x + __smns, __anti_gen):
                __rpnn, __pnlb = len(__rpnx), (_pln in set([__rn['d'] for __ri, __rn in __rpnx])) or _pn
                if len(_pdg) == 1:
                    if __pnlb:
                        yield __rpnx
                else:
                    for __spn in self.NPlGen(_pdg[1:], _pln, Interval(__rpn_min, __rpn_max) & (_pnr - __rpnn), __pnlb):
                        yield __spn + __rpnx

    def __height_update(self, _td, _ni):
        for __pi in _td.nodes[_ni]['P']:
            __temp_h = max([_td.nodes[__si]['h'] for __si in _td.nodes[__pi]['S']]) + 1
            assert _td.nodes[__pi]['h'] <= __temp_h
            if _td.nodes[__pi]['h'] < __temp_h:
                _td.nodes[__pi]['h'] = __temp_h
                self.__height_update(_td, __pi)

    def NxAnti(self, _ns: tuple, _nd: int, _nu: int):
        if 1 <= _nd <= _nu:
            __equ_buff = set()
            for __nid in _ns:
                if __nid[1]['Equ'] not in __equ_buff:
                    __equ_buff.add(__nid[1]['Equ'])
                    if _nd == 1:
                        yield (__nid,)
                    if _nu > 1:
                        __rns = tuple([__sid for __sid in _ns if
                                       __sid[0] > __nid[0] and __sid[0] not in set(__nid[1]['A'] | __nid[1]['D'])])
                        for __s_an in self.NxAnti(__rns, max(1, _nd - 1), _nu - 1):
                            yield (__nid,) + __s_an

    """ (2) MC-DAG """

    def MCGen(self, _s: tuple, _mr: Interval, _wr: Interval, _id: int, _od: int, _jlh: int, _jld: int, _sl: int):
        __n, __l = sum(_s), len(_s)
        if _sl == __l == 1 and 0 in _mr and __n in _wr:
            yield self.NullDAG(__n)

        if 1 < __l and 1 <= _id and 1 <= _od and 1 <= _jld:
            for __lgn in range(_s[0] if _sl > 1 else 1, _s[0] + 1):
                __ogn = _s[0] - __lgn
                for __tlg_d in self.number_combine(__lgn, (1 if _s[0] > __lgn else 2, min(_s[1], __lgn, min(_s[:_sl]))),
                                                   (1, __lgn)):
                    __tlg = tuple(((__gd,), __gn) for __gd, __gn in Counter(__tlg_d).items())
                    for __lg in self.LG(__tlg, _s[1:], _od, _sl):
                        __mc_sg = tuple(_x for _x, _i in __lg for _ in range(_i))
                        __mc_mr, __mc_wr = Interval(0, 0), Interval(0, 0)
                        for __mcs in __mc_sg:
                            __mc_mr += self.__edge_constrin(__mcs)
                            __mc_wr += self.__width_constrin(__mcs)
                        for __rdag_mc in self.RG(tuple(Counter(__mc_sg).items()), __mc_mr & _mr,
                                                 __mc_wr & (_wr - __ogn), _id, _od, _jld, _jlh, _sl):
                            yield self.EM(__rdag_mc + (self.NullDAG(__ogn),)) if __ogn > 0 else self.EM(__rdag_mc)

    def number_combine(self, _n: int, _lr: tuple, _sr: tuple):
        if 1 <= _lr[0] <= _lr[1] and 1 <= _sr[0] <= _sr[1] and _sr[0] * _lr[0] <= _n <= _sr[1] * _lr[1]:
            for __gxn in range(max(1, _sr[0]), min(_sr[1], _n - (_lr[0] - 1) * _sr[0]) + 1):
                __grn = _n - __gxn
                if __grn == 0:
                    yield (__gxn,)
                else:
                    for __sub_comb in self.number_combine(__grn, (max(_lr[0] - 1, 1), _lr[1] - 1), (__gxn, _sr[1])):
                        yield __sub_comb + (__gxn,)

    def LS(self, _sg, _sn):
        # lo grouping
        __tsd, __tsn = _sg[0]
        for __tlsn in range(min(_sn, __tsn) + 1):
            __rlsg = ((__tsd, __tlsn),) if __tlsn > 0 else tuple()
            __rnsg = ((__tsd, __tsn - __tlsn),) if __tsn > __tlsn else tuple()
            if __tlsn == _sn or len(_sg) == 1:
                yield __rlsg, __rnsg + _sg[1:]
            else:
                for __srlsg, __srnsg in self.LS(_sg[1:], _sn - __tlsn):
                    yield __rlsg + __srlsg, __rnsg + __srnsg

    def NA(self, _sg, _an: int, _od: int):
        __min_an = sum(__n for __s, __n in _sg)
        if __min_an == _an:
            yield tuple((__sx + (1,), __sn) for __sx, __sn in _sg)
        elif __min_an < _an:
            __tsd, __tsn = _sg[0]
            if len(_sg) == 1:
                for __lg_data in self.number_combine(_an, (__tsn, __tsn), (1, min(_an, _od * __tsd[-1]))):
                    yield tuple((__tsd + (__nn,), __sn) for __nn, __sn in Counter(__lg_data).items())
            else:
                for __t_an in range(__tsn, _an - (__min_an - __tsn) + 1):
                    for __lg_data in self.number_combine(__t_an, (__tsn, __tsn), (1, min(__t_an, _od * __tsd[-1]))):
                        __t_sg = tuple((__tsd + (__nn,), __sn) for __nn, __sn in Counter(__lg_data).items())
                        for _new_sg in self.NA(_sg[1:], _an - __t_an, _od):
                            yield _new_sg + __t_sg
        else:
            assert False

    def LG(self, _psg, _s: tuple, _od: int, _sl: int):
        # 将_s0的结点加入_psg中，不再重新分组合并 IO都是分组格式的（shape, num）
        for __newsg in self.NA(_psg, _s[0], _od):  # (1) s0结点分配,每个组必须有一个后继;
            if len(_s) == 1:  # (2) s数据处理完成直接返回；
                yield __newsg
            else:  # (3) s数据还有剩余，继续递归；
                # __tl =                                # 当前shape的长度；
                # (3.1) __newsg的L/O分组,返回__lsg, __nsg; lsg的范围[1, min(len(__newsg), _s[1])]，剩下的给__nsg，可以为空；
                if len(__newsg[0][0]) < _sl:
                    for __tu_sg in self.LG(__newsg, _s[1:], _od, _sl):
                        yield __tu_sg
                else:
                    for __lsg, __nsg in self.LS(__newsg, min(sum(__nsn for __nss, __nsn in __newsg), _s[1])):
                        if len(__lsg) > 0:
                            if sum(__lsn for __lss, __lsn in __lsg) == 1:
                                yield ((__lsg[0][0] + _s[1:], 1),) + __nsg
                            else:
                                for __sub_sg in self.LG(__lsg, _s[1:], _od, _sl):
                                    yield __sub_sg + __nsg

    def RG(self, _sg: tuple, _mr: Interval, _wr: Interval, _id: int, _od: int, _jld: int, _jlh: int, _sl: int):
        __sgmr, __sgwr = Interval(0, 0), Interval(0, 0)
        for __sgs, __sgn in _sg:
            __sgmr += self.__edge_constrin(__sgs) * __sgn
            __sgwr += self.__width_constrin(__sgs) * __sgn

        __n_mr = _mr & __sgmr
        __n_wr = _wr & __sgwr
        if __n_mr.lower <= __n_mr.upper and __n_wr.lower <= __n_wr.upper:
            __lgs, __lgn = _sg[0]
            __rgmr = __n_mr - self.__edge_constrin(__lgs) * __lgn
            __rgwr = __n_wr - self.__width_constrin(__lgs) * __lgn

            __TncGen = self.DG(__lgs, Interval(int((_mr.lower - __rgmr.upper) / __lgn), _mr.upper - __rgmr.lower),
                               Interval(int((_wr.lower - __rgwr.upper) / __lgn), _wr.upper - __rgwr.lower),
                               _id, _od, _jld, _jlh, _sl, _ms=True, _mc=False)

            for __ldagg in combinations_with_replacement(__TncGen, __lgn):
                __tmn, __twn = sum(__ldx.number_of_edges() for __ldx in __ldagg), sum(
                    __ldx.graph['w'] for __ldx in __ldagg)

                if len(_sg) > 1:
                    for __odag in self.RG(_sg[1:], _mr - __tmn, _wr - __twn, _id, _od, _jld, _jlh, _sl):
                        yield __ldagg + __odag
                else:
                    if __tmn in _mr and __twn in _wr:
                        yield __ldagg

    """ (3) MS-DAG ：主生成部分可分为: p-segement & l-segment """

    def MSGen(self, _s: tuple, _mr: Interval, _wr: Interval, _id: int, _od: int, _jld: int, _jlh: int, _sl: int):
        __n, __l = sum(_s), len(_s)

        # (1) Chain DAG;
        if _sl <= __n == __l and 1 in _wr and __n - 1 in _mr:
            yield self.ChainDAG(__n)

        # (2) Other Case;
        else:
            for __lx in range(1, __l):  # ls层数：从小到大；
                __p_s, __s_s = _s[:__lx], _s[__lx:]
                __p_l, __s_l = len(__p_s), len(__s_s)

                # (1) pl 不能是MS-DAG 根图一定是MS DAG: S[0] == 1 AND L > 1 逆否命题：S[0] > 1 or L == 1
                if (1 < __p_s[0] or __p_l == 1) and (__s_s[0] <= _od) and (__p_s[-1] <= _id):
                    __s_mr = self.__edge_constrin(__s_s)
                    __s_wr = self.__width_constrin(__s_s)
                    __p_wr = self.__width_constrin(__p_s)
                    __p_wr = Interval(
                        __p_wr.lower if _wr.lower <= __s_wr.upper and __s_wr.lower <= _wr.upper else max(_wr.lower,
                                                                                                         __p_wr.lower),
                        _wr.upper)
                    __p_mr = Interval(_mr.lower - __s_mr.upper - __p_wr.upper * __s_s[0],
                                      _mr.upper - __s_mr.lower - __p_s[-1] * __s_s[0])
                    # (2.1) P-Segment (pd are not MS-DAG)：lg的w可以达则上限不超标即可，否则pg必须达标；
                    for __ps_dg in self.DG(__p_s, __p_mr, __p_wr, _id, _od, _jld, _jlh, __p_l - _jld + 1, _ms=False,
                                           _mc=True if _id > 1 else False):
                        __psd_sn = len({__i for __i, __d in __ps_dg.nodes(data=True) if len(__d['S']) == 0})
                        if __psd_sn <= _id:
                            __nm_num = __psd_sn * __s_s[0] + __ps_dg.number_of_edges()
                            __s_mr = _mr - __nm_num
                            if not __ps_dg.graph['w'] in _wr:
                                __s_wr.lower = max(_wr.lower, __s_wr.lower)
                            __s_wr.upper = min(_wr.upper, __s_wr.upper)
                            # (2.2) S-segment: pg达标则上限不超标即可，否则lg必须达标；
                            # width of pd sat in constrain, otherwise(pg sat up constrain) width of sd must sat in constrain.
                            for __ls_dg in self.DG(__s_s, __s_mr, __s_wr, _id, _od, _jld, _jld, max(1, _sl - __s_l),
                                                   _ms=True, _mc=True):  # l-segment
                                __rpl = self.FC(__ps_dg, __ls_dg)
                                if __rpl.graph['w'] in _wr:
                                    yield __rpl

    """ """

    def EM(self, _dg):
        __rd = nx.DiGraph()
        __rd_equ = [None, ]
        __rs, __rc = list(), list()
        __rd_ns, __rd_es = list(), list()

        for __ti, __td in enumerate(_dg):
            __rd_n = len(__rd_ns)

            for __i, __d in __td.nodes(data=True):
                __nel = None if __td.number_of_nodes() == 1 else f"{__ti}_{__d['Equ']}"
                if __nel not in __rd_equ:
                    __rd_equ.append(__nel)

                __rd_ns.append((__i + __rd_n,
                                {'d': __d['d'], 'h': __d['h'], 'Equ': __rd_equ.index(__nel),
                                 'P': set(__p_i + __rd_n for __p_i in __d['P']),
                                 'S': set(__s_i + __rd_n for __s_i in __d['S']),
                                 'A': set(__a_i + __rd_n for __a_i in __d['A']),
                                 'D': set(__d_i + __rd_n for __d_i in __d['D'])}))

            __rd_es += [(__ex + __rd_n, __ey + __rd_n) for __ex, __ey in __td.edges()]

            __rs.append(__td.graph['s'])
            __rc.append(__td.graph['c'])

        __rd.add_nodes_from(__rd_ns)
        __rd.add_edges_from(__rd_es)
        __rd.graph = {'dt': 'MC', 'c': f"({'|'.join(sorted(__rc))})",
                      's': tuple(sum([0 if __x == None else __x for __x in __rx]) for __rx in zip_longest(*__rs)),
                      'w': sum([__td.graph['w'] for __td in _dg]),
                      'm': sum([__td.graph['m'] for __td in _dg]),
                      'sl': min([__td.graph['sl'] for __td in _dg]),
                      'id': max([__td.graph['id'] for __td in _dg]),
                      'od': max([__td.graph['od'] for __td in _dg]),
                      'jld': max([__td.graph['jld'] for __td in _dg]),
                      'jlh': max([__td.graph['jlh'] for __td in _dg])}
        # DAG_Det(__rd)
        return __rd

    def FC(self, _pd, _sd):
        __rd, __rd_equ = nx.DiGraph(), list()

        __pd_n, __pd_l = _pd.number_of_nodes(), len(_pd.graph['s'])
        __sd_n, __sd_l = _sd.number_of_nodes(), len(_sd.graph['s'])

        __pd_sink = set(__i for __i, __d in _pd.nodes(data=True) if len(__d['S']) == 0)
        __sd_sour = set(__i + __pd_n for __i, __d in _sd.nodes(data=True) if len(__d['P']) == 0)

        __rd_ns, __rd_es = list(), [(__ex, __ey) for __ex, __ey in _pd.edges()] + [(__ex + __pd_n, __ey + __pd_n) for
                                                                                   __ex, __ey in _sd.edges()] + list(
            product(__pd_sink, __sd_sour))

        for __i, __d in _sd.nodes(data=True):
            __nel = f"s_{__d['Equ']}"
            __nel not in __rd_equ and __rd_equ.append(__nel)

            __t_n = (__i + __pd_n,
                     {'d': __d['d'] + __pd_l, 'h': __d['h'], 'Equ': __rd_equ.index(__nel),
                      'P': set(__p_i + __pd_n for __p_i in __d['P']),
                      'S': set(__s_i + __pd_n for __s_i in __d['S']),
                      'A': set(__a_i + __pd_n for __a_i in __d['A']) | {__pd_i for __pd_i in _pd.nodes()},
                      'D': set(__d_i + __pd_n for __d_i in __d['D'])})

            if len(__d['P']) == 0:
                __t_n[1]['P'] |= __pd_sink

            __rd_ns.append(__t_n)

        for __i, __d in _pd.nodes(data=True):
            __nel = f"p_{__d['Equ']}"
            __nel not in __rd_equ and __rd_equ.append(__nel)

            __t_n = (__i,
                     {'d': __d['d'], 'h': __d['h'] + __sd_l, 'Equ': __rd_equ.index(__nel),
                      'P': set() | __d['P'],
                      'S': set() | __d['S'],
                      'A': set() | __d['A'],
                      'D': set() | __d['D'] | {__sd_i + __pd_n for __sd_i in _sd.nodes()}})

            if len(__d['S']) == 0:
                __t_n[1]['S'] |= __sd_sour

            __rd_ns.append(__t_n)

        __rd.add_nodes_from(__rd_ns)
        __rd.add_edges_from(__rd_es)

        __rd.graph |= {'dt': 'MS',
                       's': _pd.graph['s'] + _sd.graph['s'],
                       'c': f"{_pd.graph['c']}-{_sd.graph['c']}",
                       'm': _pd.graph['m'] + _sd.graph['m'] + len(__pd_sink) * len(__sd_sour),
                       'w': max(_pd.graph['w'], _sd.graph['w']),
                       'sl': _pd.graph['sl'] + _sd.graph['sl'],
                       'id': max(_pd.graph['id'], _sd.graph['id'], len(__pd_sink)),
                       'od': max(_pd.graph['od'], _sd.graph['od'], len(__sd_sour)),
                       'jld': max(_pd.graph['jld'], _sd.graph['jld'],
                                  __pd_l + 1 - min(__rd.nodes[__pi]['d'] for __pi in __pd_sink)),
                       'jlh': max(_pd.graph['jlh'], _sd.graph['jlh'],
                                  __sd_l + 1 - min(__rd.nodes[__si]['h'] for __si in __sd_sour))}
        # DAG_Det(__rd)
        return __rd

    def ChainDAG(self, _n: int):
        assert _n > 1
        __rd = nx.DiGraph()
        __rd.add_nodes_from([(_i, {'d': 1 + _i, 'h': _n - _i, 'Equ': _i,
                                   'P': set([_i - 1, ]) if 0 < _i else set(),
                                   'A': set(range(_i)) if 0 < _i else set(),
                                   'S': set([_i + 1, ]) if _i < _n - 1 else set(),
                                   'D': set(range(_i + 1, _n)) if _i < _n - 1 else set()}) for _i in range(_n)])
        __rd.add_edges_from([(_i - 1, _i) for _i in range(1, _n)])
        __rd.graph = {'s': tuple([1 for _ in range(_n)]), 'c': f"({'-'.join('t' * _n)})", 'dt': 'MS',
                      'id': 1, 'od': 1, 'jld': 1, 'jlh': 1, 'sl': _n, 'm': _n - 1, 'w': 1}
        return __rd

    def NullDAG(self, _n: int):
        # assert _n > 1
        __rd = nx.DiGraph()
        __rd.add_nodes_from([(_i, {'d': 1, 'h': 1, 'Equ': 0,
                                   'P': set(),
                                   'S': set(),
                                   'A': set(),
                                   'D': set()}) for _i in range(_n)])
        __rd.graph = {'s': (_n,), 'c': f"({'|'.join('t' * _n)})", 'dt': 'MC',
                      'id': 0, 'od': 0, 'jld': 0, 'jlh': 0, 'sl': 1, 'm': 0, 'w': _n}
        return __rd

    def certi_wadjm(self, _wadjm: np.ndarray):
        _wdiag = np.diagonal(_wadjm).T
        __cal = np.argsort(ig.Graph.Adjacency(_wadjm, mode=ig.ADJ_DIRECTED).canonical_permutation(color=_wdiag))
        return _wadjm[__cal][:, __cal]


if __name__ == "__main__":
    # LocalPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), '')
    LocalPath = "C:\\Users\James\PycharmProjects\FT-DAG-main\FT-DAG-main\data"
    # print(''.join(['# ' for _ in range(50)]))
    # print(f"#\tCurrent time:\t{datetime.now()}\n#\tCPU_NUM:\t{cpu_count()}\n#\tRoot path:\t{LocalPath}")
    # print(''.join(['# ' for _ in range(50)]))

    # parser = ArgumentParser('Run FT-DAG')
    # parser.add_argument('-db', '--db', type=str, required=True, help="input the name of DB", default='DATA.db')
    # args = parser.parse_args()
    # print(args.src)

    FilePath = os.path.join(LocalPath, 'DATA.db')
    print(123)
    print(FilePath)
    print(123)
    if os.path.exists(FilePath):
        os.remove(FilePath)
        # TCONN = sqlite3.connect(':memory:')
    TCONN = sqlite3.connect(FilePath)
    TCONN.row_factory = sqlite3.Row
    TCURS = TCONN.cursor()
    TCURS.execute(f"DROP TABLE IF EXISTS GTable;")
    TCURS.execute(
        f"CREATE TABLE IF NOT EXISTS GTable ({', '.join([f'{_i} {_t[0]} {_t[1]}' for _i, _t in GTableHead.items()])});")

    # print('E3.4')
    for _n in Interval(3, 8).iter():
        _dn = 0
        _rt = 0

        # Default *
        _lr = Interval(1, _n)

        # E1.1
        # _mr = Interval(0, _n)
        # E1.2
        # _mr = Interval(0, int(_n / 3))
        # E1.3
        # _mr = Interval(math.ceil(_n / 3), int(2 * _n / 3))
        # E1.4
        # _mr = Interval(math.ceil(2 * _n / 3), _n)

        # E2.1
        # _lr, _sr = Interval(1, 4), Interval(1, math.ceil(_n / 2))
        # E2.2
        # _lr, _sr = Interval(_n - 3, _n), Interval(1, math.ceil(_n / 2))

        # E3.1  I1O2_L N-4
        # _id, _od, _lr = 1, 2, Interval(_n - 4, _n)
        # E3.2  I1O2_L 5
        # _id, _od, _lr = 1, 2, Interval(1, 5)
        # E3.3  I2O1_L N-4
        # _id, _od, _lr = 2, 1, Interval(_n - 4, _n)
        # E3.4  I2O1_L 5
        # _id, _od, _lr = 2, 1, Interval(1, 5)

        # E4.1
        # _lr, _wr = Interval(1, 4), Interval(_n - 3, _n)
        # E4.2
        # _lr, _wr = Interval(_n - 3, _n), Interval(1, 4)

        for _l in _lr.iter():

            _sr = Interval(1, _n - _l + 1)

            for _s in SG(_n, _l, _sr):
                _g = FT_Dag_Generator(_s)
                # _g = FT_Dag_Generator(_s, _mr=_mr)
                # _g = FT_Dag_Generator(_s, _wr=_wr)
                # _g = FT_Dag_Generator(_s, _id=_id, _od=_od)
                _st = time.time()
                for _d in _g.gen():
                    _dn += 1
                    # print(_d)
                _rt += time.time() - _st
                del _g
        # print(f"n:{n}\n\tdn:{dn}\n\time:{_rt:.6f}")
        print(f"{_n},\t{_dn},\t{_rt :.6f}")

    TCURS.close()
    TCONN.commit()
    TCONN.close()
