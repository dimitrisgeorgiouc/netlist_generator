#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
RTSS 2023 DAG Generator
"""

import os
import sys
import copy
import time
import networkx as nx

from itertools import combinations, combinations_with_replacement

sys.path.append(os.path.dirname(__file__))

from _Interval import *

""" (1) combination A004250 """
def combination(_n:int, 
                _l:int,
                _sr:tuple=(0, float('Inf'))):

    if _l * _sr[0] <= _n <= _l * _sr[1] and 1 <= _l <= _n:

        # case 1 null dag
        if _l == 1:
            yield (_n, )

        # case 2 chain dag
        elif _l == _n:
            yield tuple(1 for _ in range(_l))

        # case 3 shape deter+
        # elif _sr[0] == _sr[1]:
        #     yield tuple(int(_n / _l) for _ in range(_l))
        
        # case 4 other：
        else:
            for __Xl in range(_sr[0], _sr[1] + 1):
                __Xs = _n - __Xl
                __ls = _l - 1
                if __ls * _sr[0] <= __Xs <= __ls * _sr[1]:
                    for __sub_s in combination(__Xs, __ls, (_sr[0], __Xl)):
                        yield (__Xl,) + __sub_s

""" (2) permutation xx """
def permutation(_c:tuple):
    __cset = set(_c)
    if len(__cset) == 1:
        yield _c
    else:
        for __Xl in __cset:
            __index = _c.index(__Xl)
            for __sub_c in permutation(_c[:__index] + _c[__index + 1:]):
                yield (__Xl,) + __sub_c

""" (3) connection xx """
def shape_list_trance(_s:tuple):
    __n = sum(_s)
    __ni_list = list(range(__n))
    __rs_list = []
    for __sx in _s:
        __rs_list.append(__ni_list[:__sx])
        del __ni_list[:__sx]
    return __rs_list

def gen_mine_new(_s:tuple):
    shape_list = shape_list_trance(_s)
    dag_list = [nx.DiGraph()]
    for level_id, self_node_list in enumerate(shape_list):
        if (level_id + 1) == len(shape_list):
            for dag_x in dag_list:
                p_nodes = [nodex for nodex in dag_x.nodes() if len(list(dag_x.successors(nodex))) == 0]
                dag_x.add_nodes_from([(self_node_list[0], {'level_num': level_id})])
                for p_node_x in p_nodes:
                    dag_x.add_edge(p_node_x, self_node_list[0])
        else:
            temp_dag_list = []
            for dag_x in dag_list:
                for __rdx in shape_dag_generator(dag_x, self_node_list, level_id):
                    temp_dag_list.append(__rdx)
            dag_list = temp_dag_list
    for ret_dagx in dag_list:
        yield ret_dagx

def shape_dag_generator(dag_x, self_node_list, level_num):
    for level_id in range(level_num):
        if level_id == level_num - 1:
            up_same_level_node_iso_label_comput(dag_x, level_id)
        down_same_level_node_iso_label_comput(dag_x, level_id)
    total_label_dict = {}
    for node_x in dag_x.nodes(data=True):
        node_x_label = (node_x[1]['level_num'], node_x[1]['up_iso_label'],node_x[1]['down_iso_label'])
        if node_x_label in total_label_dict:
            total_label_dict[node_x_label].append(node_x[0])
        else:
            total_label_dict[node_x_label] = [node_x[0]]

    last_level_node_list = [node_x[0] for node_x in dag_x.nodes(data=True) if node_x[1]['level_num'] == level_num - 1]
    last_level_node_id_enumerate_list = []
    for sn_num in range(len(last_level_node_list)):
        temp_id_enumerate_list_1 = list(combinations(last_level_node_list, sn_num + 1))
        temp_label_enumerate_list = list(set([tuple([(dag_x.nodes[temp_id_x]['level_num'], dag_x.nodes[temp_id_x]['up_iso_label'],dag_x.nodes[temp_id_x]['down_iso_label']) for temp_id_x in temp_id_list])
                                            for temp_id_list in temp_id_enumerate_list_1]))
        temp_id_enumerate_list_2 = []
        for temp_label_list in temp_label_enumerate_list:
            temp_total_label_dict = copy.deepcopy(total_label_dict)
            temp_id_enumerate_list_2.append([temp_total_label_dict[temp_label_x].pop(0) for temp_label_x in temp_label_list])
        last_level_node_id_enumerate_list += temp_id_enumerate_list_2
    pnode_list_enumerate = []
    for last_level_node_enumerate_x in last_level_node_id_enumerate_list:
        sample_dag = copy.deepcopy(dag_x)
        rem_set = set(last_level_node_list)
        for last_level_node_x in last_level_node_enumerate_x:
            rem_set.update(nx.ancestors(sample_dag, last_level_node_x))
        sample_dag.remove_nodes_from(rem_set)
        pred_node_opt_list = list(nx.antichains(sample_dag, topo_order=None))
        for pred_node_opt_x in pred_node_opt_list:
            pred_node_opt_x += last_level_node_enumerate_x
        pnode_list_enumerate += pred_node_opt_list
    temp_dag_x = copy.deepcopy(dag_x)
    temp_dag_x.add_nodes_from([(self_node_x, {'level_num':level_num}) for self_node_x in self_node_list])
    if len(pnode_list_enumerate) == 0:
        yield temp_dag_x
    else:
        edge_p_list = list(combinations_with_replacement(pnode_list_enumerate, len(self_node_list)))
        for edge_p_list_x in edge_p_list:
            temp_dag_list_x = copy.deepcopy(temp_dag_x)
            for self_node_id, edges_p_x in enumerate(edge_p_list_x):
                for edge_p_x in edges_p_x:
                    temp_dag_list_x.add_edge(edge_p_x, self_node_list[self_node_id])
            yield temp_dag_list_x

def up_same_level_node_iso_label_comput(dag_x, level_id):
    self_level_node_list = [node_x[0] for node_x in dag_x.nodes(data=True) if node_x[1]['level_num'] == level_id]
    up_iso_node_list = [[self_level_node_list.pop()]]
    for node_x in self_level_node_list:
        t_step = True
        sn_subg = dag_x.subgraph(list(nx.ancestors(dag_x, node_x)) + [node_x])
        for node_id_list in up_iso_node_list:
            tsn_subg = dag_x.subgraph(list(nx.ancestors(dag_x, node_id_list[0])) + [node_id_list[0]] )
            # if nx.isomorphism.GraphMatcher(sn_subg, tsn_subg).is_isomorphic():
            if nx.isomorphism.DiGraphMatcher(sn_subg, tsn_subg).is_isomorphic():
                node_id_list.append(node_x)
                t_step = False
                break
        if t_step:
            up_iso_node_list.append([node_x])
    for up_iso_label, node_list in enumerate(up_iso_node_list):
        for node_x in node_list:
            dag_x.nodes[node_x]['up_iso_label'] = up_iso_label

def down_same_level_node_iso_label_comput(dag_x, level_id):
    self_level_node_list = [node_x[0] for node_x in dag_x.nodes(data=True) if node_x[1]['level_num'] == level_id]
    down_iso_node_list = [[self_level_node_list.pop(0)]]
    for node_x in self_level_node_list:
        t_step = True
        sn_subg = dag_x.subgraph(list(dag_x.successors(node_x)) + [node_x])

        for node_id_list in down_iso_node_list:
            tsn_subg = dag_x.subgraph(list(dag_x.successors(node_id_list[0])) + [node_id_list[0]] )
            # if nx.isomorphism.GraphMatcher(sn_subg, tsn_subg).is_isomorphic():
            if nx.isomorphism.DiGraphMatcher(sn_subg, tsn_subg).is_isomorphic():
                node_id_list.append(node_x)
                t_step = False
                break
        if t_step:
            down_iso_node_list.append([node_x])
    for down_iso_label, node_list in enumerate(down_iso_node_list):
        for node_x in node_list:
            dag_x.nodes[node_x]['down_iso_label'] = down_iso_label

def width_of_dag(_cdag):
    temp_G2 = nx.DiGraph()  # 创建一个二分图
    for edge_x in _cdag.edges():
        temp_G2.add_node('p' + str(edge_x[0]), bipartite=0)
        temp_G2.add_node('d' + str(edge_x[1]), bipartite=1)
        temp_G2.add_edge('p' + str(edge_x[0]), 'd' + str(edge_x[1]))
    u = [n for n in temp_G2.nodes if temp_G2.nodes[n]['bipartite'] == 0]
    matching = nx.algorithms.bipartite.matching.hopcroft_karp_matching(temp_G2, top_nodes=u)
    return _cdag.number_of_nodes() - int(len(matching) / 2)


if __name__ == "__main__":
    print('E4.1')
    for _n in range(3, 100):
        # Source
        # _lr, _wr = Interval(_n - 2, _n), Interval(3, 3)
        # Experiment 4.1 
        _lr, _wr = Interval(1, 4), Interval(_n - 3, _n)
        # Experiment 4.2
        # _lr, _wr = Interval(_n - 3, _n), Interval(1, 4)

        _cn, _sn, _dn = 0, 0, 0
        _st = time.time()

        for _l in _lr.iter():
            for _cx in combination(_n=_n, _l=_l, _sr=(1, _n)):                
                _cn += 1
                for _sx in permutation(_cx):
                    _sn += 1
                    for _dx in gen_mine_new(_sx):
                        _dc = nx.transitive_closure_dag(_dx)  
                        _w = width_of_dag(_dc)
                        if _w in _wr:
                            _dn += 1

        _et = time.time()
        # print(f"{_n}_\t{_cn}_\t{_sn}_\t{_dn}_\t{et - st :.8f}")
        print(f"{_n},\t{_dn},\t{_et - _st :.8f}")
