# FT-DAG

FT-DAG is an efficient and formally verified full-topology DAG generator that is able to control all major parameters, including:

- in-degree: _id; 
- out-degree: _od;
- jump level: _jlh;
- jump layer: _jld; 
- width range: _wr;
- longest length: _ll;
- shortest length: _ls;
- shape value range: _sr;
- the number of nodes: _n;
- the number of edges range: _mr;

Experiments show that when the number of nodes is larger than 20, FT-DAG provides at least two orders of magnitude speedup compared to the state of the art and more orders to other generators. FT-DAG scales to 100 nodes in a typical industrial case study within hours.


This work also reproduces some variants of existing DAG generators as a control group, as follows:
- GNM[1];
- LBL(layer-by-layer)[2];
- FIO(fan-in/fan-out)[3];
- SoTA(state of the art)[4]:

 
## References
[1] P. Erdos and A. Renyi. On random graphs I.Publicationes Mathematicae Debrecen, 6:290–297, 1959.

[2] Tobita T, Kasahara H. A standard task graph set for fair evaluation of multiprocessor scheduling algorithms[J]. Journal of Scheduling, 2002, 5(5): 379-394.

[3] R. P. Dick, D. L. Rhodes, and W. Wolf. TGFF: Task Graphs For Free. In Proceedings of the 6th International Workshop on Hardware/Software Codesign, pages 97–101, Washington, DC, USA, Mar. 1998. IEEE Computer Society.

[4] Fang Y, Zhao S, Guo Y, et al. Brief Industry Paper: A DAG Generator with Full Topology Coverage[C]//2023 IEEE Real-Time Systems Symposium (RTSS). IEEE, 2023: 506-511.

[5] McKay B D, Piperno A. Practical graph isomorphism, II[J]. Journal of symbolic computation, 2014, 60: 94-112.
