import os
import sqlite3
import networkx as nx
import matplotlib.pyplot as plt
import csv

# Ensure these modules are in the same directory
from _Shape import SG
from _Interval import Interval
import FT_DAG
from FT_DAG import FT_Dag_Generator


# =========================
# DAG Drawing Function
# =========================
def draw_dag(dag, count):
    plt.figure(figsize=(10, 6))

    topo_nodes = list(nx.topological_sort(dag))
    mapping = {old: i for i, old in enumerate(topo_nodes)}
    relabel_dag = nx.relabel_nodes(dag, mapping)

    # Layer comes from 'd'
    for n in relabel_dag.nodes():
        relabel_dag.nodes[n]['layer'] = relabel_dag.nodes[n]['d']

    pos = nx.multipartite_layout(relabel_dag, subset_key="layer")

    node_colors = [
        '#ff9999' if relabel_dag.nodes[n]['d'] == 1 else '#66b3ff'
        for n in relabel_dag.nodes()
    ]

    labels = {n: n + 1 for n in relabel_dag.nodes()}

    nx.draw(
        relabel_dag,
        pos,
        with_labels=True,
        labels=labels,
        node_color=node_colors,
        node_size=800,
        font_weight='bold',
        arrows=True,
        arrowsize=15
    )

    sinks = [n + 1 for n in relabel_dag.nodes() if relabel_dag.out_degree(n) == 0]
    plt.title(f"DAG #{count} | Sinks: {sinks} | Shape: {dag.graph['s']}")
    plt.show()


def get_next_dag_id(csv_filename):
    if not os.path.exists(csv_filename):
        return 1

    max_id = 0
    with open(csv_filename, mode='r', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            try:
                max_id = max(max_id, int(row[0]))
            except ValueError:
                continue

    return max_id + 1

# =========================
# Main Task
# =========================
def run_constrained_task():
    TCONN = sqlite3.connect(':memory:')
    TCONN.row_factory = sqlite3.Row

    FT_DAG.TCONN = TCONN
    FT_DAG.TCURS = TCONN.cursor()

    cols = [
        f"{col_name} {col_type[0]} {col_type[1]}"
        for col_name, col_type in FT_DAG.GTableHead.items()
    ]

    FT_DAG.TCURS.execute(
        f"CREATE TABLE IF NOT EXISTS GTable ({', '.join(cols)});"
    )

    # Parameters
    N_TOTAL = 12
    INPUT_SIZE = 6
    IN_DEGREE = 2
    REQUIRED_SINKS = 2

    csv_filename = "dag_database_wide_62.csv"
    next_dag_id = get_next_dag_id(csv_filename)

    print(f"Appending DAGs starting from ID {next_dag_id}")

    with open(csv_filename, mode='a', newline='') as f:
        writer = csv.writer(f)

        for l in range(2, N_TOTAL - INPUT_SIZE + 2):
            for s in SG(N_TOTAL, l, Interval(1, N_TOTAL)):

                if s[0] != INPUT_SIZE:
                    continue

                generator = FT_Dag_Generator(_s=s, _id=IN_DEGREE)

                for dag in generator.gen():

                    if any(dag.degree(n) == 0 for n in dag.nodes()):
                        continue

                    internal_nodes = [
                        n for n, d in dag.nodes(data=True) if d['d'] > 1
                    ]
                    if not internal_nodes:
                        continue

                    if not all(len(dag.pred[n]) == IN_DEGREE for n in internal_nodes):
                        continue

                    sinks = [n for n in dag.nodes() if dag.out_degree(n) == 0]
                    if len(sinks) != REQUIRED_SINKS:
                        continue

                    # ===== VALID DAG =====
                    dag_id = next_dag_id
                    next_dag_id += 1

                    print(f"✔ Appending DAG #{dag_id}")

                    # draw_dag(dag, dag_id)

                    topo_nodes = list(nx.topological_sort(dag))
                    mapping = {old: i + 1 for i, old in enumerate(topo_nodes)}

                    row_data = [dag_id]
                    for node in topo_nodes:
                        preds = list(dag.pred[node])
                        if len(preds) == 2:
                            row_data.extend([
                                mapping[preds[0]],
                                mapping[preds[1]],
                                mapping[node],
                            ])

                    writer.writerow(row_data)

    print(f"\nFinished. Last DAG ID written: {next_dag_id - 1}")


# =========================
# Entry Point
# =========================
if __name__ == "__main__":
    run_constrained_task()
