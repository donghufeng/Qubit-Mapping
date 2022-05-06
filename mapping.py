from functools import partial
from mindquantum import *
import numpy as np
from typing import List
import subprocess
from mindquantum.utils.type_value_check import _check_int_type, _check_value_should_not_less
from mindquantum.utils.type_value_check import _check_value_should_between_close_set, _check_input_type
from mindquantum.algorithm.nisq.qaoa.max_cut_ansatz import _check_graph, _get_graph_act_qubits_num
from mindquantum.algorithm.nisq.qaoa.max_cut_ansatz import _get_graph_act_qubits


class Result:
    pass


def a_k(init_a, k, nstep):
    alpha = 0.602
    return init_a / (k + 1 + 0.01 * nstep)**alpha


def c_k(init_c, k, nstep):
    gamma = 0.101
    return init_c / (k + 1)**gamma


def d_k(**kwargs):
    return np.random.choice([-1, 1], kwargs['size'])


def spsa(fun, x0, nstep, c, a, eps=1e-5, callback=None):
    res = Result()
    res.x = np.array(x0)
    i = 0
    prev = None
    this = None
    while True:
        ak = a_k(a, i, nstep)
        ck = c_k(c, i, nstep)
        delta_k = d_k(size=len(x0))
        xp = res.x + ck * delta_k
        xm = res.x - ck * delta_k
        this = (fun(xp) - fun(xm)) / 2
        gk = this / ck / delta_k
        res.x = res.x - ak * gk
        if callback is not None:
            callback(i, res.x, this, gk)
        i += 1
        if i > nstep:
            break
        if prev is not None:
            if np.abs((this - prev)) < eps:
                break
        prev = this
    res.f = fun(res.x)
    res.step = i
    return res


class MaxCutAnsatz(Ansatz):
    def __init__(self, graph, depth=1):
        _check_int_type('depth', depth)
        _check_value_should_not_less('depth', 1, depth)
        _check_graph(graph)
        super(MaxCutAnsatz, self).__init__('MaxCut',
                                           _get_graph_act_qubits_num(graph),
                                           graph, depth)
        self.graph = graph
        self.all_node = set()
        for edge in self.graph:
            for node in edge:
                self.all_node.add(node)
        self.depth = depth

    def _build_hc(self, graph):
        circ = Circuit()
        for node in graph:
            circ += XX('beta').on(node)
        return circ

    def _build_hb(self, graph):
        circ = Circuit(
            [RY('alpha').on(i) for i in _get_graph_act_qubits(graph)])
        return circ

    @property
    def hamiltonian(self):
        qo = QubitOperator('', 0)
        for node in self.graph:
            qo += (QubitOperator('') -
                   QubitOperator(f'X{node[0]} X{node[1]}')) / 2
        return qo

    def get_partition(self, max_n, weight):
        _check_int_type('max_n', max_n)
        _check_value_should_between_close_set('max_n', 1,
                                              1 << self._circuit.n_qubits,
                                              max_n)
        sim = Simulator('projectq', self._circuit.n_qubits)
        sim.apply_circuit(self._circuit, weight)
        qs = sim.get_qs()
        idxs = np.argpartition(np.abs(qs), -max_n)[-max_n:]
        partitions = [
            bin(i)[2:].zfill(self._circuit.n_qubits)[::-1] for i in idxs
        ]
        res = []
        for partition in partitions:
            left = []
            right = []
            for i, j in enumerate(partition):
                if j == '0':
                    left.append(i)
                else:
                    right.append(i)
            res.append([left, right])
        return res

    def get_cut_value(self, partition):
        _check_input_type('partition', list, partition)
        if len(partition) != 2:
            raise ValueError(
                f"Partition of max-cut problem only need two parts, but get {len(partition)} parts"
            )
        for part in partition:
            _check_input_type('each part of partition', list, part)
        all_node = set()
        for part in partition:
            for node in part:
                all_node.add(node)
        if all_node != self.all_node:
            raise ValueError(
                f"Invalid partition, partition nodes are different with given graph."
            )
        cut_value = 0
        for edge in self.graph:
            node_left, node_right = edge
            for n in edge:
                if n not in partition[0] and n not in partition[1]:
                    raise ValueError(
                        f'Invalid partition, node {n} not in partition.')
                if n in partition[0] and n in partition[1]:
                    raise ValueError(
                        f'Invalid partition, node {n} in both side of cut.')
            if (node_left in partition[0] and node_right in partition[1] or
                    node_left in partition[1] and node_right in partition[0]):
                cut_value += 1
        return cut_value

    def _implement(self, graph, depth):
        """Implement of max cut ansatz."""
        self._circuit = UN(H, _get_graph_act_qubits(graph))
        for d in range(depth):
            self._circuit += CPN(self._build_hc(graph), {'beta': f'beta_{d}'})
            self._circuit += CPN(self._build_hb(graph),
                                 {'alpha': f'alpha_{d}'})


def prepare_circ_for_mapping(circ: List[BasicGate], file_name: str):
    n_qubits = circ.n_qubits
    size = len(circ)
    out = [f"{n_qubits} {size}"]
    for idx, g in enumerate(circ):
        local_qubits = [i for i in g.obj_qubits]
        local_qubits.extend(g.ctrl_qubits)
        g_qubits = len(local_qubits)

        if g_qubits == 1:
            out.append(f"X {local_qubits[0]} {idx}")
        elif g_qubits == 2:
            out.append(f"CNOT {local_qubits[0]} {local_qubits[1]} {idx}")
        else:
            raise RuntimeError(
                f"Mapping only work for gate less than two qubits, but get gate {g}"
            )
    with open(file_name, 'w') as f:
        f.write('\n'.join(out))


def prepare_coupling_graph(edges: List[List[int]], file_name: str):
    out = [""]
    n_qubits = 0
    n_edges = len(edges)
    for edge in edges:
        out.append(f"{edge[0]} {edge[1]}")
        n_qubits = max(n_qubits, max(edge))
    n_qubits += 1
    out[0] = f"{n_qubits} {n_edges}"
    with open(file_name, 'w') as f:
        f.write('\n'.join(out))


def circ_mapping(circ, coupling_graph):
    prepare_circ_for_mapping(circ, "logical_circuit.txt")
    prepare_coupling_graph(coupling_graph, "coupling_graph.txt")

    status, ret = subprocess.getstatusoutput("./a.out")
    if status:
        raise RuntimeError("Mapping failed.")
    ret = ret.split('\n')
    final_mapping = []
    t = 0
    for i in ret[-1].split(' ')[4:-1]:
        if t % 2 == 0:
            final_mapping.append([])
        t += 1
        if i.startswith('('):
            final_mapping[-1].append(int(i[1:-1]))
        else:
            final_mapping[-1].append(int(i[:-1]))
    new_circ = Circuit()
    for g in ret[:-1]:
        g = g.split(' ')
        tag = g[-1]
        q0 = int(g[1])
        q1 = int(g[2])
        if tag.startswith("SWAP"):
            new_circ += SWAP.on([q0, q1])
        else:
            ori = circ[int(tag)]
            if q0 == q1:
                new_circ += ori.on(q0)
            else:
                if ori.n_qubits == 1:
                    new_circ += ori.on(q0, q1)
                else:
                    new_circ += ori.on([q0, q1])
    return new_circ, {i: j for i, j in final_mapping}


def operator_mapping(operator: QubitOperator, final_mapping):
    out = operator.__class__()
    terms = {}
    for term, coeff in operator.terms.items():
        new_o = []
        for idx, o in term:
            new_o.append((final_mapping[idx], o))
        terms[tuple(new_o)] = coeff + 0
    out.terms = terms
    return out


def maxcut_fun(x, sim, circ, ham):
    sim.reset()
    sim.apply_circuit(circ, x)
    f = np.real(sim.get_expectation(ham))
    return f


edges = [(0, 1), (1, 2), (2, 3), (0, 3), (1, 4), (2, 4)]
maxcut = MaxCutAnsatz(edges, 2)
coupling_graph = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
]

new_circ, final_mapping = circ_mapping(maxcut.circuit, coupling_graph)
new_h = operator_mapping(maxcut.hamiltonian, final_mapping)
h = Hamiltonian(-new_h)
init = np.random.random(len(new_circ.params_name))
sim = Simulator('projectq', new_circ.n_qubits)
res = spsa(partial(maxcut_fun, sim=sim, circ=new_circ, ham=h), init, 100, 0.2,
           0.2)