# README

This program is for qubit mapping problem based on SABRE and HA algorithm.

## Input Format

### Logical Circuit

Logical circuit should contain:

1. number of qubits
2. quantum gates from left to right
3. Gate's type is a `string` such as "CNOT" or "X", every gate also has a `string` tag. 
4. Qubit's id is from `0` to `n-1`.

For example: `logical_circuit.txt`

```
3 4
CNOT 0 1 tag1
X 0 tag2
Y 1 tag3
CNOT 1 2 tag2
```

### Physical Connectivity

Physical qubit coupling graph should contaiin:

1. number of physical qubits
2. coupling graph

`graph.txt` format:

1. First line is `n m`, `n` is qubit number, `m` is edge number in coupling graph.
2. The next `m` lines' format is `u v`, which means an undirected edge between `u` and `v`.
3. Qubit's id is from `0` to `n-1`.

For example:

```
4 4
0 1
1 2
2 3
3 0
```

## Output Format

Output should contain three parts:

1. Initial mapping
2. Generated physical circuit
3. Final mapping