#include<iostream>
#include<fstream>
#include<random>
#include<ctime>
#include<cstdlib>
#include<string>

using namespace std;

void generate_logical_circuit_random(int n, int g, string filename) {
    ofstream fout;
    fout.open(filename);

    fout << n << " " << g << endl;

    for (int id = 0; id < g; ++id) {
        int type = rand() % 2 == 0 ? 1 : 2;
        if (n == 1) type = 1;
        if (type == 1) {
            int p = rand() % n;

            fout << "X" << " ";
            fout << p << " ";
            fout << "g" + to_string(id) << endl;            
        }
        else {
            int p = rand() % n;
            int q = rand() % n;
            while (p == q) {
                q = rand() % n;
            }

            fout << "CNOT" << " ";
            fout << p << " " << q << " ";
            fout << "g" + to_string(id) << endl;
        }
    }
    fout.close();
}

void generate_coupling_graph_linear(int n, string filename) {
    ofstream fout;
    fout.open(filename);
    
    fout << n << " " << n-1 << endl;
    for (int i = 0; i < n-1; ++i)
        fout << i << " " << i+1 << endl;
    
    fout.close();
}

void generate_coupling_graph_grid(int m, int n, string filename) {
    ofstream fout;
    fout.open(filename);

    auto p = [=](int x, int y) {
        return x * n + y;
    };

    fout << m * n << " " << (n-1) * m + (m-1) * n << endl;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n-1; ++j)
            fout << p(i, j) << " " << p(i, j+1) << endl;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < m-1; ++i)
            fout << p(i, j) << " " << p(i+1, j) << endl;    

    fout.close();
}

int main() {
    int num_logical = 20;
    int num_gates = 100;

    int num_physical = 20;
    int m = 4;
    int n = 6;

    srand(time(0));

    generate_logical_circuit_random(num_logical, num_gates, "logical_circuit.txt");

    // generate_coupling_graph_linear(num_physical, "coupling_graph.txt");

    generate_coupling_graph_grid(m, n, "coupling_graph.txt");

    return 0;
}