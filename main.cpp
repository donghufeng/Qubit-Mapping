#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cstdio>
#include<list>
#include<set>

using namespace std;

struct Gate {
    string type;    // gate's type, such as CNOT or X
    int q1, q2;     // logical qubit number, if 
                    // gate is single gate, q2=-1
    string tag;     // user defined tag
};

vector<vector<int>> get_circuit_DAG(int n, const vector<Gate> & gates) {
    int m = gates.size();
    vector<int> last(n, -1);
    vector<vector<int>> DAG(m);

    for (int i = 0; i < gates.size(); ++i) {
        int q1 = gates[i].q1;
        int q2 = gates[i].q2;

        if (last[q1] != -1) 
            DAG[last[q1]].push_back(i);
        last[q1] = i;

        if (q2 != -1) {
            if (last[q2] != -1) 
                DAG[last[q2]].push_back(i);
            last[q2] = i;
        }
    }

    return DAG;
}

pair<int, vector<Gate>> get_logical_circuit_from_file(string filename="logical_circuit.txt") {
    ifstream fin;
    fin.open(filename);

    int n, m;
    fin >> n >> m;
    
    vector<Gate> gates(m);
    for (int i = 0; i < m; ++i) {
        fin >> gates[i].type;
        if (gates[i].type == "CNOT")
            fin >> gates[i].q1 >> gates[i].q2;
        else {
            fin >> gates[i].q1;
            gates[i].q2 = -1;
        }
        fin >> gates[i].tag;
    }

    fin.close();

    return {n, gates};
}

void show_DAG(const vector<vector<int>> & DAG, const vector<Gate> & gates) {
    printf("show DAG:\n");
    for (int i = 0; i < DAG.size(); ++i) {
        printf("%d ", i);
        cout << gates[i].tag;
        printf(" : ");
        for (int y : DAG[i]) {
            printf("(%d, ", y);
            cout << gates[y].tag << "), ";
        }
        printf("\n");
    }
}


class SABRE {
private:
    int num_logical;
    int num_physical;
    vector<Gate> gates;         // logical circuit
    vector<vector<int>> G;      // physical coupling graph

    vector<vector<int>> DAG;    // DAG of logical circuit
    vector<vector<int>> RDAG;   // reverse graph of DAG

    vector<vector<int>> D;      // nearest neighbor distance

    inline bool is_executable(const vector<int> & pi, int g) const {
        if (gates[g].type == "CNOT") {
            int p = pi[gates[g].q1], q = pi[gates[g].q2];
            for (int x : G[p])
                if (q == x)
                    return true;
            return false;
        }
        else 
            return true;
    }

    inline vector<int> get_reverse_pi(const vector<int>& pi) const {
        vector<int> rpi(pi.size());
        for (int i = 0; i < pi.size(); ++i)
            rpi[pi[i]] = i;
        return rpi;
    }

    set<pair<int, int>> obtain_SWAPs(const list<int>& F, const vector<int>& pi) const {
        set<pair<int, int>> ret;
        // vector<int> rpi = get_reverse_pi(pi);

        for (int g : F) {
            int x = pi[gates[g].q1];
            int y = pi[gates[g].q2];
            for (int z : G[x])
                ret.insert({min(x, z), max(x, z)});
            for (int z : G[y])
                ret.insert({min(y, z), max(y, z)});
        }
        return ret;
    }

public:
    SABRE(int num_logical, vector<Gate>& gates, int num_physical, vector<vector<int>>& G)
        : num_logical(num_logical), num_physical(num_physical),
        gates(gates), G(G) {
        this->DAG = get_circuit_DAG(num_logical, gates);

        this->RDAG = vector<vector<int>> (this->DAG.size());
        for (int x; x < DAG.size(); ++x) {
            for (int y : DAG[x])
                RDAG[y].push_back(x);
        }

        // get D using Floyd algorithm
        {
            int n = num_physical;
            this->D = vector<vector<int>> (n, vector<int>(n));
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    D[i][j] = i == j ? 0 : 1e9;
            for (int i = 0; i < n; ++i)
                for (int j : G[i])
                    D[i][j] = D[j][i] = 1;
            for (int k = 0; k < n; ++n)
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        D[i][j] = min(D[i][j], D[i][k] + D[k][j]);    
        }
    }

    vector<Gate> heuristic_search(vector<int> & pi, const vector<vector<int>>& DAG) {
        vector<Gate> ans;
        vector<int> indeg(DAG.size(), 0);
        for (int i = 0; i < DAG.size(); ++i)
            for (int j : DAG[i])
                indeg[j]++;
        
        list<int> F;
        for (int i = 0; i < DAG.size(); ++i)
            if (indeg[i] == 0)
                F.push_back(i);
        
        while (!F.empty()) {
            vector<int> executable_list;
            for (auto it = F.begin(); it != F.end(); ++it) {
                if (is_executable(pi, *it)) {
                    executable_list.push_back(*it);
                }
            }

            if (!executable_list.empty()) {
                for (auto it = F.begin(); it != F.end(); ) {
                    if (is_executable(pi, *it)) {
                        int x = *it;
                        it = F.erase(it);
                        ans.push_back(gates[x]);
                        for (int y : DAG[x]) {
                            --indeg[y];
                            if (indeg[y] == 0)
                                F.push_back(y);
                        }
                    }
                    else 
                        ++it;
                }
            }
            else {
                auto candidate_SWAPs = obtain_SWAPs(F, pi);
                auto rpi = get_reverse_pi(pi);

                int min_score = 0x7fffffff;
                pair<int, int> min_SWAP;
                for (auto SWAP : candidate_SWAPs) {
                    int x = SWAP.first, y = SWAP.second;
                    int p = rpi[x], q = rpi[y];
                    
                    auto tmp = pi;
                    swap(tmp[p], tmp[q]);
                    int score = H(F, DAG, tmp, D, SWAP);
                    if (score < min_score) {
                        min_score = score;
                        min_SWAP = SWAP;
                    }
                }
                int p = rpi[min_SWAP.first], q = rpi[min_SWAP.second];
                swap(pi[p], pi[q]);
            }
        }

    }

};


int main() {
    auto par = get_logical_circuit_from_file("logical_circuit_2.txt");
    auto n = par.first;
    auto gates = par.second;

    auto DAG = get_circuit_DAG(n, gates);
    show_DAG(DAG, gates);

    return 0;
}