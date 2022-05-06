// #define DEBUG
// #define BENCHMARK

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cstdio>
#include<list>
#include<set>
#include<random>
#include<chrono>
#include<algorithm>

using namespace std;

// using time as random seed
auto seed = std::chrono::system_clock::now().time_since_epoch().count();
// generate a global random engine
auto engine = std::default_random_engine(seed);

struct Gate {
    string type;    // gate's type, such as CNOT or X
    int q1, q2;     // logical qubit number, if 
                    // gate is single gate, q2=q1
    string tag;     // user defined tag
};

/**
 * @brief Get the circuit DAG object
 * 
 * @param n number of logical qubits
 * @param gates logical circuit from left to right
 * @return vector<vector<int>> DAG of logical circuit
 */
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

        if (gates[i].type == "CNOT" || gates[i].type == "SWAP") {
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
        if (gates[i].type == "CNOT" || gates[i].type == "SWAP")
            fin >> gates[i].q1 >> gates[i].q2;
        else {
            fin >> gates[i].q1;
            gates[i].q2 = gates[i].q1; // if g is single gate, g.q2=g.q1
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

    vector<vector<int>> D;      // nearest neighbor cost

    vector<int> decay;          // decay of each logical qubit

    double W;                   // parameter between F and E
    double delta1;              // decay of a single gate
    double delta2;              // decay of a CNOT gate

    /**
     * @brief judge whether gate[g] can be executable under mapping pi
     * 
     * @param pi current mapping from logical qubit to physical qubit
     * @param g id of gate
     * @return true when pi[g.q1] and pi[g.q2] is neighbor in G
     */
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

    /**
     * @brief Get the mapping from physical qubit to logical qubit
     * 
     * @param pi mapping from logical to physical
     * @return vector<int> mapping from physical to logical
     */
    inline vector<int> get_reverse_pi(const vector<int>& pi) const {
        vector<int> rpi(pi.size());
        for (int i = 0; i < pi.size(); ++i)
            rpi[pi[i]] = i;
        return rpi;
    }

    /**
     * @brief get the candidate SWAP list when there is no executable gate.
     *   If edge (x,z) in G and gate (x,y) in F, then SWAP (x,z) is possible.
     * @param F current front layer
     * @param pi current mapping from logical to physical
     * @return set<pair<int, int>> set of candidate SWAP list, containing physical id.
     */
    set<pair<int, int>> obtain_SWAPs(const list<int>& F, const vector<int>& pi) const {
        set<pair<int, int>> ret;
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

    /**
     * @brief Get the next layer of F in DAG, only considering CNOT gates.
     *   Single gates can always be executed, so there is no need to consider.
     * @param F current front layer
     * @param DAG 
     * @param indeg current in-degree of DAG
     * @return list<int> the next layer of F, ignoring single gates.
     */
    list<int> get_next_layer(const list<int>& F, const vector<vector<int>> & DAG, const vector<int>& indeg) const {
        vector<int> tmp_deg = indeg;
        list<int> ret;
        for (int x : F) {
            for (int y : DAG[x]) {
                tmp_deg[y]--;
                if (gates[y].type == "CNOT") {  // y is CNOT gate
                    if (tmp_deg[y] == 0)
                        ret.push_back(y);
                }
                else {                          // y is single gate
                    for (int z : DAG[y]) {      // find following gate
                        tmp_deg[z]--;
                        if (tmp_deg[z] == 0)
                            ret.push_back(z);
                    }
                }
            }
        }
        return ret;
    }

    /**
     * @brief Get the extended set E
     *   There are many ways to generate E. Here we just use the next layer of F.
     * @param F current front layer
     * @param DAG 
     * @param indeg current in-degree of DAG
     * @return list<int> extended set E
     */
    list<int> get_extended_set(const list<int>& F, const vector<vector<int>>& DAG, const vector<int>& indeg) const {
        return get_next_layer(F, DAG, indeg);
    }

    /**
     * @brief basic heuristic function
     *      H = \sum_{g\in F} D[pi[g.q1]][pi[g.q2]]
     * @param F set of gates' id
     * @param pi mapping from logical to physical
     * @return double 
     */
    double H_basic(const list<int>& F, const vector<int>& pi) const {
        int sum = 0;
        for (int g : F) {
            int q1 = gates[g].q1;
            int q2 = gates[g].q2;
            sum += D[pi[q1]][pi[q2]];
        }
        return sum;
    }

    /**
     * @brief heuristic function considering extended set E
     *   H = 1 / |F| * H_basic(F, pi) + W / |E| * H_basic(E, pi)
     * @param F current front layer
     * @param E extended set 
     * @param pi mapping from logical to physical
     * @return double 
     */
    double H_look_ahead(const list<int>& F, const list<int>& E, const vector<int>& pi) const {
        double s1 = H_basic(F, pi) / (double)F.size();
        if (E.size() == 0)
            return s1;
        else {
            double s2 = H_basic(E, pi) / (double)E.size();
            return s1 + W * s2;     // where 0 <= W <= 1 is a parameter
        }
    }

    /**
     * @brief heuristic function considering trade-off between circuit depth and gates number.
     * 
     * @param F current front layer
     * @param E extended set 
     * @param pi mapping from logical to physical
     * @param SWAP physical SWAP, using logical id
     * @param decay decay of logical qubits
     * @return double 
     */
    double H(const list<int> & F, const list<int>& E, const vector<int>& pi, const pair<int, int> & SWAP, const vector<double>& decay) const {
        // return H_basic(F, pi);
        // return H_look_ahead(F, E, pi);
        return max(decay[SWAP.first], decay[SWAP.second]) * H_look_ahead(F, E, pi);
    }

public:
    /**
     * @brief Construct a new SABRE object
     * 
     * @param num_logical number of logical qubits
     * @param gates logical circuit from left to right
     * @param num_physical number of physical qubits
     * @param G physical coupling graph
     */
    SABRE(int num_logical, vector<Gate>& gates, int num_physical, vector<vector<int>>& G)
        : num_logical(num_logical), num_physical(num_physical), gates(gates), G(G) {
        
        // get DAG and RDAG of logical circuit
        this->DAG = get_circuit_DAG(num_logical, gates);
        this->RDAG = vector<vector<int>> (this->DAG.size());
        for (int x; x < DAG.size(); ++x) {
            for (int y : DAG[x])
                RDAG[y].push_back(x);
        }

        #ifdef DEBUG
            show_DAG(DAG, gates);
            printf("Reverse DAG\n");
            show_DAG(RDAG, gates);
        #endif

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
            for (int k = 0; k < n; ++k)
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        D[i][j] = min(D[i][j], D[i][k] + D[k][j]);    
            
            #ifdef DEBUG
                printf("Show D:\n");
                for (int i = 0; i < n; ++i) {
                    printf("%d: ", i);
                    for (int j = 0; j < n; ++j)
                        printf("%d ", D[i][j]);
                    printf("\n");
                }
            #endif
        }
    }

    /**
     * @brief heuristic search algorithm to generate physical circuit
     * 
     * @param pi Initial mapping from logical to physical
     * @param DAG 
     * @return vector<Gate> physical circuit that can be executed on hardware,
     *      and modified pi that maps logical qubits to physical qubits.
     */
    vector<Gate> heuristic_search(vector<int> & pi, const vector<vector<int>>& DAG) {
        #ifdef DEBUG
            printf("Initial Mapping:");
            for (int i; i < pi.size(); ++i)
                printf("(%d, %d) ", i, pi[i]);
            printf("\n");
        #endif

        vector<Gate> ans;   // physical circuit
        int tot = 0;        // total number of additional SWAP gates
        
        auto rpi = get_reverse_pi(pi);  // mapping from physical to logical
        vector<double> decay(pi.size(), 1); // decay of logical qubits

        vector<int> indeg(DAG.size(), 0);   // in-degree of DAG nodes
        for (int i = 0; i < DAG.size(); ++i)
            for (int j : DAG[i])
                indeg[j]++;
        
        list<int> F;    // front layer
        for (int i = 0; i < DAG.size(); ++i)
            if (indeg[i] == 0)
                F.push_back(i);
        
        while (!F.empty()) {
            vector<int> executable_list;
            // find all executable gates in F
            for (auto it = F.begin(); it != F.end(); ++it) {
                if (is_executable(pi, *it)) {
                    executable_list.push_back(*it);
                }
            }

            if (!executable_list.empty()) { // execute all executable gates
                for (auto it = F.begin(); it != F.end(); ) {
                    if (is_executable(pi, *it)) {
                        int x = *it;
                        if (gates[x].type == "CNOT") {
                            int p = gates[x].q1;
                            int q = gates[x].q2;
                            double tmp = max(decay[p], decay[q]);
                            decay[p] = decay[q] = tmp + delta2;
                            ans.push_back({"CNOT", pi[p], pi[q], gates[x].tag});
                        }
                        else {
                            int p = gates[x].q1;
                            decay[p] += delta1;
                            ans.push_back({gates[x].type, pi[p], pi[p], gates[x].tag});
                        }

                        for (int y : DAG[x]) {
                            --indeg[y];
                            if (indeg[y] == 0)
                                F.push_back(y);
                        }
                        it = F.erase(it);
                    }
                    else 
                        ++it;
                }
            }
            else {  // If there is no executable gate, try to SWAP
                auto candidate_SWAPs = obtain_SWAPs(F, pi);
                auto E = get_extended_set(F, DAG, indeg);

                #ifdef DEBUG
                    printf("E: ");
                    for (auto e : E) printf("%d ", e);
                    printf("\n");

                    printf("decay: ");
                    for (int i = 0; i < decay.size(); ++i)
                        printf("(%d, %lf) ", i, decay[i]);
                    printf("\n");

                    printf("SWAPs: \n");
                #endif

                // find the SWAP with minimal H-score
                double min_score = __DBL_MAX__;
                pair<int, int> min_SWAP;
                for (auto SWAP : candidate_SWAPs) {
                    int x = SWAP.first, y = SWAP.second;
                    int p = rpi[x], q = rpi[y];
                    
                    auto tmp = pi;
                    swap(tmp[p], tmp[q]);
                    double score = H(F, E, tmp, {p, q}, decay);
                    if (score < min_score) {
                        min_score = score;
                        min_SWAP = SWAP;
                    }
                    #ifdef DEBUG
                        printf("(%d, %d, %lf)\n", x, y, score);
                    #endif
                }

                int x = min_SWAP.first, y = min_SWAP.second;
                int p = rpi[x], q = rpi[y];
                swap(pi[p], pi[q]);
                swap(rpi[x], rpi[y]);
                ans.push_back({"SWAP", x, y, "SWAP" + to_string(++tot)});

                double tmp = max(decay[p], decay[q]);
                decay[p] = decay[q] = tmp + delta2 * 3;
            }
        }

        #ifdef DEBUG
            printf("additional SWAP: %d\n", tot);
            printf("updated Circuit:\n");
            for (auto g : ans) {
                cout << g.type;
                printf(" %d %d ", g.q1, g.q2);
                cout << g.tag << endl;
            }
        #endif

        return ans;
    }

    /**
     * @brief one-turn iterate to update Initial Mapping.
     *      
     * @param pi Input initial mapping
     */
    void iter_one_turn(vector<int> & pi) {
        #ifdef DEBUG
            cout << "Original Circuit\n";
        #endif
        heuristic_search(pi, this->DAG);    // using original circuit to update
        #ifdef DEBUG
            cout << endl;
            cout << "Reversed Circuit\n";
        #endif
        heuristic_search(pi, this->RDAG);   // using reversed circuit to update
    }

    /**
     * @brief solve qubit mapping problem
     * 
     * @param iter_num iterate times to update random initial mapping
     * @param W parameter to look-ahead
     * @param delta1 decay of single gate
     * @param delta2 decay of CNOT gate, decay of SWAP will be 3*delta2
     * @return pair<vector<Gate>, pair<vector<int>, vector<int>>> 
     *      (gs, (pi0, pi1)), gs is generated physical circuit, 
     *                        pi0 is initial mapping from logical to physical
     *                        pi1 is final mapping from logical to physical
     */
    pair<vector<Gate>, pair<vector<int>, vector<int>>> solve(int iter_num, double W, double delta1, double delta2) {
        this->set_parameters(W, delta1, delta2);

        // generate random initial mapping
        vector<int> pi(this->num_physical);
        for (int i = 0; i < pi.size(); ++i)
            pi[i] = i;
        shuffle(pi.begin(), pi.end(), engine);

        // iterate to update initial mapping
        for (int t = 0; t < iter_num; ++t)
        {
            #ifdef DEBUG
                printf("\n\n\n");
                printf("iterate %d:\n", t);
            #endif
            iter_one_turn(pi);
        }
        auto initial_mapping = pi;
        auto gs = heuristic_search(pi, this->DAG);
        return {gs, {initial_mapping, pi}};
    }

    inline void set_parameters(double W, double delta1, double delta2) {
        this->W = W;
        this->delta1 = delta1;
        this->delta2 = delta2;
    }
};

pair<int, vector<vector<int>>> get_physical_graph(string filename) {
    ifstream fin;
    fin.open(filename);

    int n, m;
    fin >> n >> m;

    vector<vector<int>> G(n);
    while (m--) {
        int x, y;
        fin >> x >> y;
        G[x].push_back(y);
        G[y].push_back(x);
    }

    fin.close();
    return {n, G};
}


class Benchmark {
private:
    vector<Gate> gates;
public:
    /**
     * @brief Construct a new Benchmark object
     * 
     * @param gs physical circuit from left to right
     */
    Benchmark(const vector<Gate>& gs) : gates(gs) {

    }

    /**
     * @brief count number of SWAP used
     * 
     * @return int 
     */
    int num_SWAP() const {
        int tot = 0;
        for (auto g : gates)
            if (g.type == "SWAP")
                ++tot;
        return tot;
    }

    /**
     * @brief calculate circuit depth without considering single gates
     * 
     * @return int 
     */
    int circuit_depth() const {
        vector<Gate> gs;    // list of gates without single gates
        int id = 0;
        for (auto g : gates)
            if (g.type == "CNOT")
                gs.push_back(g);
            else if (g.type == "SWAP") {
                gs.push_back({"CNOT", g.q1, g.q2, "s" + to_string(id++)});
                gs.push_back({"CNOT", g.q2, g.q1, "s" + to_string(id++)});
                gs.push_back({"CNOT", g.q1, g.q2, "s" + to_string(id++)});
            }
        
        int n = 0;  // number of qubits
        for (auto g : gs)
            n = max(n, 1 + max(g.q1, g.q2));
        int m = gs.size();  // number of gates, number of DAG nodes
        auto DAG = get_circuit_DAG(n, gs);

        // show_DAG(DAG, gs);

        vector<int> indeg(m, 0);
        for (int i = 0; i < m; ++i)
            for (int j : DAG[i])
                indeg[j]++;
        
        vector<int> F;
        for (int i = 0; i < m; ++i)
            if (indeg[i] == 0)
                F.push_back(i);
        
        int dep = 0;
        while (!F.empty()) {
            ++dep;
            vector<int> E;
            for (int x : F)
                for (int y : DAG[x]) {
                    if (--indeg[y] == 0)
                        E.push_back(y);
                }
            F = E;
        }
        return dep;
    }

    /**
     * @brief calculate total circuit time delay
     * 
     * @param d1 single gate delay
     * @param d2 CNOT gate delay, SWAP is 3*d2
     * @return double 
     */
    double time_delay(double d1, double d2) const {
        int n = 0;  // number of qubits
        for (auto g : gates)
            n = max(n, 1 + max(g.q1, g.q2));
        
        vector<double> delay(n, 0); // delay on each qubit
        for (auto g : gates) {
            if (g.type == "CNOT" || g.type == "SWAP") {
                double tmp = max(delay[g.q1], delay[g.q2]);
                delay[g.q1] = delay[g.q2] = tmp + (g.type == "CNOT" ? d2 : 3.0 * d2);
            }
            else {  // single gate
                delay[g.q1] += d1;
            }
        }
        double ans = 0; // circuit delay
        for (auto d : delay) ans = max(ans, d);
        return ans;
    }
};

int main() {
    string f1 = "logical_circuit.txt";
    string f2 = "coupling_graph.txt";

    double W = 0.5;
    double delta1 = 0.1;
    double delta2 = 0.3;

    auto par = get_logical_circuit_from_file(f1);
    auto num_logical = par.first;
    auto gates = par.second;

    // auto DAG = get_circuit_DAG(num_logical, gates);
    // show_DAG(DAG, gates);

    auto tmp = get_physical_graph(f2);
    int num_physical = tmp.first;
    auto G = tmp.second;

    auto sol = SABRE(num_logical, gates, num_physical, G);

    // {
    //     sol.set_parameters(W, d1, d2);
    //     vector<int> pi = {0, 1, 2, 3};
    //     sol.iter_one_turn(pi);
    // }

    auto ans = sol.solve(5, W, delta1, delta2);
    // #ifdef DEBUG
    // {
        // printf("\n\n");
        // printf("The Solution is:\n");
        auto gs = ans.first;
        auto pi0 = ans.second.first;
        auto pi1 = ans.second.second;
        // printf("Initial Mapping: ");
        // for (int i = 0; i < pi0.size(); ++i)
        //     printf("(%d, %d) ", i, pi0[i]);
        // printf("\n");
        for (auto g : gs) {
            cout << g.type;
            printf(" %d %d ", g.q1, g.q2);
            cout << g.tag << endl;
        }
        printf("Final Mapping:   ");
        for (int i = 0; i < pi1.size(); ++i)
            printf("(%d, %d) ", i, pi1[i]);
        printf("\n");
    // }
    // #endif

    #ifdef BENCHMARK
    {
        const int T = 10;
        for (int t = 0; t < T; ++t) {
            printf("Solution %d:\n", t);
            auto ans = sol.solve(5, W, delta1, delta2);
            auto gs = ans.first;
            Benchmark bench(gs);
            double d1 = 1.0;
            double d2 = 10.0;
            printf("SWAP number: %d\n", bench.num_SWAP());
            printf("circuit depth: %d\n", bench.circuit_depth());
            printf("time delay: %lf\n", bench.time_delay(d1, d2));
            printf("\n");
        }
    }
    #endif

    return 0;
}