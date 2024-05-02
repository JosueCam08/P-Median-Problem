#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <set>
#include <queue>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <random> 
#include <sys/time.h>
#include <climits>
using namespace std;

#define INF INT_MAX

class PMEDIAN{
    public:
        int V; // Number of vertices
        int E; // Number of edges
        int p; // Number of medians

        vector< vector<int> > matAdy;  // Adjacency matrix
        vector< vector<int> > listAdy; // Adjacency list

        vector< vector< vector<int> > > DTree; // Dijkstra trees

        vector<int> best_solution;      // Best solution;
        int best_cost;                  // Best cost
        vector<int> globalBestSolution; // Global best solution
        int globalBestCost;             // Global best cost

        // List of sorted costs for each vertex in the solution
        vector< multiset< pair<int, int> > > cost_sorted;
        
        // Time variables
        double initialTime, finalTime; 

        // Pertubation variables
        int maxDistance = 0;
        double pertubationFactor = 1.0, eps = 0.5;

        // Constructor
        PMEDIAN(string file_name){
            ifstream instance;
            instance.open(file_name);

            instance >> V >> E >> p;

            matAdy.assign(V + 1, vector<int>(V + 1, INF));
            for(int i = 1; i <= V; i++)
                matAdy[i][i] = 0;

            listAdy.assign(V + 1, vector<int>());
            DTree.assign(V + 1, vector< vector<int> >());

            best_solution.assign(p, 0);

            int u, v, w;
            for(int i = 0; i < E; i++){
                instance >> u >> v >> w;
                matAdy[u][v] = matAdy[v][u] = w;
                listAdy[u].push_back(v);
                listAdy[v].push_back(u);
            }

            for(int i = 1; i <= V; i++)
                Dijkstra(i);

            globalBestCost = INF;

            struct timeval currentTime;
            gettimeofday(&currentTime, NULL);
            initialTime = (double) (currentTime.tv_sec) + (double) (currentTime.tv_usec) / 1.0e6;

            instance.close();

            return;
        }

        // Evaluate the cost of a solution
        int evaluate(vector<int> &solution){
            int cost = 0;
            // Go through each vertex
            for(int i = 1; i <= V; i++){
                // Go through each median
                for(int j = 0; j < p; j++)
                    cost_sorted[i].insert(make_pair(matAdy[i][solution[j]], solution[j]));
                // Get the lowest cost
                auto low_cost = cost_sorted[i].begin();
                if(low_cost->first < INF)
                    cost += low_cost->first;
            }
            return cost;
        }

        // Stochastic Hill Climbing
        void LocalSearchOneNeighbor(int seed){
            // Get the initial solution
            InitializeSolution();

            // Neighborhood
            vector< vector<int> > N;

            // Random number generator
            auto rng = default_random_engine(seed);

            bool flag = true;
            int cost;
            while(flag){
                GenerateNeighborhood(N);
                shuffle(N.begin(), N.end(), rng);

                flag = false;
                for(int i = 0; i < N.size(); i++){
                    cost = IncrementalEvaluation(N[i]);
                    if(cost < best_cost){
                        best_cost = cost;
                        UpdateCostSorted(N[i]);
                        best_solution = N[i];
                        flag = true;
                        break;
                    }
                }
            }

            return;
        }

        // Local search with two neighborhoods
        void LocalSearchTwoNeighboor(char op = 'r'){
            // Get the initial solution
            InitializeSolution(op);

            bool flag = true;
            while(flag){
                flag = false;

                // Generate candidate changes
                int ch1, ch2;
                GenerateChanges(ch1, ch2);

                // Generate table of improvemts
                vector< pair<int, int> > ImprovmentsTable1, ImprovmentsTable2;
                vector<int> ign; ign.push_back(ch1); ign.push_back(ch2);
                CreateImprovementsTable(ImprovmentsTable1, ch1, ign);
                CreateImprovementsTable(ImprovmentsTable2, ch2, ign);

                // Sort the table of improvements
                sort(ImprovmentsTable1.begin(), ImprovmentsTable1.end());
                sort(ImprovmentsTable2.begin(), ImprovmentsTable2.end());

                // Get the best solution comparing the table of improvements
                int best_improve = 0, break_flag = 0, current_improve = 0;
                int best_ch1, best_ch2, current_cost, delete_cost;

                // Get the cost of the current best solution without
                // consider the median at the ign indexs
                delete_cost = IncrementalEvaluationDelete(ign);

                // Testing 
                // if(delete_cost != DeleteCostBruteForce(ch1, ch2))
                //    cout << "WARNING: delete_cost no coincide (" << delete_cost << " != " << DeleteCostBruteForce(ch1, ch2) << ")" << '\n';

                // Go through the tables of improvements and get the best neighbors
                // to make the change and improve the solution
                for(int i = 0; i < ImprovmentsTable1.size(); i++){
                    break_flag = 0;
                    for(int j = 0; j < ImprovmentsTable2.size(); j++){
                        if(abs(ImprovmentsTable1[i].first) + abs(ImprovmentsTable2[j].first) < best_improve)
                            break;

                        vector<int> add_nodes;
                        add_nodes.push_back(ImprovmentsTable1[i].second);
                        add_nodes.push_back(ImprovmentsTable2[j].second);

                        current_cost = IncrementalEvaluationMerge(delete_cost, add_nodes, ign);
                        current_improve = delete_cost - current_cost;

                        if(current_improve > best_improve){
                            best_improve = current_improve;
                            best_ch1 = ImprovmentsTable1[i].second;
                            best_ch2 = ImprovmentsTable2[j].second;
                        }

                        break_flag++;
                    }
                    if(break_flag == 0) break;
                }

                // Compare the best solution with the current solution
                vector<int> add_nodes; add_nodes.push_back(best_ch1); add_nodes.push_back(best_ch2);
                int candidate_cost = IncrementalEvaluationMerge(delete_cost, add_nodes, ign);

                if(candidate_cost < best_cost){

                    best_cost = candidate_cost;
                    vector<int> candidate_solution = best_solution;
                    candidate_solution[ch1] = best_ch1;
                    candidate_solution[ch2] = best_ch2;

                    UpdateCostSorted(candidate_solution);
                    best_solution = candidate_solution;
                    flag = true;
                    
                    // Testing
                    // if(BestCostBruteForce() != best_cost)
                    //    cout << "WARNING: best_cost no coincide (" << best_cost << " != " << BestCostBruteForce() << ")" << '\n';
                }
            }
            return;
        }

        // Local search with hard reset
        void LocalSearchHardReset(ofstream &out, int seed, double finalTime){
            double currentTime, elapsedTime, rangeTime = 60.0, lastPrintTime = 0.0;
            this->finalTime = finalTime;

            do{
                cost_sorted.assign(V + 1, multiset< pair<int, int> >());
                LocalSearchTwoNeighboor();
                globalBestCost = min(globalBestCost, best_cost);

                // Get the current time
                struct timeval cTime;
                gettimeofday(&cTime, NULL);
                currentTime = (double) (cTime.tv_sec) + (double) (cTime.tv_usec) / 1.0e6;
                elapsedTime = currentTime - initialTime;

                // Print the best solution each rangeTime seconds
                if((elapsedTime - lastPrintTime) > rangeTime && elapsedTime < finalTime){
                    out << globalBestCost << '\n';
                    lastPrintTime = elapsedTime;
                }

                cost_sorted.clear();
            } while(elapsedTime < finalTime);

            return;
        }

        // Local search with iterated local search
        void ILS(ofstream &out, int seed, double finalTime){
            double currentTime, elapsedTime, rangeTime = 60.0, lastPrintTime = 0.0;
            this->finalTime = finalTime;

            // Initialize the solution with a random method
            cost_sorted.assign(V + 1, multiset< pair<int, int> >());
            LocalSearchTwoNeighboor();
            globalBestCost = best_cost;
            globalBestSolution = best_solution;

            int cont = 0;
            do{
                // Clean the structure to store the costs
                cost_sorted.clear();
                cost_sorted.assign(V + 1, multiset< pair<int, int> >());

                // Pertubate the solution and apply the local search
                LocalSearchTwoNeighboor('p');

                // Accept the new solution if it is better
                if(best_cost < globalBestCost){
                    globalBestCost = best_cost;
                    globalBestSolution = best_solution;
                }

                // Get the current time
                struct timeval cTime;
                gettimeofday(&cTime, NULL);
                currentTime = (double) (cTime.tv_sec) + (double) (cTime.tv_usec) / 1.0e6;
                elapsedTime = currentTime - initialTime;

                // Print the best solution each rangeTime seconds
                if((elapsedTime - lastPrintTime) > rangeTime && elapsedTime < finalTime){
                    out << globalBestCost << '\n';
                    lastPrintTime = elapsedTime;
                }

                // Update the pertubation factor
                pertubationFactor = 1.0 - (elapsedTime / finalTime) + eps;

            } while(elapsedTime < finalTime);

            return;
        }

        // erase this
        int BestCostBruteForce(){
            int ccost = 0;
            for(int i = 1; i <= V; i++)
                ccost += cost_sorted[i].begin()->first;

            return ccost;
        }
        
        // erase this
        int DeleteCostBruteForce(int ch1, int ch2){
            int cc = 0;
            for(int i = 1; i <= V; i++){
                int c = INF;
                for(int j = 0; j < p; j++){
                    if(j == ch1 || j == ch2) continue;
                    c = min(c, matAdy[i][best_solution[j]]);
                }
                cc += c;
            }       
            return cc;
        }

        // erase this
        int CostBestSolutionBruteForce(vector< pair<int, int> > &ImprovmentsTable1,
                                       vector< pair<int, int> > &ImprovmentsTable2,
                                       int ch1, int ch2){
            int best_cc = INF;
            for(int i = 0; i < ImprovmentsTable1.size(); i++){
                for(int j = 0; j < ImprovmentsTable2.size(); j++){
                    vector<int> force_solution = best_solution;
                    force_solution[ch1] = ImprovmentsTable1[i].second;
                    force_solution[ch2] = ImprovmentsTable2[j].second;

                    int cc = 0;
                    for(int k = 1; k <= V; k++){
                        int c = INF;
                        for(int kk = 0; kk < p; kk++)
                            c = min(c, matAdy[k][force_solution[kk]]);
                        cc += c;
                    }

                    best_cc = min(best_cc, cc);
                }
            }
            return best_cc;
        }

        void PrintSolution(vector<int> &solution){
            for(int i = 0; i < p; i++)
                cout << solution[i] << " ";
            cout << '\n';
            return;
        }

    private:

        /* > > > This section is for general useful funcions < < < */

        // Dijkstra algorithm
        void Dijkstra(int node){
            // Initialize the distance and parent vectors
            vector<int> dist(V + 1, INF);
            vector<int> parent(V + 1, -1);
            
            // Set the distance of the node to 0 and insert it into the queue
            dist[node] = 0;
            set< pair<int, int> > q;
            q.insert(make_pair(0, node));

            // Iterate through the graph
            while(!q.empty()){
                int v = q.begin()->second;
                q.erase(q.begin());

                for(auto edge : listAdy[v]){
                    int to = edge;
                    int len = matAdy[v][to];

                    if(dist[v] + len < dist[to]){
                        q.erase(make_pair(dist[to], to));
                        dist[to] = dist[v] + len;
                        parent[to] = v;
                        q.insert(make_pair(dist[to], to));
                    }
                }
            }

            // Update the adjacency matrix with the shortest paths distances
            for(int i = 1; i <= V; i++)
                matAdy[node][i] = dist[i];

            // Create the Dijkstra tree
            CreateDijkstraTree(node, parent);

            return;
        }
        
        // Solution of the assignment problem
        int HungarianAlgorithm(int n, int m, vector< vector<int> > &A){
            vector<int> u (n+1), v (m+1), p (m+1), way (m+1);
            for (int i=1; i<=n; ++i) {
                p[0] = i;
                int j0 = 0;
                vector<int> minv (m+1, INF);
                vector<bool> used (m+1, false);
                do {
                    used[j0] = true;
                    int i0 = p[j0],  delta = INF,  j1;
                    for (int j=1; j<=m; ++j)
                        if (!used[j]) {
                            int cur = A[i0][j]-u[i0]-v[j];
                            if (cur < minv[j])
                                minv[j] = cur,  way[j] = j0;
                            if (minv[j] < delta)
                                delta = minv[j],  j1 = j;
                        }
                    for (int j=0; j<=m; ++j)
                        if (used[j])
                            u[p[j]] += delta,  v[j] -= delta;
                        else
                            minv[j] -= delta;
                    j0 = j1;
                } while (p[j0] != 0);
                do {
                    int j1 = way[j0];
                    p[j0] = p[j1];
                    j0 = j1;
                } while (j0);
            }

            vector<int> ans (n+1);
            for (int j=1; j<=m; ++j)
                ans[p[j]] = j;

            return -v[0];
        }

        // Create the Dijkstra tree for a node
        void CreateDijkstraTree(int node, vector<int> &parent){
            DTree[node].assign(V + 1, vector<int>());

            for(int i = 1; i <= V; i++){
                if(parent[i] == -1) continue;
                DTree[node][parent[i]].push_back(i);
            }

            return;
        }

        // Create a random solution
        void CreateRandomSolution(vector<int> &solution){
            unordered_set<int> medians;

            while(medians.size() < p){
                int median = (rand() % V) + 1;
                medians.insert(median);
            }

            int it = 0;
            for(auto v : medians)
                solution[it++] = v;

            return;
        }

        // Create a pertubated solution considering the pertubation factor
        void CreatePertubatedSolution(vector<int> &solution){
            // int d;
            // float distancePertubation;
            
            //do{
            CreateRandomSolution(solution);
            for(int i = 0; i < p; i++)
                if(((float)rand() / (float) RAND_MAX) > pertubationFactor)
                    solution[i] = globalBestSolution[i];

            // vector< vector<int> > A(p + 1, vector<int> (p + 1, 0));
            // for(int i = 1; i <= p; i++)
            //     for(int j = 1; j <= p; j++)
            //         A[i][j] = matAdy[globalBestSolution[i - 1]][solution[j - 1]];

            // d = HungarianAlgorithm(p, p, A);
            // maxDistance = max(maxDistance, d);

            // distancePertubation = (float) d / (float) maxDistance;
            //} while(distancePertubation > pertubationFactor);

            // cout << "Distancia: " << distancePertubation << " Perturbacion: " << pertubationFactor << '\n';
            return;
        }

        // Initialize the solution with a constructive or random method
        void InitializeSolution(char method = 'r'){
            if(method == 'r')
                CreateRandomSolution(best_solution);
            else
                CreatePertubatedSolution(best_solution);
            best_cost = evaluate(best_solution);
            return;
        }

        // Update the cost_sorted list with the new solution, erasing the old median and adding the new one
        void UpdateCostSorted(vector<int> &solution){
            for(int i = 0; i < p; i++){
                if(best_solution[i] != solution[i]){
                    for(int j = 1; j <= V; j++){
                        cost_sorted[j].erase(
                            cost_sorted[j].find(make_pair(matAdy[j][best_solution[i]], best_solution[i]))
                            );
                        cost_sorted[j].insert(make_pair(matAdy[j][solution[i]], solution[i]));
                    }
                }
            }
            return;
        }


        /* > > > This section is for one neighborhood funcions < < < */

        // Generate the neighborhood of the current solution selecting a random median
        // and replacing it with a neighbor of the vertex
        void GenerateNeighborhood(vector< vector<int> > &N){
            N.clear();
            int k = rand() % p;
            int node = best_solution[k];

            for(int i = 0; i < listAdy[node].size(); i++){
                vector<int> neighbor(best_solution);
                neighbor[k] = listAdy[node][i];

                N.push_back(neighbor);
            }

            return;
        }

        // Retulize the cost of the best solution to calculate the cost of a new solution
        // with help of the cost_sorted list
        int IncrementalEvaluation(vector<int> &solution){
            int cost = best_cost;

            for(int i = 0; i < p; i++){
                if(best_solution[i] != solution[i]){
                    queue<int> q;

                    // First we are to calculate the cost of the current solution
                    // without the median that will be removed
                    q.push(best_solution[i]);

                    // Add the cost of the median that will be removed
                    auto low_cost_best_solution = cost_sorted[best_solution[i]].begin();
                    low_cost_best_solution++;
                    cost += low_cost_best_solution->first;

                    while(!q.empty()){
                        int current_node = q.front();
                        q.pop();

                        // Go through the Dijsktra tree of the median that will be removed
                        // and update the cost of the nodes that are linked to the median
                        for(auto next : DTree[best_solution[i]][current_node]){
                            auto low_cost = cost_sorted[next].begin();
                            if(low_cost->second != best_solution[i]) continue;

                            cost -= low_cost->first;
                            low_cost++;
                            cost += low_cost->first;
                            q.push(next);
                        }
                    }

                    // Now we are to calculate the cost consider the candidate median
                    // in the current solution without the median that will be removed
                    q.push(solution[i]);

                    // Erase the cost of the median that will be inserted
                    auto low_cost_solution = cost_sorted[solution[i]].begin();
                    if(low_cost_solution->second == best_solution[i]) low_cost_solution++;
                    cost -= low_cost_solution->first;

                    while(!q.empty()){
                        int current_node = q.front();
                        q.pop();
                        
                        // Go through the Dijkstra tree of the candidate median and update the cost
                        for(auto next : DTree[solution[i]][current_node]){
                            auto low_cost = cost_sorted[next].begin();

                            if(low_cost->second == best_solution[i])
                                low_cost++;

                            if(low_cost->first > matAdy[next][solution[i]]){
                                cost -= low_cost->first;
                                cost += matAdy[next][solution[i]];
                                q.push(next);
                            }
                        }
                    }
                }
            }

            return cost;
        }
        

        /* > > > This section is for two neighborhood funcions < < < */

        // Select two random medians to generate a change
        void GenerateChanges(int &ch1, int &ch2){
            ch1 = rand() % p;
            ch2 = rand() % p;

            while(ch1 == ch2) ch2 = rand() % p;

            return;
        }
        
        // Create the table of improvements for the neighbords of the median at the index ch
        // considering that the medians at the indexes in the vector ign are ignored
        void CreateImprovementsTable(vector< pair<int, int> > &ImprovmentsTable, int ch, vector<int> &ign){
            // Fill the table with the neighbors
            for(auto neighboor : listAdy[best_solution[ch]])
                ImprovmentsTable.push_back(make_pair(INF, neighboor));

            // Get the cost of the current solution without consider the medians at ign positions
            int cost_ign = IncrementalEvaluationDelete(ign);

            // Get the cost of improvment for each neighboor
            for(int i = 0; i < ImprovmentsTable.size(); i++){
                // Get the cost of the current solution without consider the medians at ign positions
                // and adding a neighboor as a new median
                int cost = IncrementalEvaluationAdd(cost_ign, ImprovmentsTable[i].second, ign);
                ImprovmentsTable[i].first = cost - cost_ign;
            }

            return;
        }

        // Get the cost of the current best solution without consider the median at the ign indexs
        // and adding the neighboors ch1 and ch2 as new medians
        int IncrementalEvaluationMerge(int cost, vector<int> add_nodes, vector<int> &ign){
            queue<int> q;
            map<int, int> min_cost;

            // Go through the nodes that will be added as medians
            for(int i = 0; i < add_nodes.size(); i++){
                int node = add_nodes[i];
                q.push(node);

                // Erase the cost of the median that will be inserted
                auto low_cost_node1 = cost_sorted[node].begin();
                SkipIgnoredNodes(low_cost_node1, ign);
                min_cost[node] = min(min_cost[node], -low_cost_node1->first);

                while(!q.empty()){
                    int current_node = q.front();
                    q.pop();

                    // Go through the Dijkstra tree of the candidate median and update the cost
                    for(auto next : DTree[node][current_node]){
                        auto low_cost = cost_sorted[next].begin();
                        
                        // Get the lowest cost ignoring the medians at the ign indexs
                        SkipIgnoredNodes(low_cost, ign);
                        if(matAdy[next][node] > low_cost->first) continue;

                        // Add the cost of the candidate median
                        min_cost[next] = min(min_cost[next], matAdy[next][node] - low_cost->first);
                        q.push(next);
                    }
                }
            }

            // Calculate the cost of adding the nodes
            for(auto node : min_cost)
                cost += node.second;

            return cost;
        }

        // Return the cost of add add_node for the solution corresponding to cost
        int IncrementalEvaluationAdd(int cost, int add_node, vector<int> &ign){
            queue<int> q;
            
            // Now we are going to calculate the cost consider the candidate median
            // in the current solution without the median that will be removed
            q.push(add_node);

            // Erase the cost of the median that will be inserted
            auto low_cost_add_node = cost_sorted[add_node].begin();
            SkipIgnoredNodes(low_cost_add_node, ign);
            cost -= low_cost_add_node->first;

            while(!q.empty()){
                int current_node = q.front();
                q.pop();
                
                // Go through the Dijkstra tree of the candidate median and update the cost
                for(auto next : DTree[add_node][current_node]){
                    auto low_cost = cost_sorted[next].begin();

                    SkipIgnoredNodes(low_cost, ign);

                    if(low_cost->first > matAdy[next][add_node]){
                        cost -= low_cost->first;
                        cost += matAdy[next][add_node];
                        q.push(next);
                    }
                }
            }

            return cost;
        }

        // Return the cost of the current best solution without consider
        // the median at the ign indexs
        int IncrementalEvaluationDelete(vector<int> &ign){
            int cost = best_cost, delete_node;
            queue<int> q;

            // Go through the nodes that are linked to the median that will be removed
            for(int i = 0; i < ign.size(); i++){

                // First we are going to consider the nodes matches
                // connected to the median best_solution[ign[i]]
                int delete_node = best_solution[ign[i]];
                q.push(delete_node);

                // Add the cost of linked the delete_node to any other median
                auto low_cost_delete_node = cost_sorted[delete_node].begin();
                SkipIgnoredNodes(low_cost_delete_node, ign);
                cost += low_cost_delete_node->first;

                // Now we are going to the Dijkstra Tree of the delete node
                while(!q.empty()){
                    int current_node = q.front();
                    q.pop();

                    // Go through the Dijsktra tree of the median that will be removed
                    // and update the cost of the nodes that are linked to the median
                    for(auto next : DTree[delete_node][current_node]){
                        auto low_cost = cost_sorted[next].begin();

                         if(low_cost->second != delete_node) continue;
                        SkipIgnoredNodes(low_cost, ign);
                        
                        cost -= cost_sorted[next].begin()->first;
                        cost += low_cost->first;
                        q.push(next);
                    }
                }
            }

            return cost;
        }
        
        // Skip the ignored nodes from a set of cost sorted for the current best solution
        void SkipIgnoredNodes(set< pair<int, int> >::iterator &low_cost, vector<int> &ign){
            bool flag = true;
            while(flag){
                flag = false;
                for(int i = 0; i < ign.size(); i++){
                    if(low_cost->second == best_solution[ign[i]]){
                        flag = true;
                        low_cost++;
                        break;
                    }
                }
            }
            return;
        }

};

int main(int argc, char *argv[]){

    // Check if the number of arguments is correct
    if(argc != 5){
        cout << "Error. Usage: " << argv[0] 
             << " instance=<instance_file> seed=<seed> output=<output_file> duration=<time in seconds>\n";
        exit(0);
    }

    // Variables for the arguments
    string instance_file, output_file;
    int seed;
    double finalTime;

    for(int i = 1; i < argc; i++){
        // Variables
        string arg = argv[i];
        int j = 0;

        // Find the '=' character
        while(j < arg.size() && arg[j] != '=')
            j++;
        
        // Check if the argument is valid
        if(j == arg.size()){
            cout << "Error. Argument " << i << " is invalid.\n";
            exit(0);
        }

        // Check the type of argument
        if(arg.substr(0, j) == "instance")
            instance_file = arg.substr(j + 1);
        else if(arg.substr(0, j) == "seed")
            seed = stoi(arg.substr(j + 1));
        else if(arg.substr(0, j) == "output")
            output_file = arg.substr(j + 1);
        else if(arg.substr(0, j) == "duration")
            finalTime = stod(arg.substr(j + 1));
        else{
            cout << "Error. Argument " << i << " is invalid.\n";
            exit(0);
        }
    }

    // Initialize parameters
    srand(seed);
    PMEDIAN *pmed = new PMEDIAN(instance_file);
    ofstream out;
    out.open(output_file);

    // Variables to measure time
    struct timeval start, end;

    gettimeofday(&start, NULL);
    pmed->ILS(out, seed, finalTime);
    gettimeofday(&end, NULL);

    // Calculate the time taken by the local search
    double time_taken = (end.tv_sec - start.tv_sec) * 1e6;
    time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;

    // Print the results
    out << pmed->globalBestCost << " " << time_taken << endl;

    out.close();
    return 0;
}