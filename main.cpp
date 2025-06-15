/*
Alan Sanmiguel Garay, Juan Diego Susunaga, Adrián Salazar Rodríguez & Manuel Alejandro Cruz
    15/06/2025
    Equipo 2 - Análisis y diseño de algoritmos avanzados (Gpo 602)
    Evidencia Integradora 2
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <climits>
#include <cmath>
#include <string>
#include <sstream>

// Estructura para representar una arista con peso entre dos nodos
struct Edge {
    int weight, u, v;
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

// Estructura para representar un punto en 2D
struct Point {
    int x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
    double Distance(const Point& other) const {
        return std::sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
    }
};

class NetworkOptimizer {
private:
    int n;
    std::vector<std::vector<int>> DistanceMatrix;
    std::vector<std::vector<int>> capacityMatrix;
    std::vector<Point> centrals;
    std::vector<int> parent, rank_;

    void MakeSet(int v) {
        parent[v] = v;
        rank_[v] = 0;
    }

    int FindSet(int v) {
        if (v == parent[v])
            return v;
        return parent[v] = FindSet(parent[v]);
    }

    bool UnionSets(int a, int b) {
        a = FindSet(a);
        b = FindSet(b);
        if (a != b) {
            if (rank_[a] < rank_[b])
                std::swap(a, b);
            parent[b] = a;
            if (rank_[a] == rank_[b])
                rank_[a]++;
            return true;
        }
        return false;
    }

    bool BreadthFirstSearch(std::vector<std::vector<int>>& graph, int s, int t, std::vector<int>& parent) {
        std::vector<bool> visited(n, false);
        std::queue<int> q;
        q.push(s);
        visited[s] = true;
        parent[s] = -1;

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (int v = 0; v < n; v++) {
                if (!visited[v] && graph[u][v] > 0) {
                    if (v == t) {
                        parent[v] = u;
                        return true;
                    }
                    q.push(v);
                    parent[v] = u;
                    visited[v] = true;
                }
            }
        }
        return false;
    }

    // Función auxiliar para inicializar la tabla DP
    void InitializeDPTable(std::vector<std::vector<int>>& dp, std::vector<std::vector<int>>& parent) {
        dp.resize(1 << n, std::vector<int>(n, INT_MAX));
        parent.resize(1 << n, std::vector<int>(n, -1));
        dp[1][0] = 0; // Comenzar en el nodo 0
    }

    // Función auxiliar para procesar cada máscara en el DP
    void ProcessMask(int mask, std::vector<std::vector<int>>& dp, std::vector<std::vector<int>>& parent) {
        for (int u = 0; u < n; u++) {
            if (!(mask & (1 << u)) || dp[mask][u] == INT_MAX)
                continue;

            ProcessNode(mask, u, dp, parent);
        }
    }

    // Función auxiliar para procesar cada nodo en una máscara
    void ProcessNode(int mask, int u, std::vector<std::vector<int>>& dp, std::vector<std::vector<int>>& parent) {
        for (int v = 0; v < n; v++) {
            if ((mask & (1 << v)) || DistanceMatrix[u][v] == 0)
                continue;

            UpdateDPTable(mask, u, v, dp, parent);
        }
    }

    // Función auxiliar para actualizar la tabla DP
    void UpdateDPTable(int mask, int u, int v, std::vector<std::vector<int>>& dp, std::vector<std::vector<int>>& parent) {
        int newMask = mask | (1 << v);
        int newDistance = dp[mask][u] + DistanceMatrix[u][v];

        if (newDistance < dp[newMask][v]) {
            dp[newMask][v] = newDistance;
            parent[newMask][v] = u;
        }
    }

    // Función auxiliar para encontrar el costo mínimo y último nodo
    std::pair<int, int> FindMinCostAndLastNode(const std::vector<std::vector<int>>& dp) {
        int finalMask = (1 << n) - 1;
        int minCost = INT_MAX;
        int lastNode = -1;

        for (int i = 1; i < n; i++) {
            if (DistanceMatrix[i][0] > 0 && dp[finalMask][i] != INT_MAX) {
                if (dp[finalMask][i] + DistanceMatrix[i][0] < minCost) {
                    minCost = dp[finalMask][i] + DistanceMatrix[i][0];
                    lastNode = i;
                }
            }
        }

        return {minCost, lastNode};
    }

    // Función auxiliar para reconstruir el camino
    void ReconstructPath(int lastNode, const std::vector<std::vector<int>>& parent, std::vector<int>& path) {
        std::vector<int> tempPath;
        int mask = (1 << n) - 1;
        int curr = lastNode;

        while (curr != -1) {
            tempPath.push_back(curr);
            int prev = parent[mask][curr];
            mask ^= (1 << curr);
            curr = prev;
        }

        std::reverse(tempPath.begin(), tempPath.end());
        tempPath.push_back(0); // Regresar al inicio
        path = tempPath;
    }

    int TSPDynamicProgramming(std::vector<int>& path) {
        if (n > 15) {
            return TSPNearestNeighbor(path);
        }

        std::vector<std::vector<int>> dp, parentTable;
        InitializeDPTable(dp, parentTable);

        for (int mask = 0; mask < (1 << n); mask++) {
            ProcessMask(mask, dp, parentTable);
        }

        auto [minCost, lastNode] = FindMinCostAndLastNode(dp);

        if (lastNode == -1) {
            return TSPNearestNeighbor(path);
        }

        ReconstructPath(lastNode, parentTable, path);
        return minCost;
    }

    int TSPNearestNeighbor(std::vector<int>& path) {
        std::vector<bool> visited(n, false);
        path.clear();
        path.push_back(0);
        visited[0] = true;
        int totalCost = 0;
        int current = 0;

        for (int i = 0; i < n - 1; i++) {
            int next = -1;
            int minDist = INT_MAX;

            for (int j = 0; j < n; j++) {
                if (!visited[j] && DistanceMatrix[current][j] > 0 && DistanceMatrix[current][j] < minDist) {
                    minDist = DistanceMatrix[current][j];
                    next = j;
                }
            }

            if (next != -1) {
                path.push_back(next);
                visited[next] = true;
                totalCost += minDist;
                current = next;
            }
        }

        path.push_back(0);
        if (DistanceMatrix[current][0] > 0) {
            totalCost += DistanceMatrix[current][0];
        }

        return totalCost;
    }

public:
    void ReadInput() {
        std::cin >> n;

        DistanceMatrix.resize(n, std::vector<int>(n));
        capacityMatrix.resize(n, std::vector<int>(n));
        parent.resize(n);
        rank_.resize(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cin >> DistanceMatrix[i][j];
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cin >> capacityMatrix[i][j];
            }
        }

        centrals.resize(n);
        for (int i = 0; i < n; i++) {
            std::string line;
            std::cin >> line;
            line = line.substr(1, line.length() - 2);
            size_t comma = line.find(',');
            int x = std::stoi(line.substr(0, comma));
            int y = std::stoi(line.substr(comma + 1));
            centrals[i] = Point(x, y);
        }
    }

    std::vector<std::pair<char, char>> kruskalMST() {
        std::vector<Edge> edges;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (DistanceMatrix[i][j] > 0) {
                    edges.push_back({DistanceMatrix[i][j], i, j});
                }
            }
        }

        std::sort(edges.begin(), edges.end());

        for (int i = 0; i < n; i++) {
            MakeSet(i);
        }

        std::vector<std::pair<char, char>> mstEdges;

        for (const Edge& e : edges) {
            if (UnionSets(e.u, e.v)) {
                char u = 'A' + e.u;
                char v = 'A' + e.v;
                mstEdges.push_back({u, v});
                if (mstEdges.size() == n - 1)
                    break;
            }
        }

        return mstEdges;
    }

    std::vector<char> SolveTravellingSalesman() {
        std::vector<int> path;
        TSPNearestNeighbor(path);

        std::vector<char> result;
        for (int node : path) {
            result.push_back('A' + node);
        }

        return result;
    }

    int FordFulkerson() {
        std::vector<std::vector<int>> graph = capacityMatrix;
        std::vector<int> parent(n);
        int maxFlow = 0;
        int source = 0, sink = n - 1;

        while (BreadthFirstSearch(graph, source, sink, parent)) {
            int pathFlow = INT_MAX;

            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                pathFlow = std::min(pathFlow, graph[u][v]);
            }

            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                graph[u][v] -= pathFlow;
                graph[v][u] += pathFlow;
            }

            maxFlow += pathFlow;
        }

        return maxFlow;
    }

    Point findClosestCentral(const Point& query) {
        double minDistance = 1e9;
        Point closest;

        for (const Point& central : centrals) {
            double dist = query.Distance(central);
            if (dist < minDistance) {
                minDistance = dist;
                closest = central;
            }
        }

        return closest;
    }
};

int main() {
    NetworkOptimizer optimizer;
    optimizer.ReadInput();

    // 1. Kruskal MST
    std::vector<std::pair<char, char>> mstEdges = optimizer.kruskalMST();
    std::cout << "1.\n";
    for (size_t i = 0; i < mstEdges.size(); i++) {
        if (i > 0)
            std::cout << "\n";
        std::cout << "(" << mstEdges[i].first << ", " << mstEdges[i].second << ")";
    }
    std::cout << std::endl;

    // 2. TSP
    std::vector<char> tspRoute = optimizer.SolveTravellingSalesman();
    std::cout << "2.\n";
    for (size_t i = 0; i < tspRoute.size(); i++) {
        if (i > 0)
            std::cout << " ";
        std::cout << tspRoute[i];
    }
    std::cout << std::endl;

    // 3. Ford-Fulkerson
    int maxFlow = optimizer.FordFulkerson();
    std::cout << "3.\n" << maxFlow << std::endl;

    // 4. Central más cercana
    std::string queryLine;
    std::cin >> queryLine;
    queryLine = queryLine.substr(1, queryLine.length() - 2);
    size_t comma = queryLine.find(',');
    Point query;
    query.x = std::stoi(queryLine.substr(0, comma));
    query.y = std::stoi(queryLine.substr(comma + 1));

    Point closest = optimizer.findClosestCentral(query);
    std::cout << "4.\n(" << closest.x << ", " << closest.y << ")" << std::endl;

    return 0;
}