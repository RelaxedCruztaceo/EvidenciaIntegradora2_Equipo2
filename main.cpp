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
#include <memory>

// Estructura para representar un punto en 2D
struct Point {
    int x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
    double Distance(const Point& other) const {
        return std::sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
    }
};

// Clase base para algoritmos de grafos
class GraphAlgorithm {
protected:
    int n;
    std::vector<std::vector<int>> distanceMatrix;

public:
    GraphAlgorithm(int nodes, const std::vector<std::vector<int>>& matrix)
        : n(nodes), distanceMatrix(matrix) {}
    virtual ~GraphAlgorithm() = default;
};

// Clase para manejar Union-Find (Disjoint Set Union)
class UnionFind {
private:
    std::vector<int> parent, rank;

public:
    UnionFind(int size) {
        parent.resize(size);
        rank.resize(size);
        for (int i = 0; i < size; ++i) {
            MakeSet(i);
        }
    }

    void MakeSet(int v) {
        parent[v] = v;
        rank[v] = 0;
    }

    int FindSet(int v) {
        if (v == parent[v]) return v;
        return parent[v] = FindSet(parent[v]);
    }

    bool UnionSets(int a, int b) {
        a = FindSet(a);
        b = FindSet(b);
        if (a != b) {
            if (rank[a] < rank[b]) std::swap(a, b);
            parent[b] = a;
            if (rank[a] == rank[b]) rank[a]++;
            return true;
        }
        return false;
    }
};

// Clase para Kruskal MST
class KruskalMST : public GraphAlgorithm {
private:
    struct Edge {
        int weight, u, v;
        bool operator<(const Edge& other) const {
            return weight < other.weight;
        }
    };

public:
    KruskalMST(int nodes, const std::vector<std::vector<int>>& matrix)
        : GraphAlgorithm(nodes, matrix) {}

    std::vector<std::pair<char, char>> FindMST() {
        std::vector<Edge> edges;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (distanceMatrix[i][j] > 0) {
                    edges.push_back({distanceMatrix[i][j], i, j});
                }
            }
        }

        std::sort(edges.begin(), edges.end());
        UnionFind uf(n);
        std::vector<std::pair<char, char>> mstEdges;

        for (const Edge& e : edges) {
            if (uf.UnionSets(e.u, e.v)) {
                mstEdges.push_back({'A' + e.u, 'A' + e.v});
                if (mstEdges.size() == n - 1) break;
            }
        }

        return mstEdges;
    }
};

// Clase para TSP
class TSPSolver : public GraphAlgorithm {
private:
    struct TSPResult {
        std::vector<int> path;
        int totalCost;
    };

    TSPResult SolveDynamicProgramming() {
        std::vector<std::vector<int>> dp(1 << n, std::vector<int>(n, INT_MAX));
        std::vector<std::vector<int>> parent(1 << n, std::vector<int>(n, -1));
        dp[1][0] = 0;

        for (int mask = 0; mask < (1 << n); mask++) {
            for (int u = 0; u < n; u++) {
                if (!(mask & (1 << u)) || dp[mask][u] == INT_MAX) continue;

                for (int v = 0; v < n; v++) {
                    if ((mask & (1 << v)) || distanceMatrix[u][v] == 0) continue;

                    int newMask = mask | (1 << v);
                    if (dp[newMask][v] > dp[mask][u] + distanceMatrix[u][v]) {
                        dp[newMask][v] = dp[mask][u] + distanceMatrix[u][v];
                        parent[newMask][v] = u;
                    }
                }
            }
        }

        int finalMask = (1 << n) - 1;
        int minCost = INT_MAX;
        int lastNode = -1;

        for (int i = 1; i < n; i++) {
            if (distanceMatrix[i][0] > 0 && dp[finalMask][i] != INT_MAX) {
                if (dp[finalMask][i] + distanceMatrix[i][0] < minCost) {
                    minCost = dp[finalMask][i] + distanceMatrix[i][0];
                    lastNode = i;
                }
            }
        }

        TSPResult result;
        result.totalCost = minCost;

        if (lastNode != -1) {
            std::vector<int> path;
            int mask = finalMask;
            int curr = lastNode;

            while (curr != -1) {
                path.push_back(curr);
                int prev = parent[mask][curr];
                mask ^= (1 << curr);
                curr = prev;
            }

            std::reverse(path.begin(), path.end());
            path.push_back(0);
            result.path = path;
        }

        return result;
    }

    TSPResult SolveNearestNeighbor() {
        std::vector<bool> visited(n, false);
        std::vector<int> path;
        path.push_back(0);
        visited[0] = true;
        int totalCost = 0;
        int current = 0;

        for (int i = 0; i < n - 1; i++) {
            int next = -1;
            int minDist = INT_MAX;

            for (int j = 0; j < n; j++) {
                if (!visited[j] && distanceMatrix[current][j] > 0 &&
                    distanceMatrix[current][j] < minDist) {
                    minDist = distanceMatrix[current][j];
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
        if (distanceMatrix[current][0] > 0) {
            totalCost += distanceMatrix[current][0];
        }

        return {path, totalCost};
    }

public:
    TSPSolver(int nodes, const std::vector<std::vector<int>>& matrix)
        : GraphAlgorithm(nodes, matrix) {}

    std::vector<char> Solve() {
        TSPResult result;
        if (n > 15) {
            result = SolveNearestNeighbor();
        } else {
            result = SolveDynamicProgramming();
            if (result.path.empty()) {
                result = SolveNearestNeighbor();
            }
        }

        std::vector<char> charPath;
        for (int node : result.path) {
            charPath.push_back('A' + node);
        }
        return charPath;
    }
};

// Clase para Ford-Fulkerson
class MaxFlowSolver : public GraphAlgorithm {
private:
    std::vector<std::vector<int>> capacityMatrix;

    bool BFS(const std::vector<std::vector<int>>& graph, int s, int t,
             std::vector<int>& parent) {
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

public:
    MaxFlowSolver(int nodes, const std::vector<std::vector<int>>& distMatrix,
                  const std::vector<std::vector<int>>& capMatrix)
        : GraphAlgorithm(nodes, distMatrix), capacityMatrix(capMatrix) {}

    int CalculateMaxFlow() {
        std::vector<std::vector<int>> graph = capacityMatrix;
        std::vector<int> parent(n);
        int maxFlow = 0;
        int source = 0, sink = n - 1;

        while (BFS(graph, source, sink, parent)) {
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
};

// Clase para manejar centrales
class CentralManager {
private:
    std::vector<Point> centrals;

public:
    CentralManager(const std::vector<Point>& points) : centrals(points) {}

    Point FindClosestCentral(const Point& query) const {
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

// Clase principal que coordina todo
class NetworkOptimizer {
private:
    int n;
    std::vector<std::vector<int>> distanceMatrix;
    std::vector<std::vector<int>> capacityMatrix;
    std::vector<Point> centrals;

public:
    void ReadInput() {
        std::cin >> n;

        distanceMatrix.resize(n, std::vector<int>(n));
        capacityMatrix.resize(n, std::vector<int>(n));
        centrals.resize(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cin >> distanceMatrix[i][j];
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cin >> capacityMatrix[i][j];
            }
        }

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

    std::vector<std::pair<char, char>> GetMST() {
        KruskalMST mstSolver(n, distanceMatrix);
        return mstSolver.FindMST();
    }

    std::vector<char> SolveTSP() {
        TSPSolver tspSolver(n, distanceMatrix);
        return tspSolver.Solve();
    }

    int CalculateMaxFlow() {
        MaxFlowSolver flowSolver(n, distanceMatrix, capacityMatrix);
        return flowSolver.CalculateMaxFlow();
    }

    Point FindClosestCentral(const Point& query) {
        CentralManager centralMgr(centrals);
        return centralMgr.FindClosestCentral(query);
    }
};

int main() {
    NetworkOptimizer optimizer;
    optimizer.ReadInput();

    // 1. Kruskal MST
    std::vector<std::pair<char, char>> mstEdges = optimizer.GetMST();
    std::cout << "1.\n";
    for (size_t i = 0; i < mstEdges.size(); i++) {
        if (i > 0) std::cout << "\n";
        std::cout << "(" << mstEdges[i].first << ", " << mstEdges[i].second << ")";
    }
    std::cout << std::endl;

    // 2. TSP
    std::vector<char> tspRoute = optimizer.SolveTSP();
    std::cout << "2.\n";
    for (size_t i = 0; i < tspRoute.size(); i++) {
        if (i > 0) std::cout << " ";
        std::cout << tspRoute[i];
    }
    std::cout << std::endl;

    // 3. Ford-Fulkerson
    int maxFlow = optimizer.CalculateMaxFlow();
    std::cout << "3.\n" << maxFlow << std::endl;

    // 4. Central más cercana
    std::string queryLine;
    std::cin >> queryLine;
    queryLine = queryLine.substr(1, queryLine.length() - 2);
    size_t comma = queryLine.find(',');
    Point query;
    query.x = std::stoi(queryLine.substr(0, comma));
    query.y = std::stoi(queryLine.substr(comma + 1));

    Point closest = optimizer.FindClosestCentral(query);
    std::cout << "4.\n(" << closest.x << ", " << closest.y << ")" << std::endl;

    return 0;
}
