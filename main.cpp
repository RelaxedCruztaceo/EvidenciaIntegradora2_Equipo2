/*
Alan Sanmiguel Garay, Juan Diego Susunaga, Adrián Salazar Rodríguez y Manuel Alejandro Cruz
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
    int weight, u, v; // peso, nodo origen, nodo destino
    bool operator<(const Edge& other) const {
        return weight < other.weight; // Para ordenar las aristas por peso
    }
};

// Estructura para representar un punto en 2D
struct Point {
    int x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
    // Calcula la distancia euclidiana entre dos puntos
    double Distance(const Point& other) const {
        return std::sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
    }
};

class NetworkOptimizer {
private:
    int n; // Número de nodos/colonias
    std::vector<std::vector<int>> DistanceMatrix; // Matriz de distancias
    std::vector<std::vector<int>> capacityMatrix; // Matriz de capacidades
    std::vector<Point> centrals; // Coordenadas de las centrales

    // Estructuras para Union-Find (usado en Kruskal)
    std::vector<int> parent, rank_;

    // Inicializa un conjunto para Union-Find - O(1)
    void MakeSet(int v) {
        parent[v] = v;
        rank_[v] = 0;
    }

    // Encuentra el representante de un conjunto con compresión de ruta - O(α(n)) (inverse Ackermann)
    int FindSet(int v) {
        if (v == parent[v])
            return v;
        return parent[v] = FindSet(parent[v]);
    }

    // Une dos conjuntos - O(α(n))
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

    // BFS para Ford-Fulkerson - O(n^2)
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

    // TSP con programación dinámica - O(n^2 * 2^n)
    int TravellingSalesmanDynamicProgramming(std::vector<int>& path) {
        if (n > 15) {
            return TravellingSalesmanNearestNeighbor(path);
        }

        std::vector<std::vector<int>> dp(1 << n, std::vector<int>(n, INT_MAX));
        std::vector<std::vector<int>> parent(1 << n, std::vector<int>(n, -1));

        dp[1][0] = 0; // Comenzar en el nodo 0

        for (int mask = 0; mask < (1 << n); mask++) {
            for (int u = 0; u < n; u++) {
                if (!(mask & (1 << u)) || dp[mask][u] == INT_MAX)
                    continue;

                for (int v = 0; v < n; v++) {
                    if (mask & (1 << v) || DistanceMatrix[u][v] == 0)
                        continue;

                    int newMask = mask | (1 << v);
                    if (dp[newMask][v] > dp[mask][u] + DistanceMatrix[u][v]) {
                        dp[newMask][v] = dp[mask][u] + DistanceMatrix[u][v];
                        parent[newMask][v] = u;
                    }
                }
            }
        }

        // Encontrar el mejor final
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

        if (lastNode == -1) {
            return TravellingSalesmanNearestNeighbor(path);
        }

        // Reconstruir camino
        std::vector<int> tempPath;
        int mask = finalMask;
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
        return minCost;
    }

    // TSP con heurística del vecino más cercano - O(n^2)
    int TravellingSalesmanNearestNeighbor(std::vector<int>& path) {
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

        path.push_back(0); // Regresar al inicio
        if (DistanceMatrix[current][0] > 0) {
            totalCost += DistanceMatrix[current][0];
        }

        return totalCost;
    }

public:
    // Lee la entrada del problema - O(n^2)
    void ReadInput() {
        std::cin >> n;

        DistanceMatrix.resize(n, std::vector<int>(n));
        capacityMatrix.resize(n, std::vector<int>(n));
        parent.resize(n);
        rank_.resize(n);

        // Leer matriz de distancias
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cin >> DistanceMatrix[i][j];
            }
        }

        // Leer matriz de capacidades
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cin >> capacityMatrix[i][j];
            }
        }

        // Leer coordenadas de centrales
        centrals.resize(n);
        for (int i = 0; i < n; i++) {
            std::string line;
            std::cin >> line;
            // Parsear (x,y)
            line = line.substr(1, line.length() - 2); // Quitar paréntesis
            size_t comma = line.find(',');
            int x = std::stoi(line.substr(0, comma));
            int y = std::stoi(line.substr(comma + 1));
            centrals[i] = Point(x, y);
        }
    }

    // Calcula el MST usando Kruskal - O(E log E) = O(n^2 log n)
    std::vector<std::pair<char, char>> kruskalMST() {
        std::vector<Edge> edges;

        // Crear lista de aristas - O(n^2)
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (DistanceMatrix[i][j] > 0) {
                    edges.push_back({DistanceMatrix[i][j], i, j});
                }
            }
        }

        // Ordenar aristas - O(n^2 log n^2) = O(n^2 log n)
        std::sort(edges.begin(), edges.end());

        // Inicializar Union-Find - O(n)
        for (int i = 0; i < n; i++) {
            MakeSet(i);
        }

        std::vector<std::pair<char, char>> mstEdges;

        // Construir MST - O(n^2 α(n))
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

    // Resuelve el TSP usando el vecino más cercano - O(n^2)
    std::vector<char> SolveTravellingSalesman() {
        std::vector<int> path;
        TravellingSalesmanNearestNeighbor(path);

        std::vector<char> result;
        for (int node : path) {
            result.push_back('A' + node);
        }

        return result;
    }

    // Calcula el flujo máximo con Ford-Fulkerson - O(E * max_flow) = O(n^2 * f)
    int FordFulkerson() {
        // Crear grafo residual
        std::vector<std::vector<int>> graph = capacityMatrix;
        std::vector<int> parent(n);
        int maxFlow = 0;
        int source = 0, sink = n - 1;

        // Cada BFS cuesta O(n^2) y se ejecuta O(f) veces
        while (BreadthFirstSearch(graph, source, sink, parent)) {
            int pathFlow = INT_MAX;

            // Encontrar capacidad mínima del camino - O(n)
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                pathFlow = std::min(pathFlow, graph[u][v]);
            }

            // Actualizar capacidades residuales - O(n)
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                graph[u][v] -= pathFlow;
                graph[v][u] += pathFlow;
            }

            maxFlow += pathFlow;
        }

        return maxFlow;
    }

    // Encuentra la central más cercana a un punto - O(n)
    Point FindClosestCentral(const Point& query) {
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
    // Parsear (x,y)
    queryLine = queryLine.substr(1, queryLine.length() - 2);
    size_t comma = queryLine.find(',');
    Point query;
    query.x = std::stoi(queryLine.substr(0, comma));
    query.y = std::stoi(queryLine.substr(comma + 1));

    Point closest = optimizer.findClosestCentral(query);
    std::cout << "4.\n(" << closest.x << ", " << closest.y << ")" << std::endl;
    
    return 0;
}
