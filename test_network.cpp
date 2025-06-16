#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include "main.cpp"
using namespace std;

void test_Point_Distance() {
    Point p1(0, 0), p2(3, 4);
    assert(abs(p1.Distance(p2) - 5.0) < 1e-6);
    cout << "test_Point_Distance passed.\n";
}

void test_UnionFind() {
    UnionFind uf(5);
    assert(uf.FindSet(1) == 1);
    uf.UnionSets(1, 2);
    assert(uf.FindSet(1) == uf.FindSet(2));
    uf.UnionSets(2, 3);
    assert(uf.FindSet(1) == uf.FindSet(3));
    cout << "test_UnionFind passed.\n";
}

void test_KruskalMST() {
    vector<vector<int>> dist = {{0, 2, 0}, {2, 0, 3}, {0, 3, 0}};
    KruskalMST mst(3, dist);
    auto edges = mst.FindMST();
    assert(edges.size() == 2);
    cout << "test_KruskalMST passed.\n";
}

void test_TSP() {
    vector<vector<int>> dist = {{0, 10, 15, 20}, {10, 0, 35, 25}, {15, 35, 0, 30}, {20, 25, 30, 0}};
    TSPSolver tsp(4, dist);
    auto path = tsp.Solve();
    assert(path.front() == 'A' && path.back() == 'A');
    cout << "test_TSP passed.\n";
}

void test_MaxFlow() {
    vector<vector<int>> dist = {{0, 10, 0, 0}, {0, 0, 5, 15}, {0, 0, 0, 10}, {0, 0, 0, 0}};
    vector<vector<int>> cap = dist;
    MaxFlowSolver mf(4, dist, cap);
    assert(mf.CalculateMaxFlow() == 10);
    cout << "test_MaxFlow passed.\n";
}

void test_FindClosestCentral() {
    vector<Point> points = {Point(0, 0), Point(3, 4), Point(10, 10)};
    CentralManager cm(points);
    Point q(2, 2);
    Point closest = cm.FindClosestCentral(q);
    assert(closest.x == 0 && closest.y == 0);
    cout << "test_FindClosestCentral passed.\n";
}

int main() {
    test_Point_Distance();
    test_UnionFind();
    test_KruskalMST();
    test_TSP();
    test_MaxFlow();
    test_FindClosestCentral();
    cout << "\nAll unit tests passed.\n";
    return 0;
}
