#pragma once

#include <algorithm>
#include <set>
#include <sstream>
#include <tuple>

std::string itos(int i) { std::stringstream s; s << i; return s.str(); }

// Returns the set of edges from set1 to set2
std::set<std::tuple<int, int>> edges_from_to(const std::set<std::tuple<int, int>>& edges,
    const std::set<int>& set1,
    const std::set<int>& set2) {
    std::set<std::tuple<int, int>> edges_set;
    for (const auto& edge : edges) {
        int i = std::get<0>(edge);
        int j = std::get<1>(edge);
        if ((set1.find(i) != set1.end() && set2.find(j) != set2.end()) ||
            (set2.find(i) != set2.end() && set1.find(j) != set1.end())) {
            edges_set.insert(edge);
        }
    }
    return edges_set;
}


// Returns the set of edges in a set
std::set<std::tuple<int, int>> edges_in(const std::set<std::tuple<int, int>>& edges,
    const std::set<int>& set) {
    std::set<std::tuple<int, int>> edges_set;
    for (const auto& edge : edges) {
        int i = std::get<0>(edge);
        int j = std::get<1>(edge);
        if (set.find(i) != set.end() && set.find(j) != set.end()) {
            edges_set.insert(edge);
        }
    }
    return edges_set;
}


// Returns the set of edges going out of a set
std::set<std::tuple<int, int>> delta_edges(const std::set<std::tuple<int, int>>& edges,
    const std::set<int>& vertices,
    const std::set<int>& set) {
    std::set<int> set_bar;
    std::set_difference(vertices.begin(), vertices.end(), set.begin(), set.end(), std::inserter(set_bar, set_bar.begin()));
    return edges_from_to(edges, set, set_bar);
}
