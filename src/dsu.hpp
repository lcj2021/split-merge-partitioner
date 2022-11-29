#pragma once
#include <vector>
#include <numeric>


class DSU {
    private:
        std::vector<int> p, sz;
    public: 
        DSU(int n): p(n + 1), sz(n + 1, 1) {std::iota(p.begin(), p.end(), 0); }
        DSU() {}
        int find(int x) {
            return p[x] == x ? p[x] : p[x] = find(p[x]);
        }
        bool same(int x, int y) {return find(x) == find(y); }
        bool merge(int x, int y) {
            x = find(x), y = find(y);
            if (x == y) return false;
            p[y] = x;
            sz[x] += sz[y];
            return true;
        }
        int size(int x) {return sz[find(x)]; }
};