#include <algorithm>

struct KM {     // KM max-weighted matching
    const static int N = 1010;
    int matchb[N], vb[N], ka[N], kb[N], p[N], c[N], w[N][N], n;
    int inf = 2e9;
    void bfs(int u)
    {
        int a, v = 0, vl = 0, d;
        for (int i = 1; i <= n; i++)    p[i] = 0, c[i] = inf;
        matchb[v] = u;
        do {
            a = matchb[v], d = inf, vb[v] = 1;
            for (int b = 1; b <= n; b++)
                if (!vb[b]) {
                    if (c[b] > ka[a] + kb[b] - w[a][b])
                        c[b] = ka[a] + kb[b] - w[a][b], p[b] = v;
                    if (c[b] < d)    d = c[b], vl = b;
                }
            for (int b = 0; b <= n; b++)
                if (vb[b])      ka[matchb[b]] -= d, kb[b] += d;
                else            c[b] -= d;
            v = vl;
        } while (matchb[v]);
        while (v)               matchb[v] = matchb[p[v]], v = p[v];
    }
    int solve()
    {
        for (int i = 1; i <= n; i++)    matchb[i] = ka[i] = kb[i] = 0;
        for (int a = 1; a <= n; a++)
            std::fill (vb, vb + 1 + n, 0),   bfs(a);
        int res = 0;
        for (int b = 1; b <= n; b++)    res += (w[matchb[b]][b] == -inf ? 0 : w[matchb[b]][b]);
        return res;
    }
};