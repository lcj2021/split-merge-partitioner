#pragma once
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unordered_set>
#include <vector>
using namespace std;
class kcore_t
{
    public:
        kcore_t( const int n );
        ~kcore_t();
        void AddEdge( const int u, const int v );
        unordered_set<int> solve( const int k );
    private:
        std::vector<std::vector<int>> adj_list;
        std::vector<int> degree;
        int num_vertices;
};


kcore_t::kcore_t( const int n ): num_vertices( n )
{
    adj_list.resize( n );
    degree.resize( n, 0 );
}

kcore_t::~kcore_t()
{

}

void kcore_t::AddEdge( const int u, const int v )
{
    if ( u == v ) return;
    adj_list[u].push_back( v );
    degree[u]++;
}

/* Batagelj and Zaversnik Algorithm */
unordered_set<int> kcore_t::solve( const int k )
{
    vector<int> vertex( num_vertices, 0 );  // vertices sorted by their degrees in increasing order
    vector<int> deg( num_vertices, 0 );
    vector<int> bin( num_vertices, 0 );     // contains, for each possible degree, the position of the first vertex of that degree in array "vertex"
    vector<int> pos( num_vertices, 0 );     // stores positions of vertices in array "vertex"

    /* Initialization */
    for ( int i = 0; i < num_vertices; ++i ) {
        deg[i] = degree[i];
        bin[deg[i]]++;
    }
    int start = 0;
    for ( int i = 0; i < num_vertices; ++i ) {
        int num = bin[i];
        bin[i] = start;
        start += num;
    }
    // Bin sort vertices by their degree
    for ( int i = 0; i < num_vertices; ++i ) {
        pos[i] = bin[deg[i]];
        vertex[pos[i]] = i;
        bin[deg[i]]++;
    }
    // Recover bin[]
    for ( int i = num_vertices - 1; i >= 1; --i ) bin[i] = bin[i - 1];
    bin[0] = 0;
    
    /* Main Algorithm */
    for ( int i = 0; i < num_vertices; ++i ) {
        int v = vertex[i];
        for ( auto u : adj_list[v] ) {
            if ( deg[u] > deg[v] ) {
                int du = deg[u], pu = pos[u];
                int pw = bin[du];
                int w = vertex[pw];
                if ( u != w ) {
                    pos[u] = pw;
                    pos[w] = pu;
                    vertex[pu] = w;
                    vertex[pw] = u;
                }
                bin[du]++;
                deg[u]--;
            }
        }
    }
    
    /* Output the Result */
    ofstream outputFile( "kcore.txt" );
    if ( !outputFile ) {
        cerr << "Output file could not be opened!" << endl;
        exit( EXIT_FAILURE );
    }
    unordered_set<int> S;
    for ( int i = 0; i < num_vertices; ++i ) {
        if ( deg[i] < k ) continue;
        // outputFile << i << ' ' << deg[i] << '\n';
        S.insert(i);
    }
    outputFile.close();
    return S;
}