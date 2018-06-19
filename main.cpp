#include <string.h>
#include <ctime>
#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#define TESTMODE

using namespace std;

const int maxInteger = 2 << 20;
const int NMax = 20000;

int maximum(int a, int b)
{
    return a > b ? a : b;
}

class edge
{
protected:
    int vertex;
    int edgeNumber;
    int color;

public:
    edge(int v)
    {
        vertex = v;
        color = -1;
    }

    edge(int v, int c)
    {
        vertex = v;
        color = c;
    }

    int getEdgeNumber()
    {
        return edgeNumber;
    }

    void setEdgeNumber(int x)
    {
        edgeNumber = x;
    }

    void setColor(int c)
    {
        color = c;
    }

    int getVertex()
    {
        return vertex;
    }

    int getColor()
    {
        return color;
    }
};

class Graph
{
    int n, m;
    vector <edge> graf[NMax];


public:
    #ifdef TESTMODE
    void generateTree(int siz, int maxRand)
    {
        setn(siz);
        setm(siz - 1);
        srand (time(NULL));
        int node = 0;
        int connectedTo = 1;
        while(connectedTo < siz)
        {
            int upTo = connectedTo + rand() % maxRand;
            while(connectedTo < siz && connectedTo <= upTo)
            {
                addEdge(node, connectedTo, connectedTo - 1);
                ++connectedTo;
            }
            ++node;
        }

    }

    void generateGraph(int siz, int muchii)
    {
        if(muchii < siz - 1)
        {
            cerr << "Parametri imposibili";
        }
        setn(siz);
        setm(muchii);
        srand(time(NULL));
        int remain = muchii;
        int edgeNumber = n - 1;
        for(int i = 0; i < n - 1; ++i)
        {
            addEdge(i, i + 1, i);
            --remain;
        }
        while(remain)
        {
            int a = rand() % n;
            int b = rand() % n;
            int ok = 0;
            for(int i = 0; i < graf[a].size(); ++i)
                if(graf[a][i].getVertex() == b)
                    ok = 1;
            if(!ok)
            {
                --remain;
                addEdge(a, b, edgeNumber);
                ++edgeNumber;
            }
        }
    }
    #endif
    int maxColorGroup()
    {
        vector<int> colors(m + 2, 0);
        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < graf[i].size(); ++j)
            {
                ++colors[graf[i][j].getColor() + 1];
            }
        }
        int result = 0;
        for(int i = 0; i <= m; ++i)
        {
            if(result < colors[i])
                result = colors[i];
        }
        return result / 2;
    }

    void addEdge(int a, int b, int x)
    {
        edge ab(b), ba(a);
        ab.setEdgeNumber(x);
        ba.setEdgeNumber(x);
        graf[a].push_back(ab);
        graf[b].push_back(ba);

    }

    void setn(int x)
    {
        n = x;
    }

    int getn()
    {
        return n;
    }

    void setm(int x)
    {
        m = x;
    }

    int getm()
    {
        return m;
    }

    void afisare()
    {
        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < graf[i].size(); ++j)
                cout<<graf[i][j].getVertex()<<' ';
            cout<<'\n';
        }
    }

    void colorEdges(vector<int> &coloring)
    {
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < graf[i].size(); ++j)
                graf[i][j].setColor(coloring[graf[i][j].getEdgeNumber()]);
    }

    void afisareCuCulori()
    {
        string colorName[20] = {"blue", "red", "orange", "green", "yellow", "black", "white"};

        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < graf[i].size(); ++j)
                cout<< i << ' ' << graf[i][j].getColor() << ' ' << graf[i][j].getVertex()<<'\n';
            cout<<'\n';
        }
    }

    int isTree()
    {
        vector<int> parent(n, -1);
        vector<int> visited(n, 0);
        queue<int> q;
        q.push(0);
        while(!q.empty())
        {
            int currentVertex = q.front();
            visited[currentVertex] = 1;
            for(int i = 0; i < graf[currentVertex].size(); ++i)
            {
                int vertex2 = graf[currentVertex][i].getVertex();
                if(parent[currentVertex] != vertex2 && visited[vertex2])
                    return 0;
                if(parent[currentVertex] != vertex2)
                {
                    parent[vertex2] = currentVertex;
                    q.push(vertex2);
                }
            }
            q.pop();
        }
        return 1;
    }

    void knapsack(vector<int> &indices, vector<int> &w, vector<int> &v, int &W, vector<int> &solution)
    {
        int number = indices.size();
        int dp[number + 1][W + 1];
        for(int j = 0; j < w[0]; ++j)
            dp[0][j] = 0;
        for(int j = w[0]; j <= W; ++j)
            dp[0][j] = v[0];
        for(int i = 1; i < number; ++i)
            for(int j = 0; j <= W; ++j)
                if(w[i] > j)
                {
                    dp[i][j] = dp[i - 1][j];
                }
                else
                    dp[i][j] = maximum(dp[i - 1][j], dp[i - 1][j - w[i]] + v[i]);
      //  vector<int> solution;
        int i = number - 1;
        int j = W;
        while(i >= 0)
        {
            if(i == 0)
            {
                if(dp[i][j] > 0)
                    solution.push_back(indices[i]);
                break;
            }
            if(dp[i][j] == dp[i - 1][j])
                --i;
            else
            {
                solution.push_back(indices[i]);
                j -= w[i];
                --i;
            }
        }
   //     return solution;
    }

    void algorithm4()
    {
        //delta e gradul maxim dintre nodurile grafului. Pasul 1: m = delta - 1
        int delta = 0;
        int mm = 0;
        for(int i = 0; i < n; ++i)
            if(delta < graf[i].size())
                delta = graf[i].size();
        mm = delta - 1;

        vector<int> U(m, -1);
        vector<int> M(m, -1);
        vector<int> coloring(m, -1);
        int colored = 0;
        for(int root = 0; root < n; ++root)
        {
            int nextColor = 0;
            vector<int> isLeaf(n, 1);
            vector<int> v(n, 0);
            v[root] = 1;
            vector<vector <int> > layers;
            vector<int> layer(1, root);
            vector<int> distance(n, 0);

            while(layer.size() != 0)
            {
                vector<int> aux;
                for(int i = 0; i < layer.size(); ++i)
                {
                    int currentVertex = layer[i];
                    for(int j = 0; j < graf[currentVertex].size(); ++j)
                    {
                        int vertex2 = graf[currentVertex][j].getVertex();
                        if(!v[vertex2])
                        {
                            v[vertex2] = 1;
                            aux.push_back(vertex2);
                            distance[vertex2] = distance[currentVertex] + 1;
                        }
                        else
                            isLeaf[vertex2] = 0;
                    }
                }
                layers.push_back(layer);
                layer = aux;
            }
            int l = delta / 2;
            int u = mm;

            while(true)
            {
                vector<int> currentColoring(m, -1);
                vector<int> residual(n, 1);
                residual[root] = 0;
                int c = (l + u) / 2;
                for(int d = layers.size() - 2; d >= 0; --d)
                {
                    for(int i = 0; i < layers[d].size(); ++i)
                    {
                        int currentVertex = layers[d][i];
                        if(isLeaf[currentVertex])
                        {
                            continue;
                        }
                        vector<int> indices;
                        vector<int> weights;
                        int S1 = 0, S2 = 0;
                        for(int j = 0; j < graf[currentVertex].size(); ++j)
                        {
                            int vertex2 = graf[currentVertex][j].getVertex();
                            if(distance[vertex2] > distance[currentVertex])
                            {
                                indices.push_back(vertex2);
                                weights.push_back(residual[vertex2]);
                                S1 += residual[vertex2];
                            }
                        }

                        vector<int> S;
                        knapsack(indices, weights, weights, c, S);

                        for(int j = 0; j < S.size(); ++j)
                            S2 += residual[S[j]];
                        if(S1 - S2 + residual[currentVertex] > c)
                        {
                            goto step16;
                        }
                        //BFS in jos de la currentVertex
                        queue<int> BFS;
                        vector<int> viz(n, 0);
                        for(int j = 0; j < S.size(); ++j)
                        {
                            viz[S[j]] = 1;
                            BFS.push(S[j]);
                        }
                        viz[currentVertex] = 1;
                        while(!BFS.empty())
                        {
                            int x = BFS.front();
                            for(int j = 0; j < graf[x].size(); ++j)
                            {
                                if(currentColoring[graf[x][j].getEdgeNumber()] == -1)
                                    currentColoring[graf[x][j].getEdgeNumber()] = nextColor;
                                if(!viz[graf[x][j].getVertex()])
                                {
                                    BFS.push(graf[x][j].getVertex());
                                    viz[graf[x][j].getVertex()] = 1;
                                }
                            }
                            BFS.pop();
                        }
                        ++nextColor;
                        residual[currentVertex] = residual[currentVertex] + S1 - S2;
                    }
                }

                for(int i = 0; i < graf[root].size(); ++i)
                {
                    if(currentColoring[graf[root][i].getEdgeNumber()] == -1)
                        currentColoring[graf[root][i].getEdgeNumber()] = nextColor;
                }

                if(l == u)
                {
                    U = currentColoring;
                    break;
                }
                u = c;
                continue;

                step16:
                if(l == u)
                {
                    M = U;
                    break;
                }
                l = c + 1;
            }

            if(!colored)
            {
                colored = 1;
                M = U;
            }

            if(u < mm)
            {
                M = U;
                mm = u;
            }
        }
        colorEdges(M);
    }

    void algorithm5()
    {
        int mmin = maxInteger;
        int vmin;
        for(int i = 0; i < n; ++i)
        {
            vector<int> layer(1, i);
            vector<int> v(m, 0);
            int cv = -1; //numarul maxim de muchii la aceasi distanta fata de i
            while(layer.size() != 0)
            {
                vector<int> aux;
                int cd = 0; //contor dupa distanta
                for(int j = 0; j < layer.size(); ++j)
                {
                    for(int k = 0; k < graf[layer[j]].size(); ++k)
                    {
                        int e = graf[layer[j]][k].getEdgeNumber();
                        if(!v[e])
                        {
                            ++cd;
                            aux.push_back(graf[layer[j]][k].getVertex());
                            v[e] = 1;
                        }

                    }
                }
                if(cv < cd)
                    cv = cd;
                layer = aux;
            }
            if(cv < mmin)
            {
                mmin = cv;
                vmin = i;
            }
        }
        //coloram muchiile dupa distanta de la vmin
        vector<int> layer(1, vmin);
        vector<int> v(m, 0);
        int nextColor = 0;
        vector<int> colors(m, -1);
        while(layer.size() != 0)
        {
            vector<int> aux;
            for(int j = 0; j < layer.size(); ++j)
            {
                for(int k = 0; k < graf[layer[j]].size(); ++k)
                {
                    int e = graf[layer[j]][k].getEdgeNumber();
                    if(!v[e])
                    {
                        colors[e] = nextColor;
                        aux.push_back(graf[layer[j]][k].getVertex());
                        v[e] = 1;
                    }

                }
            }
            layer = aux;
            ++nextColor;
        }
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < graf[i].size(); ++j)
                graf[i][j].setColor(colors[graf[i][j].getEdgeNumber()]);
    }



};

std::istream& operator >> (std::istream& f, Graph& G)
{
    int n, m;
    f >> n >> m;
    G.setn(n);
    G.setm(m);
    for(int i = 0; i < m; ++i)
    {
        int a, b;
        f >> a >> b;
        G.addEdge(a, b, i);
    }
}

int main()
{
    freopen("testlog.out", "w", stdout);
//    cout << "ARBORI BINARI:\n";
//    for(int n = 100; n <= 5000; n += 100)
//    {
//        Graph G;
//        G.generateTree(n, 2);
//        clock_t before = clock();
//        G.algorithm4();
//        clock_t after = clock();
//        double duration = (after - before) / (double)CLOCKS_PER_SEC;
//        int maxColorGroup = G.maxColorGroup();
//        cout << n << " noduri: " << duration << ' ' << maxColorGroup << '\n';
//    }
//    cout << "ARBORI NE-BINARI:\n";
//    for(int n = 4800; n <= 5000; n += 200)
//    {
//        for(int maxRand = 5; maxRand <= 105; maxRand += 10)
//        {
//            Graph G;
//            G.generateTree(n, maxRand);
//            clock_t before = clock();
//            G.algorithm4();
//            clock_t after = clock();
//            double duration = (after - before) / (double)CLOCKS_PER_SEC;
//            int maxColorGroup = G.maxColorGroup();
//            cout << n << " " << maxRand << " " << duration << " " << maxColorGroup << '\n';
//        }
//    }
    {
        for(int nrNoduri = 823; nrNoduri <= 1000; nrNoduri *= 1.5f )
        {
            for(int nrMuchii = 106575 * 1.5f; nrMuchii < (nrNoduri * (nrNoduri - 1)) / 2; nrMuchii *= 1.5f)
            {
                Graph G;
                G.generateGraph(nrNoduri, nrMuchii);
                clock_t before = clock();
                G.algorithm5();
                clock_t after = clock();
                double duration = (after - before) / (double)CLOCKS_PER_SEC;
                int maxColorGroup = G.maxColorGroup();
                cout << nrNoduri << " & " << nrMuchii << " & " << duration << " & " << maxColorGroup << " \\\\ \n";
            }
        }
    }

    return 0;
}
