#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const ll INF = 1e18;
struct Edge { int to, cap, flow; ll cost; };

int N;
vector<Edge> edges;
vector<int> g[105];

void addEdge(int u, int v, int cap, ll cost) {
    g[u].push_back(edges.size());
    edges.push_back({v, cap, 0, cost});
    g[v].push_back(edges.size());
    edges.push_back({u, 0, 0, -cost}); // back-edge
}

vector<ll> dijkstra(int s, vector<ll>& pot) {
    vector<ll> dist(N, INF);
    priority_queue<pair<ll,int>, vector<pair<ll,int>>, greater<>> pq;
    dist[s] = 0; pq.push({0, s});
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        for (int id : g[u]) {
            auto& e = edges[id];
            if (e.cap > e.flow) {
                ll w = e.cost + pot[u] - pot[e.to];
                if (dist[u] + w < dist[e.to]) {
                    dist[e.to] = dist[u] + w;
                    pq.push({dist[e.to], e.to});
                }
            }
        }
    }
    return dist;
}

int edmondsKarp(int s, int t, vector<ll>& dist, vector<ll>& pot_old,
                vector<int>& prevv, vector<int>& preve) {
    fill(prevv.begin(), prevv.end(), -1);
    fill(preve.begin(), preve.end(), -1);
    vector<bool> vis(N, false);
    queue<int> q;
    q.push(s); vis[s] = true;

    while (!q.empty() && !vis[t]) {
        int u = q.front(); q.pop();
        for (int id : g[u]) {
            auto& e = edges[id];
            if (!vis[e.to] && e.cap > e.flow) {
                ll cost_red = e.cost + pot_old[u] - pot_old[e.to];
                if (dist[u] + cost_red == dist[e.to]) { // aresta no caminho ótimo
                    vis[e.to] = true;
                    prevv[e.to] = u; preve[e.to] = id;
                    q.push(e.to);
                }
            }
        }
    }
    if (!vis[t]) return 0;

    int flow = INT_MAX;
    for (int v = t; v != s; v = prevv[v])
        flow = min(flow, edges[preve[v]].cap - edges[preve[v]].flow);
    for (int v = t; v != s; v = prevv[v]) {
        edges[preve[v]].flow   += flow;
        edges[preve[v]^1].flow -= flow;
    }
    return flow;
}

pair<int,ll> solve(int s, int t, int need) {
    int totalFlow = 0; ll totalCost = 0;

    vector<ll> pot(N, INF);
    pot[s] = 0;
    for (int iter = 0; iter < N - 1; iter++)
        for (int u = 0; u < N; u++)
            if (pot[u] < INF)
                for (int id : g[u]) {
                    auto& e = edges[id];
                    if (e.cap > e.flow && pot[u] + e.cost < pot[e.to])
                        pot[e.to] = pot[u] + e.cost;
                }

    vector<int> prevv(N), preve(N);
    while (totalFlow < need) {
        vector<ll> dist = dijkstra(s, pot);
        if (dist[t] == INF) break;

        ll pathCost = dist[t] + pot[t];
        vector<ll> pot_old = pot;
        for (int i = 0; i < N; i++)
            if (dist[i] < INF) pot[i] += dist[i];

        while (totalFlow < need) {
            int pushed = edmondsKarp(s, t, dist, pot_old, prevv, preve);
            if (pushed == 0) break;
            int send = min(pushed, need - totalFlow);
            if (pushed > send) {
                int excess = pushed - send;
                for (int v = t; v != s; v = prevv[v]) {
                    edges[preve[v]].flow   -= excess;
                    edges[preve[v]^1].flow += excess;
                }
            }
            totalFlow += send;
            totalCost += (ll)send * pathCost;
            if (pushed > send) break;
        }
    }
    return {totalFlow, totalCost};
}

int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int M, inst = 1;
    while (cin >> N >> M) {
        for (int i = 0; i < N; i++) g[i].clear();
        edges.clear();
        vector<tuple<int,int,int>> routes(M);
        for (int i = 0; i < M; i++) {
            int a, b, c; cin >> a >> b >> c;
            routes[i] = {a-1, b-1, c};
        }
        int D, K; cin >> D >> K;
        for (auto& [a, b, c] : routes) { addEdge(a, b, K, (ll)c); addEdge(b, a, K, (ll)c); }
        auto [flow, cost] = solve(0, N-1, D);
        cout << "Instancia " << inst++ << "\n";
        if (flow < D) cout << "impossivel\n"; else cout << cost << "\n";
        cout << "\n";
    }
}