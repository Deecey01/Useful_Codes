#include <bits/stdc++.h>
using namespace std;
 
#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt")
 
typedef long long int ll;
#define f(i, a, b) for (ll i = a; i < b; i++)
#define rev(v) reverse(v.begin(), v.end())
#define srt(v) sort(v.begin(), v.end())
#define all(v) v.begin(), v.end()
#define mnv(v) *min_element(v.begin(), v.end())
#define mxv(v) *max_element(v.begin(), v.end())
#define vll vector<ll>
#define pb push_back
#define F first 
#define S second 
#define pll pair<ll,ll>
// ----Policy based Data structure----//
//---order_of_key (k) : Number of items strictly smaller than k.---//
//---find_by_order(k) : K-th element in a set (counting from zero).---//
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
template <class T> using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template <class T> using ordered_multiset = tree<T , null_type ,  less_equal<T> , rb_tree_tag , tree_order_statistics_node_update>;
void myerase(ordered_multiset<ll> &t, ll v){
    ll rank = t.order_of_key(v); // Number of elements that are less than v in t
    auto it = t.find_by_order(rank); // Iterator that points to the (rank+1)th element in t
    t.erase(it);
}
ll binpow(ll a, ll b, ll m)
{
    a %= m;
    ll res = 1;
    while (b > 0)
    {
        if (b & 1)
            res = res * a % m;
        a = a * a % m;
        b >>= 1;
    }
    return res;
}

bool sortbysec(const pair<int, int>& a,const pair<int, int>& b){
    return (a.second < b.second);
}

ll maxSubArraySum(vector<ll>a, ll size)//kadane's algorithm
{
    ll max_so_far = INT_MIN, max_ending_here = 0;
 
    for (ll i = 0; i < size; i++) {
        max_ending_here = max_ending_here + a[i];
        if (max_so_far < max_ending_here)
            max_so_far = max_ending_here;
 
        if (max_ending_here < 0)
            max_ending_here = 0;
    }
    return max_so_far;
}

////////////////Centrid of a Tree///////////////////////
vector<ll> Centroid(vll g[],ll n) {
    vector<ll> centroid;
    vector<ll> sz(n);
    function<void (ll, ll)> dfs = [&](ll u,ll prev) {
      sz[u] = 1;
      bool is_centroid = true;
      for (auto v : g[u]) if (v != prev) {
        dfs(v, u);
        sz[u] += sz[v];
        if (sz[v] > n / 2) is_centroid = false;
      }
      if (n - sz[u] > n / 2) is_centroid = false;
      if (is_centroid) centroid.push_back(u);
    };
    dfs(0, -1);
    return centroid;
}

////////////////////////Segment Tree////////////////////
class SGTree {
	vector<ll> seg;
public:
	SGTree(ll n) {
		seg.resize(4 * n + 1);
	}

	void build(ll ind, ll low, ll high, ll arr[]) {
		if (low == high) {
			seg[ind]+= arr[low];
			return;
		}

		ll mid = (low + high) / 2;
		build(2 * ind + 1, low, mid, arr);
		build(2 * ind + 2, mid + 1, high, arr);
		seg[ind] = (seg[2 * ind + 1]+seg[2 * ind + 2]);
	}

	ll query(ll ind, ll low,ll high,ll l, ll r) {
		// no overlap
		// l r low high or low high l r
		if (r < low || high < l) return 0;

		// complete overlap
		// [l low high r]
		if (low >= l && high <= r) return seg[ind];

		ll mid = (low + high) >> 1;
		ll left = query(2 * ind + 1, low, mid, l, r);
		ll right = query(2 * ind + 2, mid + 1, high, l, r);
		return (left+ right);
	}
	void update(ll ind, ll low, ll high, ll i, ll val) {
		if (low == high) {
			seg[ind] = val;
			return;
		}

		ll mid = (low + high) >> 1;
		if (i <= mid) update(2 * ind + 1, low, mid, i, val);
		else update(2 * ind + 2, mid + 1, high, i, val);
		seg[ind] = (seg[2 * ind + 1]+ seg[2 * ind + 2]);
	}
};

void dfs(ll v, ll par, ll h, vector<ll> &d, vll adj[]) {
  d[v] = h;
  for (ll i : adj[v]) {
    if (i != par) {
      dfs(i, v, h + 1, d,adj);
    }
  }
  // queue<ll>q;
  // vll dist(n,-1);
  // q.push(0);
  // dist[0]=0;
  // while(!q.empty()){
  //   ll node=q.front();q.pop();
  //   for(auto x:adj[node]){
  //     if(dist[x]==-1){
  //       dist[x]=1+dist[node];
  //       q.push(x);
  //     }
  //   }
  // }
}

//////////////////////////////////Bipartite//////////////////////////////
bool bipartite(int n,vector<int>adj[], int col[]){
  queue<pair<int,int>>q;
  f(i,0,n){
    if(col[i]==-1){
      q.push({i,1});
      col[i]=1;
      pair<int,int>p;
      while (!q.empty()){
        p=q.front();
        q.pop();
        int v=p.first;
        int c=p.second;
        for(auto x: adj[v]){
          if (col[x]==c) return false;
          if(col[x]==-1){
            if (c==2)col[x]=1;
            else col[x]=2;
            q.push({x,col[x]});
          }
        }
      }
    }
  }
  return true;
}

////////////////////Prim's/////////////////////
priority_queue<pll, vector<pll> , greater<pll> > pq;
ll src = 0; 
vector<ll> key(n,1e18);
vector<ll> parent(n, -1);
vector<bool> inMST(n, false);
pq.push({0, src});
key[src] = 0;
while (!pq.empty()){
  ll u = pq.top().second;
  pq.pop();
  if(inMST[u] == true){
    continue;
  }
  inMST[u] = true;
  for(auto x:adj[u]){
    ll v = x.F;
    ll weight = x.S;
    if (inMST[v] == false && key[v] > weight){
      key[v] = weight;
      pq.push({key[v],v});
      parent[v] = u;
    }
  }
}
f(i,1,n){
  cout<<ind[{parent[i],i}]<<" ";
}

///////////////////////////////Bellman Ford//////////////////////////////
ll bellman_ford(ll V, vector<vll>& edges, ll S) {
	vll dist(V, 1e18);
	dist[S] = 0;
	for (ll i = 0; i < V - 1; i++) {
		for (auto it : edges) {
			ll u = it[0];
			ll v = it[1];
			ll wt = it[2];
			if (dist[u] != 1e18 && dist[u] + wt < dist[v]) {
				dist[v] = dist[u] + wt;
			}
		}
	}
	// Nth relaxation to check negative cycle
  f(i,0,V){
  	for (auto it : edges) {
  		ll u = it[0];
  		ll v = it[1];
  		ll wt = it[2];
  		if (dist[u] != 1e18 && dist[u] + wt < dist[v]){
  			  return -1;
  		}
  	}
  }
}

/////////////////////////////////////Disjoint set Union///////////////////////////
struct DSU {
	vector<int> e;
	void init(int n) { e = vector<int>(n, -1); }
	int get(int x) { return (e[x] < 0 ? x : e[x] = get(e[x])); }
	bool sameSet(int x, int y) { return get(x) == get(y); }
	int size(int x) { return -e[get(x)]; }
	bool unite(int x, int y) {
		x = get(x), y = get(y);
		if (x == y) return 0;
		if (e[x] > e[y]) swap(x, y);
		e[x] += e[y];
		e[y] = x;
		return 1;
	}
};

////////////////////Topological Sort/////////////////////
vll topoSort(ll V, vll adj[]){
	ll indegree[V] = {0};
	for (ll i = 0; i < V; i++) {
		for (auto it : adj[i]) {
			indegree[it]++;
		}
	}
 
	queue<ll> q;
	for (ll i = 0; i < V; i++) {
		if (indegree[i] == 0) {
			q.push(i);
		}
	}
	vll topo;
	while (!q.empty()) {
		ll node = q.front();
		q.pop();
		topo.push_back(node);
 
		for (auto it : adj[node]) {
			indegree[it]--;
			if (indegree[it] == 0) q.push(it);
		}
	}
	return topo;
}
