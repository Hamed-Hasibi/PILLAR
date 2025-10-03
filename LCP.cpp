#include <bits/stdc++.h>
#include <atcoder/string>   // ACL suffix_array, lcp_array
using namespace std;

struct LCPQuery {
    string A, B, S;
    int nA;
    vector<int> sa, lcp, rank;
    vector<vector<int>> st; // Sparse table
    vector<int> log;

    LCPQuery(const string& A_, const string& B_) {
        A = A_;
        B = B_;
        nA = A.size();
        char sep = '#'; // must not appear in A,B
        S = A + sep + B;

        // SA-IS suffix array (ACL)
        sa = atcoder::suffix_array(S);
        lcp = atcoder::lcp_array(S, sa);

        // rank array
        int n = S.size();
        rank.resize(n);
        for (int i = 0; i < n; i++) rank[sa[i]] = i;

        // Build sparse table on LCP
        buildSparseTable(lcp);
    }

    void buildSparseTable(const vector<int>& arr) {
        int n = arr.size();
        log.assign(n+1,0);
        for (int i=2;i<=n;i++) log[i] = log[i/2]+1;

        int K = log[n]+1;
        st.assign(K, vector<int>(n));
        st[0] = arr;
        for (int k=1;k<K;k++) {
            for (int i=0;i+(1<<k)<=n;i++) {
                st[k][i] = min(st[k-1][i], st[k-1][i+(1<<(k-1))]);
            }
        }
    }

    int rmqQuery(int l, int r) {
        int j = log[r-l+1];
        return min(st[j][l], st[j][r-(1<<j)+1]);
    }

    // LCP of A[i:] and B[j:]
    int query(int i, int j) {
        int posA = i;
        int posB = nA + 1 + j;
        int r1 = rank[posA], r2 = rank[posB];
        if (r1 > r2) swap(r1, r2);
        if (r1 == r2) return (int)S.size() - sa[r1];
        return rmqQuery(r1+1, r2);
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string A = "aba";
    string B = "bab";
    LCPQuery lcpq(A, B);

    cout << lcpq.query(0,1) << "\n"; // "aba" vs "ab" -> 2
    cout << lcpq.query(1,0) << "\n"; // "ba"  vs "bab" -> 2
    cout << lcpq.query(0,0) << "\n"; // "aba" vs "bab" -> 0
}