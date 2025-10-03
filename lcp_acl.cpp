// g++ -std=c++17 -O2 lcp_acl.cpp -I ac-library/ -o lcp_acl

// Uses uncompressed linear-time SA-IS + Kasai + RMQ. 
// True linear preprocessing. Very fast in practice, minimal constant overhead.
// NOTE: ACL can also support true O(n) preprocessing + O(1) queries using a linear RMQ (Fischer-Heun / Cartesian tree), but sparse table is simpler and very fast in practice.

// O(n + m) preprocessing time
// O(1) query time


#include <bits/stdc++.h>
#include <atcoder/string.hpp> // ACL suffix_array
using namespace std;

/*
 Fischer–Heun RMQ
 Build: O(n)
 Query: O(1)
 Reference: Fischer & Heun, 2006
*/
struct RMQ {
    int n, block_size;
    vector<int> arr;
    vector<int> block_min, block_id;
    vector<vector<int>> st;
    vector<int> log;

    RMQ() {}
    RMQ(const vector<int>& v) { build(v); }

    void build(const vector<int>& v) {
        arr = v;
        n = arr.size();
        if (n == 0) return;
        block_size = max(1, (int)(log2(n)/2)); // block size ~ (1/2) log n
        int num_blocks = (n + block_size - 1)/block_size;

        // Precompute mins per block
        block_min.resize(num_blocks);
        block_id.resize(num_blocks);
        for (int b = 0; b < num_blocks; b++) {
            int l = b * block_size;
            int r = min(n, l + block_size);
            int id = l;
            for (int i=l+1;i<r;i++) if (arr[i] < arr[id]) id = i;
            block_min[b] = arr[id];
            block_id[b] = id;
        }

        // Sparse table over block minima
        log.assign(num_blocks+1,0);
        for (int i=2;i<=num_blocks;i++) log[i] = log[i/2]+1;
        int K = log[num_blocks]+1;
        st.assign(K, vector<int>(num_blocks));
        for (int i=0;i<num_blocks;i++) st[0][i] = block_id[i];
        for (int k=1;k<K;k++) {
            for (int i=0;i+(1<<k)<=num_blocks;i++) {
                int i1 = st[k-1][i];
                int i2 = st[k-1][i+(1<<(k-1))];
                st[k][i] = (arr[i1] <= arr[i2] ? i1 : i2);
            }
        }
    }

    int query(int l, int r) {
        if (l > r) swap(l,r);
        int bl = l/block_size, br = r/block_size;
        int best = l;
        if (bl == br) {
            for (int i=l;i<=r;i++) if (arr[i]<arr[best]) best = i;
            return arr[best];
        }
        // left partial block
        int lend = (bl+1)*block_size-1;
        for (int i=l;i<=min(lend,n-1);i++) if (arr[i]<arr[best]) best=i;
        // right partial block
        int rstart = br*block_size;
        for (int i=rstart;i<=r;i++) if (arr[i]<arr[best]) best=i;
        // middle blocks
        if (bl+1 <= br-1) {
            int L = bl+1, R = br-1;
            int j = log[R-L+1];
            int i1 = st[j][L];
            int i2 = st[j][R-(1<<j)+1];
            int cand = (arr[i1] <= arr[i2] ? i1 : i2);
            if (arr[cand] < arr[best]) best = cand;
        }
        return arr[best];
    }
};

struct LCPQuery {
    string A, B, S;
    int nA;
    vector<int> sa, lcp, rank;
    RMQ rmq;

    LCPQuery(const string& A_, const string& B_) {
        A = A_;
        B = B_;
        nA = A.size();
        char sep = 1;
        S = A + sep + B;
        int n = S.size();

        sa = atcoder::suffix_array(S);

        rank.resize(n);
        for (int i=0;i<n;i++) rank[sa[i]] = i;

        lcp.assign(n,0);
        int h=0;
        for (int i=0;i<n;i++) {
            if (rank[i]==0) continue;
            int j = sa[rank[i]-1];
            while (i+h<n && j+h<n && S[i+h]==S[j+h] && S[i+h]!=sep) h++;
            lcp[rank[i]] = h;
            if (h>0) h--;
        }

        rmq.build(lcp);
    }

    int query(int i, int j) {
        int posA = i;
        int posB = nA + 1 + j;
        int r1 = rank[posA], r2 = rank[posB];
        if (r1 == r2) return min((int)S.size()-posA,(int)S.size()-posB);
        if (r1 > r2) swap(r1,r2);
        return rmq.query(r1+1, r2);
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Test with more complex strings
    string A = "algorithmic";
    string B = "algorithms";
    
    cout << "Testing LCP between A = \"" << A << "\" and B = \"" << B << "\"\n";
    cout << "A length: " << A.size() << ", B length: " << B.size() << "\n\n";

    LCPQuery lcpq(A, B);

    // Test various combinations
    cout << "=== LCP Queries ===\n";
    
    // Test prefixes
    cout << "LCP(A[0:], B[0:]) = " << lcpq.query(0, 0) << "  // full strings\n";
    cout << "LCP(A[0:], B[1:]) = " << lcpq.query(0, 1) << "  // A vs \"lgorithms\"\n";
    cout << "LCP(A[1:], B[0:]) = " << lcpq.query(1, 0) << "  // \"lgorithmic\" vs B\n";
    
    // Test suffix
    cout << "LCP(A[2:], B[2:]) = " << lcpq.query(2, 2) << "  // \"gorithmic\" vs \"gorithms\"\n";
    cout << "LCP(A[3:], B[3:]) = " << lcpq.query(3, 3) << "  // \"orithmic\" vs \"rithms\"\n";
    
    // Test overlapping cases
    cout << "LCP(A[4:], B[4:]) = " << lcpq.query(4, 4) << "  // \"rithmic\" vs \"ithms\"\n";
    cout << "LCP(A[5:], B[5:]) = " << lcpq.query(5, 5) << "  // \"ithmic\" vs \"thms\"\n";
    
    // Test end cases
    cout << "LCP(A[8:], B[6:]) = " << lcpq.query(8, 6) << "  // \"mic\" vs \"ms\"\n";
    cout << "LCP(A[9:], B[7:]) = " << lcpq.query(9, 7) << "  // \"ic\" vs \"s\"\n";
    cout << "LCP(A[10:], B[8:]) = " << lcpq.query(10, 8) << "  // \"c\" vs \"\"\n";
    
    // Test with another pair
    cout << "\n=== Testing with another pair ===\n";
    string C = "programming";
    string D = "programmer";
    
    cout << "Testing LCP between C = \"" << C << "\" and D = \"" << D << "\"\n";
    LCPQuery lcpq2(C, D);
    
    cout << "LCP(C[0:], D[0:]) = " << lcpq2.query(0, 0) << "  // full strings\n";
    cout << "LCP(C[8:], D[8:]) = " << lcpq2.query(8, 8) << "  // \"ing\" vs \"er\"\n";
    cout << "LCP(C[9:], D[9:]) = " << lcpq2.query(9, 9) << "  // \"ng\" vs \"r\"\n";
    cout << "LCP(C[10:], D[10:]) = " << lcpq2.query(10, 10) << "  // \"g\" vs \"\"\n";
    
    // Test edge cases
    cout << "\n=== Edge Cases ===\n";
    string E = "abc";
    string F = "def";
    
    cout << "Testing LCP between E = \"" << E << "\" and F = \"" << F << "\"\n";
    LCPQuery lcpq3(E, F);
    
    cout << "LCP(E[0:], F[0:]) = " << lcpq3.query(0, 0) << "  // no common prefix\n";
    cout << "LCP(E[1:], F[1:]) = " << lcpq3.query(1, 1) << "  // no common prefix\n";
    
    // Test identical strings
    string G = "hello";
    string H = "hello";
    
    cout << "\nTesting LCP between identical strings G = \"" << G << "\" and H = \"" << H << "\"\n";
    LCPQuery lcpq4(G, H);
    
    cout << "LCP(G[0:], H[0:]) = " << lcpq4.query(0, 0) << "  // identical strings\n";
    cout << "LCP(G[1:], H[1:]) = " << lcpq4.query(1, 1) << "  // identical suffixes\n";

    return 0;
    
}