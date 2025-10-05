#ifndef LCP_LIBSAIS_H
#define LCP_LIBSAIS_H

#include <bits/stdc++.h>
extern "C" {
    #include "libsais.h"
}
using namespace std;

/*
 Fischer–Heun RMQ (O(n) build, O(1) query)
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
        block_size = max(1, (int)(log2(n)/2));
        int num_blocks = (n + block_size - 1) / block_size;

        block_min.resize(num_blocks);
        block_id.resize(num_blocks);
        for (int b = 0; b < num_blocks; b++) {
            int l = b * block_size;
            int r = min(n, l + block_size);
            int id = l;
            for (int i = l + 1; i < r; i++) if (arr[i] < arr[id]) id = i;
            block_min[b] = arr[id];
            block_id[b] = id;
        }

        log.assign(num_blocks + 1, 0);
        for (int i = 2; i <= num_blocks; i++) log[i] = log[i / 2] + 1;
        int K = log[num_blocks] + 1;
        st.assign(K, vector<int>(num_blocks));
        for (int i = 0; i < num_blocks; i++) st[0][i] = block_id[i];
        for (int k = 1; k < K; k++) {
            for (int i = 0; i + (1 << k) <= num_blocks; i++) {
                int i1 = st[k - 1][i];
                int i2 = st[k - 1][i + (1 << (k - 1))];
                st[k][i] = (arr[i1] <= arr[i2] ? i1 : i2);
            }
        }
    }

    int query(int l, int r) {
        if (l > r) swap(l, r);
        int bl = l / block_size, br = r / block_size;
        int best = l;
        if (bl == br) {
            for (int i = l; i <= r; i++) if (arr[i] < arr[best]) best = i;
            return arr[best];
        }
        int lend = (bl + 1) * block_size - 1;
        for (int i = l; i <= min(lend, n - 1); i++) if (arr[i] < arr[best]) best = i;
        int rstart = br * block_size;
        for (int i = rstart; i <= r; i++) if (arr[i] < arr[best]) best = i;
        if (bl + 1 <= br - 1) {
            int L = bl + 1, R = br - 1;
            int j = log[R - L + 1];
            int i1 = st[j][L];
            int i2 = st[j][R - (1 << j) + 1];
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
        char sep = 1; // separator smaller than any normal char
        S = A + sep + B;
        int n = S.size();

        // --- Build SA with LibSAIS ---
        sa.assign(n, 0);
        vector<uint8_t> T(n);
        for (int i = 0; i < n; i++) T[i] = (uint8_t)(unsigned char)S[i];
        libsais(T.data(), sa.data(), n, 0, nullptr);

        // --- Rank ---
        rank.resize(n);
        for (int i = 0; i < n; i++) rank[sa[i]] = i;

        // --- Kasai LCP ---
        lcp.assign(n, 0);
        int h = 0;
        for (int i = 0; i < n; i++) {
            if (rank[i] == 0) continue;
            int j = sa[rank[i] - 1];
            while (i + h < n && j + h < n && S[i + h] == S[j + h] && S[i + h] != sep) h++;
            lcp[rank[i]] = h;
            if (h > 0) h--;
        }

        // --- RMQ ---
        rmq.build(lcp);
    }

    int query(int i, int j) {
        int posA = i;
        int posB = nA + 1 + j;
        int r1 = rank[posA], r2 = rank[posB];
        if (r1 == r2) return min((int)S.size() - posA, (int)S.size() - posB);
        if (r1 > r2) swap(r1, r2);
        return rmq.query(r1 + 1, r2);
    }
};

#endif // LCP_LIBSAIS_H
