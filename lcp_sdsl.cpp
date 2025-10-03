// compile with:
// g++ -std=c++17 -O2 -I /u20/hasibih/include lcp_sdsl.cpp -o lcp_sdsl -L /u20/hasibih/lib -lsdsl -ldivsufsort -ldivsufsort64

//Uses compressed suffix arrays / LCP arrays. 
//Preprocessing slightly super-linear due to internal data structures. 
// Extremely fast for medium-sized data, more memory-efficient for large texts.


// O((n + m) log(n + m)) preprocessing time
// O(1) query time


#include <bits/stdc++.h>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp.hpp>

using namespace std;
using namespace sdsl;

class LCPQuery {
private:
    string concat;                 // concatenated string: A + $ + B + #
    size_t nA;
    csa_wt<> csa;                  // compressed suffix array
    lcp_wt<> lcp;                  // LCP support structure

public:
    // constructor: builds SA and LCP
    LCPQuery(const string& A, const string& B) {
        nA = A.size();
        concat = A + "$" + B + "#";  // unique separators not in A/B

        // Write to temporary file
        string temp_file = "/tmp/lcp_temp_" + to_string(time(nullptr));
        ofstream fout(temp_file);
        fout << concat;
        fout.close();

        // build suffix array
        construct(csa, temp_file, 1);

        // build LCP support
        construct(lcp, temp_file, 1);

        // cleanup
        remove(temp_file.c_str());
    }

    // query LCP between A[i:] and B[j:]
    size_t query(size_t i, size_t j) const {
        size_t posA = i;
        size_t posB = nA + 1 + j;  // offset due to separator $

        // Find ranks by searching through suffix array
        size_t rankA = 0, rankB = 0;
        for (size_t k = 0; k < csa.size(); k++) {
            if (csa[k] == posA) rankA = k;
            if (csa[k] == posB) rankB = k;
        }

        if (rankA == rankB) return concat.size() - posA;

        if (rankA > rankB) swap(rankA, rankB);

        // Range minimum query on LCP array
        size_t min_lcp = lcp[rankA + 1];
        for (size_t k = rankA + 2; k <= rankB; k++) {
            min_lcp = min(min_lcp, (size_t)lcp[k]);
        }
        return min_lcp;
    }
};

// ----------------- complex example usage -----------------
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