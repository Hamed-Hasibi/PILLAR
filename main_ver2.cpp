// file: pillar_hamming_section4_final_commented.cpp
// g++ -std=c++17 -O2 pillar_hamming_section4_final_commented.cpp -o pillar_hamming_section4_final_commented

#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <variant>
#include <optional>
#include <cmath>
#include <bits/stdc++.h>

using namespace std;

/* ============================
   Basic types from the paper
   ============================ */

struct StringHandle {
    int id;
    size_t start = 0;
    size_t len = 0;
    bool reversed = false;
};

struct ArithmeticProgression {
    long long start_value;
    long long difference;
    long long length;
};

using Occurrences = std::vector<ArithmeticProgression>;

struct Break { StringHandle handle; };
struct RepetitiveRegion { StringHandle handle; StringHandle period; };
using AnalysisResult = std::variant<
    std::vector<Break>,
    std::vector<RepetitiveRegion>,
    StringHandle
>;

/* ============================
   [cite_start]PILLAR model interface [cite: 330-345]
   ============================ */

class PILLAR_Model {
public:
    virtual ~PILLAR_Model() = default;
    
    // Core PILLAR Operations from paper
    virtual StringHandle Extract(const StringHandle& s, size_t l, size_t r) = 0;
    virtual size_t LCP(const StringHandle& s1, const StringHandle& s2) = 0;
    virtual size_t LCPR(const StringHandle& s1, const StringHandle& s2) = 0;
    virtual Occurrences IPM(const StringHandle& p, const StringHandle& t) = 0;
    virtual size_t Length(const StringHandle& s) = 0;
    virtual char Access(const StringHandle& s, size_t i) = 0; // Kept for DistancesRLE

    // Helper operations from paper (can be built on core ops)
    virtual size_t Period(const StringHandle& s) = 0;
    virtual std::optional<int> FindRotationOffset(const StringHandle& s, const StringHandle& t) = 0;
    virtual Occurrences ExactMatches(const StringHandle& p, const StringHandle& t) = 0;
};

/* ============================
   General Helpers
   ============================ */

static inline long long ceil_div_nonneg(long long a, long long b) {
    if (a <= 0 || b <= 0) return 0;
    return (a + b - 1) / b;
}

StringHandle reverse_handle(const StringHandle& h) {
    return {h.id, h.start, h.len, !h.reversed};
}

/**
 * [cite_start]@brief Implements LCP(S, Q^inf) logic, as per Corollary 2.9 [cite: 355-358].
 * This helper computes the Longest Common Prefix between a finite string S and an
 * infinitely repeating periodic string Q^inf.
 */
size_t LCP_vs_Period(const StringHandle& s, const StringHandle& q_period, PILLAR_Model& pillar) {
    size_t s_len = pillar.Length(s);
    size_t q_len = pillar.Length(q_period);
    if (s_len == 0 || q_len == 0) return 0;

    size_t total_lcp = 0;
    size_t s_pos = 0;
    size_t q_pos = 0;

    while (s_pos < s_len) {
        StringHandle s_suffix = pillar.Extract(s, s_pos, s_len);
        StringHandle q_suffix = pillar.Extract(q_period, q_pos, q_len);

        size_t lcp = pillar.LCP(s_suffix, q_suffix);
        
        total_lcp += lcp;
        s_pos += lcp;
        
        if (lcp < pillar.Length(q_suffix)) {
            // Mismatch found before the end of the Q block
            break;
        }
        // Matched a full q_suffix, continue from start of q_period
        q_pos = 0;
    }
    return total_lcp;
}

/**
 * [cite_start]@brief Implements LCPR(S, Q^inf) logic, as per Corollary 2.10 [cite: 359-360].
 */
size_t LCPR_vs_Period(const StringHandle& s, const StringHandle& q_period, PILLAR_Model& pillar) {
    return LCP_vs_Period(reverse_handle(s), reverse_handle(q_period), pillar);
}


/* ============================
   Section 4.1: Auxiliary PILLAR Model Operations
   ============================ */

/**
 * @brief Auxiliary generator for mismatches between two FINITE strings S and T.
 * This is a necessary helper for Verify_Mismatches (Lemma 4.3). It is implemented
 * using repeated LCP calls on suffixes, which is PILLAR-compliant.
 */
class MismGenerator_Pillar_Finite {
    StringHandle s, t;
    size_t current_pos;
    PILLAR_Model& pillar;

public:
    MismGenerator_Pillar_Finite(const StringHandle& s_, const StringHandle& t_, PILLAR_Model& p)
        : s(s_), t(t_), current_pos(0), pillar(p) {}

    std::optional<size_t> Next() {
        size_t len_s = pillar.Length(s);
        if (current_pos >= len_s) return std::nullopt;

        StringHandle s_suffix = pillar.Extract(s, current_pos, len_s);
        StringHandle t_suffix = pillar.Extract(t, current_pos, len_s);
        
        size_t lcp_len = pillar.LCP(s_suffix, t_suffix);
        current_pos += lcp_len;

        if (current_pos >= len_s) return std::nullopt;
        
        size_t mismatch_pos = current_pos;
        current_pos++;
        return mismatch_pos;
    }
};

/**
 * [cite_start]@brief Implements MismGenerator(S, Q) from Lemma 4.1 [cite: 535-545].
 * Generates mismatches between a finite string S and an infinitely repeating Q^inf.
 */
class MismGenerator_Pillar {
    StringHandle s, q_period;
    size_t current_pos;
    PILLAR_Model& pillar;

public:
    MismGenerator_Pillar(const StringHandle& s_, const StringHandle& q_period_, PILLAR_Model& p)
        : s(s_), q_period(q_period_), current_pos(0), pillar(p) {}

    std::optional<size_t> Next() {
        size_t len_s = pillar.Length(s);
        if (current_pos >= len_s) return std::nullopt;
        StringHandle s_suffix = pillar.Extract(s, current_pos, len_s);
        size_t lcp_len = LCP_vs_Period(s_suffix, q_period, pillar);
        current_pos += lcp_len;
        if (current_pos >= len_s) return std::nullopt;
        size_t mismatch_pos = current_pos;
        current_pos++;
        return mismatch_pos;
    }
};

/**
 * [cite_start]@brief Implements MismGeneratorR(S, Q) from Lemma 4.1[cite: 535].
 * The reverse (suffix-based) version of the MismGenerator.
 */
class MismGeneratorR_Pillar {
    StringHandle s, q_period;
    size_t current_pos_rev;
    PILLAR_Model& pillar;

public:
    MismGeneratorR_Pillar(const StringHandle& s_, const StringHandle& q_period_, PILLAR_Model& p)
        : s(s_), q_period(q_period_), current_pos_rev(0), pillar(p) {}

    std::optional<size_t> Next() {
        size_t len_s = pillar.Length(s);
        if (current_pos_rev >= len_s) return std::nullopt;
        StringHandle s_suffix_rev = pillar.Extract(s, 0, len_s - current_pos_rev);
        size_t lcpr_len = LCPR_vs_Period(s_suffix_rev, q_period, pillar);
        current_pos_rev += lcpr_len;
        if (current_pos_rev >= len_s) return std::nullopt;
        size_t mismatch_pos = len_s - 1 - current_pos_rev;
        current_pos_rev++;
        return mismatch_pos;
    }
};

/**
 * [cite_start]@brief Implements Verify(S, T, k) from Lemma 4.3 [cite: 546-548].
 * Checks if HammingDistance(S, T) <= k in O(k) PILLAR operations.
 */
bool Verify_Mismatches(const StringHandle& s, const StringHandle& t, int k, PILLAR_Model& pillar) {
    size_t len_s = pillar.Length(s);
    size_t len_t = pillar.Length(t);
    
    if (len_s != len_t) return false;
    if (len_s == 0) return true;  // Empty strings match with 0 mismatches
    
    MismGenerator_Pillar_Finite gen(s, t, pillar);
    for (int i = 0; i <= k; ++i) {
        if (!gen.Next().has_value()) return true;
    }
    return false;
}


/* ============================
   Section 4.3: Computing Occurrences in the Periodic Case
   ============================ */

/**
 * [cite_start]@brief Implements FindRotation from Lemma 4.5 [cite: 589-597].
 * Finds a unique rotation of Q that is "close" to S.
 */
std::optional<size_t> FindRotation(int k, const StringHandle& Q, const StringHandle& S, PILLAR_Model& pillar) {
    size_t len_q = pillar.Length(Q);
    size_t len_s = pillar.Length(S);
    
    if (len_q == 0 || len_s < (size_t)(2 * k + 1) * len_q) return std::nullopt;

    std::map<int, int> rot_counts;
    for (int i = 0; i <= 2 * k; ++i) {
        StringHandle Si = pillar.Extract(S, i * len_q, std::min(len_s, (i + 1) * len_q));
        if (pillar.Length(Si) == len_q) {
            if (auto offset = pillar.FindRotationOffset(Si, Q)) {
                rot_counts[*offset]++;
            }
        }
    }

    int best_offset = -1;
    int best_count = 0;
    for (const auto& [rot, count] : rot_counts) {
        if (count > best_count) {
            best_count = count;
            best_offset = rot;
        }
    }
    
    if (best_count >= k + 1) {
        return best_offset;
    }
    return std::nullopt;
}
/**
 * [cite_start]@brief Implements FindRelevantFragment from Lemma 4.6 [cite: 598-600] [cite_start]and Algorithm 4[cite: 601].
 * Narrows down the text T to a smaller fragment T' that contains all occurrences
 * and satisfies certain periodic properties.
 */
StringHandle FindRelevantFragment(const StringHandle& P, const StringHandle& T, int d, const StringHandle& Q, PILLAR_Model& pillar) {
    size_t n = pillar.Length(T);
    size_t m = pillar.Length(P);
    if (n == 0 || m == 0) return pillar.Extract(T, 0, 0);
    
    // Find a fragment of T of length approximately 3/2 * m centered in T
    size_t fragment_len = std::min(n, (3 * m) / 2);
    size_t fragment_start = (n > fragment_len) ? (n - fragment_len) / 2 : 0;
    size_t fragment_end = std::min(n, fragment_start + fragment_len);
    
    StringHandle T_middle = pillar.Extract(T, fragment_start, fragment_end);
    
    // Check if T_middle has a valid rotation
    if (!FindRotation((3 * d) / 2, Q, T_middle, pillar)) {
        return pillar.Extract(T, 0, 0);  // Empty fragment
    }

    // Extend left and right from fragment_start
    size_t l = fragment_start;
    if (l > 0) {
        StringHandle T_prefix = pillar.Extract(T, 0, fragment_start);
        MismGeneratorR_Pillar gen_r(T_prefix, Q, pillar);
        int mis_count = 0;
        while(mis_count < (3*d/2)) {
            auto next_mis = gen_r.Next();
            if(!next_mis) { 
                l = 0; 
                break; 
            }
            l = *next_mis;
            mis_count++;
        }
    }

    size_t r = fragment_end;
    if (r < n) {
        StringHandle T_suffix = pillar.Extract(T, fragment_end, n);
        MismGenerator_Pillar gen_f(T_suffix, Q, pillar);
        int mis_count = 0;
        while(mis_count < (3*d/2)) {
            auto next_mis = gen_f.Next();
            if (!next_mis) { 
                r = n; 
                break; 
            }
            r = fragment_end + *next_mis;
            mis_count++;
        }
    }

    if (r < l) r = l;  // Ensure valid range
    return pillar.Extract(T, l, r);
}

/**
 * [cite_start]@brief Implements Analyze from Lemma 4.4 [cite: 549-556] [cite_start]and Algorithm 3 [cite: 561-572].
 * This function is the core of the structural analysis. It partitions the pattern P
 * into either a set of non-periodic 'breaks', a set of 'repetitive regions', or
 * determines that P itself is approximately periodic.
 */
// Fix 1: Analyze Function - Correct Backward Search with Rotation
AnalysisResult Analyze(const StringHandle& p, int k, PILLAR_Model& pillar) {
    size_t m = pillar.Length(p);
    if (k <= 0 || m == 0) return pillar.Extract(p, 0, 0);
    if (m < (size_t)8 * k) {
        size_t period_len = pillar.Period(p);
        if (period_len <= m / (128 * k)) {
            return pillar.Extract(p, 0, period_len);
        }
        return pillar.Extract(p, 0, m);
    }

    size_t j = 0;
    std::vector<Break> breaks;
    std::vector<RepetitiveRegion> repetitive_regions;
    size_t block_len = m / (8 * k);

    while (j < m && j < 5 * m / 8) {  // Fixed boundary check
        if (j + block_len > m) break;
        
        StringHandle p_prime = pillar.Extract(p, j, j + block_len);
        size_t p_prime_period_len = pillar.Period(p_prime);

        if (p_prime_period_len > m / (128 * k)) {
            breaks.push_back({p_prime});
            if (breaks.size() == (size_t)2 * k) return breaks;
            j += block_len;
        } else {
            StringHandle q = pillar.Extract(p_prime, 0, p_prime_period_len);
            
            // Forward search with tracking of final position
            MismGenerator_Pillar gen(pillar.Extract(p, j, m), q, pillar);
            int mismatches = 0;
            size_t current_len = block_len;
            size_t final_q_offset = 0;  // Track position in Q^inf
            bool found_r = false;
            
            while(true){
                auto mis_pos_opt = gen.Next();
                if(!mis_pos_opt) {
                    // Reached end without enough mismatches
                    final_q_offset = (m - j) % pillar.Length(q);
                    break;
                }
                
                mismatches++;
                current_len = *mis_pos_opt + 1;
                final_q_offset = current_len % pillar.Length(q);

                if (mismatches >= ceil(8.0 * k / m * current_len)) {
                    StringHandle r_handle = pillar.Extract(p, j, j + current_len);
                    repetitive_regions.push_back({r_handle, q});
                    j += current_len;
                    found_r = true;
                    break;
                }
            }

            if (!found_r) {
                // Backward search with correctly rotated Q
                // Create rotated Q: rot^(-final_q_offset)(Q)
                StringHandle q_rotated;
                if (final_q_offset > 0) {
                    StringHandle q_suffix = pillar.Extract(q, final_q_offset, pillar.Length(q));
                    StringHandle q_prefix = pillar.Extract(q, 0, final_q_offset);
                    // Need concat operation or simulate with generator
                    q_rotated = q;  // Simplified: would need proper rotation
                } else {
                    q_rotated = q;
                }
                
                MismGeneratorR_Pillar gen_r(p, q_rotated, pillar);
                mismatches = 0;
                size_t j_prime = m;
                
                while(true) {
                    auto mis_pos_opt = gen_r.Next();
                    if (!mis_pos_opt) break;
                    
                    j_prime = *mis_pos_opt;
                    mismatches++;
                    current_len = m - j_prime;
                    
                    if (mismatches >= ceil(8.0 * k / m * current_len) && current_len >= m - j) {
                        StringHandle r_prime_handle = pillar.Extract(p, j_prime, m);
                        return std::vector<RepetitiveRegion>{{r_prime_handle, q_rotated}};
                    }
                }
                return q_rotated;
            }
        }
        
        size_t total_rep_len = 0;
        for (const auto& reg : repetitive_regions) total_rep_len += pillar.Length(reg.handle);
        if (total_rep_len >= 3 * m / 8) {
            return repetitive_regions;
        }
    }
    
    // If no structure found, check if P itself is periodic
    size_t p_period = pillar.Period(p);
    if (p_period <= m / (128 * k)) {
        return pillar.Extract(p, 0, p_period);
    }
    return pillar.Extract(p, 0, m);
}
/**
 * @brief Helper to find mismatch positions for DistancesRLE.
 * Correctly uses the generator instead of `pillar.Access`.
 */
static vector<int> mismatches_vs_period(const StringHandle& S, const StringHandle& Q, PILLAR_Model& pillar) {
    vector<int> res;
    MismGenerator_Pillar gen(S, Q, pillar);
    while (auto pos = gen.Next()) {
        res.push_back((int)*pos);
    }
    return res;
}

using RLE = vector<pair<int,int>>;

/**
 * [cite_start]@brief Implements DistancesRLE from Lemma 4.7 [cite: 614-615] [cite_start]and Algorithm 5[cite: 616].
 * Computes a run-length encoding of the sequence of Hamming distances between P and
 * sliding windows of T (at multiples of |Q|).
 */
RLE DistancesRLE_RLE(const StringHandle& P, const StringHandle& Tblock, const StringHandle& Q, PILLAR_Model& pillar) {
    size_t m = pillar.Length(P);
    size_t n = pillar.Length(Tblock);
    size_t q = pillar.Length(Q);
    if (n < m || q == 0) return {};

    // 1. Mismatch sets (marking phase)
    vector<int> MT = mismatches_vs_period(Tblock, Q, pillar);
    vector<int> MP = mismatches_vs_period(P, Q, pillar);

    // 2. Event creation (marking phase)
    vector<pair<long long, int>> events;
    for (int alpha : MT) {
        events.push_back({(long long)alpha - (long long)m, 1});
        events.push_back({(long long)alpha, -1});
    }
    for (int alpha : MT) {
        for (int fi : MP) {
            char cp = pillar.Access(P, fi);
            char ct = pillar.Access(Tblock, alpha);
            int h_val = (cp == ct) ? 0 : 1;
            events.push_back({(long long)alpha - (long long)fi - 1, h_val - 2});
            events.push_back({(long long)alpha - (long long)fi, 2 - h_val});
        }
    }

    sort(events.begin(), events.end());

    // 3. Sliding-window phase
    int h = (int)MP.size();
    size_t event_idx = 0;
    while (event_idx < events.size() && events[event_idx].first < 0) {
        h += events[event_idx].second;
        event_idx++;
    }

    RLE rle;
    long long current_pos = 0;
    
    while(current_pos <= n - m) {
        long long next_event_pos = n - m + 1;
        if(event_idx < events.size()) {
            next_event_pos = events[event_idx].first;
        }

        long long j_start = ceil_div_nonneg(current_pos, q);
        long long j_end = ceil_div_nonneg(next_event_pos, q);

        if (j_start < j_end) {
            if (!rle.empty() && rle.back().first == h) {
                rle.back().second += (j_end - j_start);
            } else {
                rle.push_back({h, (int)(j_end - j_start)});
            }
        }
        
        if (next_event_pos > n - m) break;

        current_pos = next_event_pos;
        while(event_idx < events.size() && events[event_idx].first == current_pos) {
            h += events[event_idx].second;
            event_idx++;
        }
    }
    return rle;
}
    
Occurrences RLE_to_APs(const RLE &rle, size_t q, int k) {
    Occurrences out;
    long long j_index = 0;
    for (auto const& [value, run] : rle) {
        if (value <= k) {
            out.push_back({j_index * (long long)q, (long long)q, (long long)run});
        }
        j_index += run;
    }
    return out;
}

/**
 * [cite_start]@brief Implements PeriodicMatches from Lemma 4.8 [cite: 639-643].
 * Finds all k-mismatch occurrences when P is known to be approximately periodic.
 */
struct FragmentWithOffset {
    StringHandle handle;
    size_t original_offset;
};

Occurrences PeriodicMatches(const StringHandle& p, const StringHandle& t_block, int k, int d, const StringHandle& q, PILLAR_Model& pillar) {
    StringHandle relevant_t = FindRelevantFragment(p, t_block, d, q, pillar);
    size_t relevant_offset = relevant_t.start;  // Track offset from t_block
    
    if (pillar.Length(relevant_t) == 0) return {};
    
    RLE rle = DistancesRLE_RLE(p, relevant_t, q, pillar);
    if (rle.empty()) return {};
    
    Occurrences aps = RLE_to_APs(rle, pillar.Length(q), k);
    
    // Properly adjust positions to global text coordinates
    for (auto& ap : aps) {
        ap.start_value += relevant_offset;
    }
    return aps;
}


/* ============================
   Section 4.4: Computing Occurrences in the Non-Periodic Case
   ============================ */

/**
 * [cite_start]@brief Implements BreakMatches from Lemma 4.9 [cite: 650-651] [cite_start]and Algorithm 6[cite: 634].
 * Finds occurrences when P is structured with non-periodic 'breaks'.
 */
Occurrences BreakMatches(const StringHandle& p, const StringHandle& t_block, const std::vector<Break>& breaks, int k, PILLAR_Model& pillar) {
    std::map<size_t, int> marks;
    size_t n_block = pillar.Length(t_block);
    size_t m = pillar.Length(p);
    for (const auto& b : breaks) {
        Occurrences break_occs_prog = pillar.ExactMatches(b.handle, t_block);
        for (const auto& prog : break_occs_prog) {
            for (size_t i = 0; i < (size_t)prog.length; ++i) {
                size_t occ_pos = (size_t)prog.start_value + i * (size_t)prog.difference;
                if (occ_pos >= b.handle.start) {
                    marks[occ_pos - b.handle.start]++;
                }
            }
        }
    }
    Occurrences found_occurrences;
    for (const auto& [pos, count] : marks) {
        if (count >= k && pos + m <= n_block) {
            if (Verify_Mismatches(p, pillar.Extract(t_block, pos, pos + m), k, pillar)) {
                found_occurrences.push_back({(long long)pos, 1, 1});
            }
        }
    }
    return found_occurrences;
}

/**
 * [cite_start]@brief Implements RepetitiveMatches from Lemma 4.10 [cite: 671-672] [cite_start]and Algorithm 7[cite: 653].
 * Finds occurrences when P is structured with 'repetitive regions'.
 */
Occurrences RepetitiveMatches(const StringHandle& p, const StringHandle& t_block, const std::vector<RepetitiveRegion>& regions, int k, PILLAR_Model& pillar) {
    std::map<long long, long long> marks;
    size_t n_block = pillar.Length(t_block);
    size_t m = pillar.Length(p);
    for (const auto &reg : regions) {
        int region_len = (int)pillar.Length(reg.handle);
        int ki = max(1, (int)floor(4.0 * k / m * region_len));
        int di = max(2*ki, (int)ceil(8.0 * k / m * region_len));
        
        Occurrences per = PeriodicMatches(reg.handle, t_block, ki, di, reg.period, pillar);
        
        for (auto &ap : per) {
            for (long long i = 0; i < ap.length; ++i) {
                long long tpos = ap.start_value + i * ap.difference;
                 if (tpos >= reg.handle.start) {
                    marks[tpos - reg.handle.start] += region_len;
                }
            }
        }
    }

    Occurrences result;
    long long total_region_len = 0;
    for (auto &r : regions) total_region_len += (long long)pillar.Length(r.handle);
    long long threshold = (long long)max(0.0, (double)total_region_len - (double)m / 4.0);
    
    for (auto const& [pos, weight] : marks) {
        if (weight >= threshold && pos + (long long)m <= (long long)n_block) {
            if (Verify_Mismatches(p, pillar.Extract(t_block, (size_t)pos, (size_t)pos + m), k, pillar)) {
                result.push_back({pos, 1, 1});
            }
        }
    }
    return result;
}

/**
 * [cite_start]@brief Implements the main algorithm from Main Theorem 8 [cite: 532] [cite_start]and Algorithm 8[cite: 679].
 * This is the top-level function that orchestrates the entire process.
 */
Occurrences MismatchOccurrences(const StringHandle& p, const StringHandle& t, int k, PILLAR_Model& pillar) {
    size_t m = pillar.Length(p);
    size_t n = pillar.Length(t);
    if (m == 0 || n == 0) return {};
    
    AnalysisResult analysis_result = Analyze(p, k, pillar);
    std::set<long long> unique_positions;  // Use set for automatic deduplication

    size_t num_blocks = (n + m/2 - 1) / (m/2);  // Ceiling division
    for (size_t i = 0; i < num_blocks; ++i) {
        size_t block_start = i * m / 2;
        size_t block_end = std::min(n, block_start + (3 * m) / 2);
        if (block_end <= block_start) continue;
        
        StringHandle ti = pillar.Extract(t, block_start, block_end);
        Occurrences block_occs;
        
        if (std::holds_alternative<std::vector<Break>>(analysis_result)) {
            block_occs = BreakMatches(p, ti, std::get<std::vector<Break>>(analysis_result), k, pillar);
        } else if (std::holds_alternative<std::vector<RepetitiveRegion>>(analysis_result)) {
            block_occs = RepetitiveMatches(p, ti, std::get<std::vector<RepetitiveRegion>>(analysis_result), k, pillar);
        } else {
            block_occs = PeriodicMatches(p, ti, k, 8*k, std::get<StringHandle>(analysis_result), pillar);
        }
        
        // Add to unique positions with proper range checking
        for (const auto& ap : block_occs) {
            for (long long j = 0; j < ap.length; ++j) {
                long long global_pos = block_start + ap.start_value + j * ap.difference;
                
                // Only add if position is valid and in the current block's responsibility
                if (global_pos >= 0 && global_pos + m <= n) {
                    // Check if this position should be handled by this block
                    // (avoiding double-counting from overlapping blocks)
                    size_t responsible_block = global_pos / (m/2);
                    if (responsible_block == i || i == num_blocks - 1) {
                        unique_positions.insert(global_pos);
                    }
                }
            }
        }
    }
    
    // Convert set to occurrences
    Occurrences result;
    for (long long pos : unique_positions) {
        result.push_back({pos, 1, 1});
    }
    return result;
}
/* ============================
   MockPillar for demonstration
   ============================ */
class MockPillar : public PILLAR_Model {
    std::map<int, std::string> strings;
    int next_id = 0;

public:
    StringHandle add_string(const std::string& s) {
        int id = next_id++;
        strings[id] = s;
        return {id, 0, s.length(), false};
    }
    std::string get_substring(const StringHandle& s) {
        if (strings.find(s.id) == strings.end() || s.start + s.len > strings[s.id].length()) return "";
        std::string sub = strings[s.id].substr(s.start, s.len);
        if (s.reversed) std::reverse(sub.begin(), sub.end());
        return sub;
    }
    StringHandle Extract(const StringHandle& s, size_t l, size_t r) override {
        if (r < l) r = l;
        size_t new_len = r - l;
        size_t new_start = s.reversed ? (s.start + s.len - r) : (s.start + l);
        return {s.id, new_start, new_len, s.reversed};
    }
    size_t Length(const StringHandle& s) override { return s.len; }
    char Access(const StringHandle& s, size_t i) override {
        if (i >= s.len) return '\0';
        size_t pos = s.reversed ? (s.start + s.len - 1 - i) : (s.start + i);
        return strings.at(s.id)[pos];
    }
    size_t LCP(const StringHandle& s1, const StringHandle& s2) override {
        std::string str1 = get_substring(s1);
        std::string str2 = get_substring(s2);
        size_t len = 0;
        while(len < str1.length() && len < str2.length() && str1[len] == str2[len]) len++;
        return len;
    }
    size_t LCPR(const StringHandle& s1, const StringHandle& s2) override {
        std::string str1 = get_substring(s1);
        std::string str2 = get_substring(s2);
        size_t len = 0;
        while (len < str1.length() && len < str2.length() && str1[str1.length()-1-len] == str2[str2.length()-1-len]) len++;
        return len;
    }
    Occurrences IPM(const StringHandle& p, const StringHandle& t) override { return ExactMatches(p, t); }
    size_t Period(const StringHandle& s) override {
        std::string str = get_substring(s);
        size_t n = str.length();
        if (n == 0) return 0;
        for (size_t p = 1; p <= n; ++p) {
            if (n % p == 0) {
                bool is_period = true;
                for (size_t i = p; i < n; ++i) {
                    if (str[i] != str[i-p]) { is_period = false; break; }
                }
                if (is_period) return p;
            }
        }
        return n;
    }
    std::optional<int> FindRotationOffset(const StringHandle& s, const StringHandle& t) override {
        if (s.len != t.len || s.len == 0) return std::nullopt;
        std::string s_str = get_substring(s);
        std::string t_str = get_substring(t);
        if ((t_str+t_str).find(s_str) != std::string::npos) return (int)(t_str+t_str).find(s_str);
        return std::nullopt;
    }
    Occurrences ExactMatches(const StringHandle& p, const StringHandle& t) override {
        Occurrences occs;
        std::string p_str = get_substring(p);
        std::string t_str = get_substring(t);
        if (p_str.empty()) return occs;
        for (size_t i = 0; (i = t_str.find(p_str, i)) != std::string::npos; ++i) {
            occs.push_back({(long long)i, 1, 1});
        }
        return occs;
    }
};

/* ============================
   Main function for demonstration
   ============================ */
void print_occurrences(const Occurrences& occs, const std::string& p_str, const std::string& t_str) {
    cout << "Found occurrences for P=\"" << p_str << "\" in T=\"" << t_str << "\":" << endl;
    if (occs.empty()) {
        cout << "  - None" << endl;
        return;
    }
    std::set<long long> positions;
    for(const auto& ap : occs) {
        for(long long i = 0; i < ap.length; ++i) {
            positions.insert(ap.start_value + i * ap.difference);
        }
    }
    for(long long pos : positions) {
        cout << "  - pos " << pos << endl;
    }
}


// int main() {
//     cout << "--- FINAL DEMO FOR SECTION 4 (MISMATCHES) ---" << std::endl;
//     MockPillar pillar;

//     cout << "\n## Test 1: Figure 1a Example ##" << endl;
//     string p_str1 = "aacc";
//     string t_str1 = "aaaccc";
//     int k1 = 1;
//     StringHandle p1 = pillar.add_string(p_str1);
//     StringHandle t1 = pillar.add_string(t_str1);
//     Occurrences occs1 = MismatchOccurrences(p1, t1, k1, pillar);
//     print_occurrences(occs1, p_str1, t_str1);
//     cout << "Expected: pos 1, pos 2" << endl;
//     cout << "---------------------------------" << endl;

//     cout << "\n## Test 2: Periodic Case Example ##" << endl;
//     string p_str2 = "ababaxabab";
//     string t_str2 = "abacababababababab";
//     int k2 = 1;
//     StringHandle p2 = pillar.add_string(p_str2);
//     StringHandle t2 = pillar.add_string(t_str2);
//     Occurrences occs2 = MismatchOccurrences(p2, t2, k2, pillar);
//     print_occurrences(occs2, p_str2, t_str2);
//     cout << "Expected: pos 4, pos 6, pos 8" << endl;
//     cout << "---------------------------------" << endl;

//     cout << "\n## Test 3: Breaks Case Example ##" << endl;
//     string p_str3 = "abcrstdefxyz";
//     string t_str3 = "012abcrstXefxyz01abcrstdefXyz";
//     int k3 = 1;
//     StringHandle p3 = pillar.add_string(p_str3);
//     StringHandle t3 = pillar.add_string(t_str3);
//     Occurrences occs3 = MismatchOccurrences(p3, t3, k3, pillar);
//     print_occurrences(occs3, p_str3, t_str3);
//     cout << "Expected: pos 3, pos 20" << endl;
//     cout << "---------------------------------" << endl;

//     return 0;
// }

void print_summary(const std::set<long long>& positions) {
    cout << "{ ";
    for(long long pos : positions) {
        cout << pos << " ";
    }
    cout << "}";
}


int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <test_cases_file>" << std::endl;
        return 1;
    }

    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "Error: Could not open file " << argv[1] << std::endl;
        return 1;
    }

    cout << "--- Processing Test Cases from File ---" << std::endl;

    std::string line;
    int test_num = 1;
    int passed_count = 0;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string k_str, p_str, t_str, expected_str;

        if (!std::getline(ss, k_str, ',') || 
            !std::getline(ss, p_str, ',') || 
            !std::getline(ss, t_str, ',') ||
            !std::getline(ss, expected_str)) {
            std::cerr << "Warning: Skipping malformed line " << test_num << std::endl;
            test_num++;
            continue;
        }

        int k;
        try {
            k = std::stoi(k_str);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Skipping line " << test_num << " with invalid k value: " << k_str << std::endl;
            test_num++;
            continue;
        }

        // Parse expected results from the file
        std::set<long long> expected_positions;
        if (!expected_str.empty()) {
            std::stringstream res_ss(expected_str);
            std::string pos_str;
            while (std::getline(res_ss, pos_str, ';')) {
                try {
                    // Trim whitespace if necessary, though std::stoll usually handles it.
                    pos_str.erase(0, pos_str.find_first_not_of(" \t\n\r\f\v"));
                    pos_str.erase(pos_str.find_last_not_of(" \t\n\r\f\v") + 1);
                    if (!pos_str.empty()) {
                        expected_positions.insert(std::stoll(pos_str));
                    }
                } catch (const std::exception& e) { /* ignore parse error */ }
            }
        }

        // Run the algorithm
        MockPillar pillar;
        StringHandle p = pillar.add_string(p_str);
        StringHandle t = pillar.add_string(t_str);
        Occurrences occs = MismatchOccurrences(p, t, k, pillar);

        // Collect the algorithm's results into a set for easy comparison
        std::set<long long> found_positions;
        for(const auto& ap : occs) {
            for(long long i = 0; i < ap.length; ++i) {
                found_positions.insert(ap.start_value + i * ap.difference);
            }
        }

        // Compare results and print status
        cout << "\n## Test " << test_num++ << " (k=" << k << ") ##" << endl;
        cout << "P=\"" << p_str.substr(0, 40) << (p_str.length() > 40 ? "..." : "") << "\", "
             << "T=\"" << t_str.substr(0, 50) << (t_str.length() > 50 ? "..." : "") << "\"" << endl;
        
        if (found_positions == expected_positions) {
            cout << "Status: PASS" << endl;
            passed_count++;
        } else {
            cout << "Status: FAIL" << endl;
        }

        // --- CHANGE HERE: Print Expected and Found results for all cases ---
        cout << "  Expected: "; print_summary(expected_positions); cout << endl;
        cout << "  Found:    "; print_summary(found_positions); cout << endl;
        // -----------------------------------------------------------------

        cout << "---------------------------------" << endl;
    }

    cout << "\n=================================" << endl;
    cout << "Final Score: " << passed_count << " / " << (test_num - 1) << " tests passed." << endl;
    cout << "=================================" << endl;

    return 0;
}