// file: pillar_hamming_section4_final_commented.cpp
// g++ -std=c++17 -O2 -I./libsais/include main.cpp ./libsais/src/libsais.c -o pillar_hamming_with_libsais

#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <variant>
#include <optional>
#include <cmath>
#include <memory>
#include <bits/stdc++.h>
extern "C" {
    #include "libsais.h"
}

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
   LibSAIS-based LCP Implementation
   ============================ */

// Include the LCP implementation from separate header file
#include "lcp_libsais.h"

/* ============================
   PILLAR model interface [cite: 330-345]
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

    virtual StringHandle Rotate(const StringHandle& s, size_t amount) = 0;
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
 * @brief Implements LCP(S, Q^inf) logic, as per Corollary 2.9 [cite: 355-358].
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
            break;
        }
        q_pos = 0;
    }
    return total_lcp;
}

/**
 * @brief Implements LCPR(S, Q^inf) logic, as per Corollary 2.10 [cite: 359-360].
 */
size_t LCPR_vs_Period(const StringHandle& s, const StringHandle& q_period, PILLAR_Model& pillar) {
    return LCP_vs_Period(reverse_handle(s), reverse_handle(q_period), pillar);
}


/* ============================
   Section 4.1: Auxiliary PILLAR Model Operations
   ============================ */

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
 * @brief Implements Verify(S, T, k) from Lemma 4.3 [cite: 546-548].
 */
bool Verify_Mismatches(const StringHandle& s, const StringHandle& t, int k, PILLAR_Model& pillar) {
    if (pillar.Length(s) != pillar.Length(t)) return false;
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
 * @brief Implements FindRotation from Lemma 4.5 [cite: 589-597].
 */
std::optional<StringHandle> FindRotation(int k, const StringHandle& Q, const StringHandle& S, PILLAR_Model& pillar) {
    size_t len_q = pillar.Length(Q);
    size_t len_s = pillar.Length(S);
    if (len_s < (size_t)(2 * k + 1) * len_q) return std::nullopt;

    std::map<int, int> rot_counts;
    for (int i = 0; i <= 2 * k; ++i) {
        StringHandle Si = pillar.Extract(S, i * len_q, (i + 1) * len_q);
        if (auto offset = pillar.FindRotationOffset(Si, Q)) {
            rot_counts[*offset]++;
        }
    }

    for (auto const& [rot, count] : rot_counts) {
        if (count >= k + 1) {
            return pillar.Rotate(Q, rot);
        }
    }
    return std::nullopt;
}

/**
 * @brief Implements FindRelevantFragment from Lemma 4.6 [cite: 598-600] and Algorithm 4 [cite: 601-605].
 */
std::pair<StringHandle, size_t> FindRelevantFragment(const StringHandle& T, const StringHandle& P, int d, const StringHandle& Q, PILLAR_Model& pillar) {
    size_t n = pillar.Length(T);
    size_t m = pillar.Length(P);
    if (n < m) return {pillar.Extract(T, 0, 0), 0};

    StringHandle T_middle = pillar.Extract(T, n - m, m);
    auto rotated_Q_opt = FindRotation(ceil(1.5 * d), Q, T_middle, pillar);

    if (!rotated_Q_opt) {
        return {pillar.Extract(T, 0, 0), 0};
    }
    StringHandle rotated_Q = *rotated_Q_opt;

    size_t l = n - m;
    size_t r = m;
    
    MismGenerator_Pillar gen_r(pillar.Extract(T, r, n), rotated_Q, pillar);
    int mis_count_r = 0;
    while(mis_count_r <= (int)ceil(1.5 * d)) {
        auto next_mis = gen_r.Next();
        if(!next_mis) { r = n; break; }
        r = r + *next_mis + 1;
        mis_count_r++;
    }

    MismGeneratorR_Pillar gen_l(pillar.Extract(T, 0, l), rotated_Q, pillar);
    int mis_count_l = 0;
    while(mis_count_l <= (int)ceil(1.5 * d)) {
        auto next_mis = gen_l.Next();
        if(!next_mis) { l = 0; break; }
        l = *next_mis;
        mis_count_l++;
    }

    if (r < l) r = l;
    return {pillar.Extract(T, l, r), l};
}

/**
 * @brief Implements Analyze from Lemma 4.4 [cite: 549-556] and Algorithm 3 [cite: 561-572].
 */
AnalysisResult Analyze(const StringHandle& p, int k, PILLAR_Model& pillar) {
    size_t m = pillar.Length(p);
    if (k <= 0 || m == 0) return pillar.Extract(p, 0, 0);
    
    // FIX: Correctly handle short patterns where block analysis is invalid.
    if (m < (size_t)8 * k) {
        // For these patterns, the block-based analysis doesn't apply.
        // Treat the whole pattern as a single break, which is a robust fallback.
        return std::vector<Break>{{p}};
    }

    size_t j = 0;
    std::vector<Break> breaks;
    std::vector<RepetitiveRegion> repetitive_regions;
    
    size_t block_len = m / (8 * k);

    while (true) {
        size_t total_processed_len = j;
        if (total_processed_len + block_len > m || total_processed_len >= 5.0/8.0 * m) break;

        StringHandle p_prime = pillar.Extract(p, j, j + block_len);
        size_t p_prime_period_len = pillar.Period(p_prime);

        if (p_prime_period_len > m / (128.0 * k)) {
            breaks.push_back({p_prime});
            if (breaks.size() == (size_t)2 * k) return breaks;
            j += block_len;
        } else {
            StringHandle q = pillar.Extract(p_prime, 0, p_prime_period_len);
            
            MismGenerator_Pillar gen(pillar.Extract(p, j, m), q, pillar);
            int mismatches = 0;
            size_t current_len = block_len;
            bool found_r = false;
            
            while(true){
                auto mis_pos_opt = gen.Next();
                if(!mis_pos_opt) { 
                    current_len = m - j;
                    break;
                }
                
                mismatches++;
                current_len = *mis_pos_opt + 1;

                if (mismatches >= ceil(8.0 * k / m * current_len)) {
                    StringHandle r_handle = pillar.Extract(p, j, j + current_len);
                    repetitive_regions.push_back({r_handle, q});
                    j += current_len;
                    found_r = true;
                    break;
                }
            }

            if (!found_r) { 
                size_t q_len = pillar.Length(q);
                size_t rot_amount = (q_len > 0) ? (q_len - (current_len % q_len)) % q_len : 0;
                StringHandle q_rotated = pillar.Rotate(q, rot_amount);

                MismGeneratorR_Pillar gen_r(p, q_rotated, pillar);
                mismatches = 0;
                
                while(true) {
                    auto mis_pos_opt = gen_r.Next();
                    if (!mis_pos_opt) break;

                    mismatches++;
                    size_t back_len = m - *mis_pos_opt;
                    
                    if (mismatches >= ceil(8.0 * k / m * back_len) && back_len >= m - j) {
                        StringHandle r_prime_handle = pillar.Extract(p, *mis_pos_opt, m);
                        return std::vector<RepetitiveRegion>{{r_prime_handle, q}};
                    }
                }
                return q;
            }
        }
        
        size_t total_rep_len = 0;
        for (const auto& reg : repetitive_regions) total_rep_len += pillar.Length(reg.handle);
        if (total_rep_len >= (size_t)floor(3.0/8.0 * m)) {
            return repetitive_regions;
        }
    }
    
    return p;
}

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
 * @brief Implements DistancesRLE from Lemma 4.7 [cite: 614-615] and Algorithm 5 [cite: 616-623].
 */
RLE DistancesRLE_RLE(const StringHandle& P, const StringHandle& Tblock, const StringHandle& Q, PILLAR_Model& pillar) {
    size_t m = pillar.Length(P);
    size_t n = pillar.Length(Tblock);
    size_t q = pillar.Length(Q);
    if (n < m || q == 0) return {};

    vector<int> MT = mismatches_vs_period(Tblock, Q, pillar);
    vector<int> MP = mismatches_vs_period(P, Q, pillar);

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

    // FIX: Sort events with a custom comparator for stability at same-position events.
    sort(events.begin(), events.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second < b.second; 
    });

    int h = (int)MP.size();
    size_t event_idx = 0;
    
    while (event_idx < events.size() && events[event_idx].first < 0) {
        h += events[event_idx].second;
        event_idx++;
    }

    RLE rle;
    long long current_pos = 0;
    
    while(current_pos <= (long long)(n - m)) {
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
        
        if (next_event_pos > (long long)(n - m)) break;

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
 * @brief Implements PeriodicMatches from Lemma 4.8 [cite: 639-643].
 */
Occurrences PeriodicMatches(const StringHandle& p, const StringHandle& t_block, int k, int d, const StringHandle& q, PILLAR_Model& pillar) {
    auto [relevant_t, relevant_t_offset] = FindRelevantFragment(t_block, p, d, q, pillar);
    if (pillar.Length(relevant_t) == 0) return {};
    
    RLE rle = DistancesRLE_RLE(p, relevant_t, q, pillar);
    if (rle.empty()) return {};

    Occurrences aps = RLE_to_APs(rle, pillar.Length(q), k);
    
    for (auto& ap : aps) {
        ap.start_value += relevant_t_offset;
    }
    return aps;
}


/* ============================
   Section 4.4: Computing Occurrences in the Non-Periodic Case
   ============================ */

/**
 * @brief Implements BreakMatches from Lemma 4.9 [cite: 650-651] and Algorithm 6 [cite: 634-639].
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
 * @brief Implements RepetitiveMatches from Lemma 4.10 [cite: 671-672] and Algorithm 7 [cite: 653-658].
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
 * @brief Implements the main algorithm from Main Theorem 8 [cite: 532] and Algorithm 8 [cite: 679-686].
 */
Occurrences MismatchOccurrences(const StringHandle& p, const StringHandle& t, int k, PILLAR_Model& pillar) {
    size_t m = pillar.Length(p);
    size_t n = pillar.Length(t);
    if (m == 0) return {};
    
    AnalysisResult analysis_result = Analyze(p, k, pillar);

    Occurrences total_occurrences;
    
    for (size_t i = 0; i * m / 2 < n; ++i) {
        size_t block_start = i * m / 2;
        size_t block_end = std::min(n, block_start + 3 * m / 2);
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
        
        for (auto& ap : block_occs) {
            for (long long j = 0; j < ap.length; ++j) {
                long long pos = ap.start_value + j * ap.difference;
                long long global_pos = pos + block_start;

                if (global_pos < (i + 1) * m / 2) {
                    total_occurrences.push_back({global_pos, 1, 1});
                }
            }
        }
    }
    
    sort(total_occurrences.begin(), total_occurrences.end(), [](auto& a, auto& b){
        return a.start_value < b.start_value;
    });
    total_occurrences.erase(unique(total_occurrences.begin(), total_occurrences.end(), [](auto& a, auto& b){
        return a.start_value == b.start_value;
    }), total_occurrences.end());

    return total_occurrences;
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
        LCPQuery lcp_query(str1, str2);
        return lcp_query.query(0, 0);
    }
    size_t LCPR(const StringHandle& s1, const StringHandle& s2) override {
        std::string str1 = get_substring(s1);
        std::string str2 = get_substring(s2);
        std::reverse(str1.begin(), str1.end());
        std::reverse(str2.begin(), str2.end());
        LCPQuery lcp_query(str1, str2);
        return lcp_query.query(0, 0);
    }
    Occurrences IPM(const StringHandle& p, const StringHandle& t) override { return ExactMatches(p, t); }
    size_t Period(const StringHandle& s) override {
        std::string str = get_substring(s);
        size_t n = str.length();
        if (n == 0) return 0;
        for (size_t p = 1; p <= n / 2; ++p) {
             if (n % p == 0 && str.substr(0, n - p) == str.substr(p)) {
                return p;
            }
        }
        return n;
    }
    std::optional<int> FindRotationOffset(const StringHandle& s, const StringHandle& t) override {
        if (s.len != t.len || s.len == 0) return std::nullopt;
        std::string s_str = get_substring(s);
        std::string t_str = get_substring(t);
        std::string temp = t_str + t_str;
        auto pos = temp.find(s_str);
        if (pos != std::string::npos) return (int)pos;
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
    StringHandle Rotate(const StringHandle& s, size_t amount) override {
        std::string str = get_substring(s);
        if (str.empty()) return s;
        size_t len = str.length();
        size_t rot = amount % len;
        std::string rotated_str = str.substr(len - rot) + str.substr(0, len - rot);
        return add_string(rotated_str);
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
//     cout << "Expected: pos 0, pos 1, pos 2" << endl;
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
//     cout << "Expected: pos 3, pos 17" << endl;
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


