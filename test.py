import random
import string
import sys

def hamming_distance(s1, s2):
    """Calculates the Hamming distance between two strings of equal length."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def find_expected_matches(k, pattern, text):
    """Finds all occurrences using a brute-force sliding window."""
    matches = []
    p_len = len(pattern)
    t_len = len(text)
    if p_len == 0 or p_len > t_len:
        return ""
    for i in range(t_len - p_len + 1):
        substring = text[i:i + p_len]
        if hamming_distance(pattern, substring) <= k:
            matches.append(str(i))
    return ";".join(matches)

def generate_random_string(length, alphabet=string.ascii_lowercase):
    """Generates a random string of a given length."""
    return ''.join(random.choice(alphabet) for _ in range(length))

def generate_string_with_errors(original, k, alphabet=string.ascii_lowercase):
    """Introduces k random character substitutions into a string."""
    s_list = list(original)
    if not s_list:
        return ""
    for _ in range(k):
        pos = random.randint(0, len(s_list) - 1)
        s_list[pos] = random.choice(alphabet)
    return "".join(s_list)

def create_test_file(filename, num_cases):
    """Generates and writes test cases to the specified file."""
    print(f"Generating {num_cases} test cases...")
    with open(filename, "w") as f:
        test_cases = []
        # Add some initial, deterministic cases
        test_cases.extend([
            (1, "aacc", "aaaccc"),
            (1, "ababaxabab", "abacababababababab"),
            (1, "abcrstdefxyz", "012abcrstXefxyz01abcrstdefXyz"),
            (0, "abc", "abcabcabc"),
            (1, "", "abc"), # Empty pattern
            (2, "short", "a very long text that does not contain the pattern at all"),
            (1, "aba", "abababa") # Overlapping matches
        ])

        # Generate the rest of the cases randomly
        while len(test_cases) < num_cases:
            case_type = random.choice(['periodic', 'random', 'embedded'])
            
            p, t, k = "", "", 0

            if case_type == 'periodic':
                base = generate_random_string(random.randint(2, 4), "ab")
                p = base * random.randint(3, 6)
                t = (base * random.randint(10, 20))
                k = random.randint(1, 3)
                t = generate_string_with_errors(t, k * 5, "abc")
            
            elif case_type == 'embedded':
                p_len = random.randint(8, 15)
                t_len = random.randint(50, 100)
                p = generate_random_string(p_len)
                t = generate_random_string(t_len)
                k = random.randint(1, p_len // 4)
                
                # Embed the pattern with errors a few times
                for _ in range(random.randint(1, 2)):
                    p_with_errors = generate_string_with_errors(p, random.randint(0, k))
                    if len(t) > p_len:
                        insert_pos = random.randint(0, len(t) - p_len -1)
                        t = t[:insert_pos] + p_with_errors + t[insert_pos+p_len:]
            
            else: # 'random' case
                p_len = random.randint(5, 20)
                t_len = random.randint(30, 80)
                p = generate_random_string(p_len)
                t = generate_random_string(t_len)
                k = random.randint(1, p_len // 3)

            test_cases.append((k, p, t))
            
        # Write all generated cases to the file
        for k_val, p_str, t_str in test_cases:
            # Clean strings of commas to not break the format
            p_clean = p_str.replace(",",";")
            t_clean = t_str.replace(",",";")
            expected_positions = find_expected_matches(k_val, p_clean, t_clean)
            f.write(f"{k_val},{p_clean},{t_clean},{expected_positions}\n")
    print(f"Successfully generated {num_cases} test cases in '{filename}'")

if __name__ == "__main__":
    create_test_file("test_cases.txt", 10000)