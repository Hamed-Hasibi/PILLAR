## PILLAR: k-mismatch pattern matching (section 4 demo)

This repository contains a C++17 implementation and demo harness for computing k-mismatch occurrences of a pattern `P` in a text `T` using the PILLAR model (P. Charalampopoulos, T. Kociumaka and P. Wellnitz, "Faster Approximate Pattern Matching: A Unified Approach," 2020 IEEE 61st Annual Symposium on Foundations of Computer Science (FOCS) (Section 4). Two executable variants are provided:

- `main` (LibSAIS-based LCP) — uses `libsais` for efficient suffix/LCP functionality within a mocked PILLAR interface.
- `main2` (standalone ver2) — a self-contained variant not requiring LibSAIS for the demo.

The `test.py` helper generates large test suites and expected answers using a brute-force verifier, enabling quick validation of the algorithm’s results.

### Repository layout
- `main.cpp`: PILLAR mismatches demo using LibSAIS-backed LCP and a `MockPillar` for string ops. Reads test cases from file.
- `main_ver2.cpp`: Alternative standalone demo (no LibSAIS requirement to build).
- `Makefile`: Build targets for the demos and optional LCP experiments.
- `test.py`: Generates `test_cases.txt` with expected results for regression testing.
- `ac-library/`, `libsais/`: Optional third-party sources used by some targets.

### Prerequisites
- g++ with C++17 support
- Make
- Python 3 (for generating test cases)
- Optional (only for some targets):
  - `libsais` available at `libsais/` with `include/` and `src/libsais.c` (used by `main` and `lcp_libsais`)
  - `ac-library` at `ac-library/` (used by `lcp_acl`)
  - SDSL and divsufsort libraries installed and discoverable (used by `lcp_sdsl`). The current `Makefile` expects headers at `/u20/hasibih/include` and libs at `/u20/hasibih/lib`.

### Build
To build everything:

```bash
make all
```

This produces the following binaries in the repo root:
- `main`
- `main2`
- `lcp_acl` (optional, requires `ac-library/`)
- `lcp_libsais` (optional, requires `libsais/`)
- `lcp_sdsl` (optional, requires `SDSL` and `divsufsort`)

Clean builds:
```bash
make clean
```

### Generate test cases
`main` and `main2` expect a CSV-like test file with lines:
```
<k>,<pattern>,<text>,<expected_positions>
```
- `<expected_positions>` is a semicolon-separated list of starting indices, e.g. `0;3;7`.

You can generate a large suite with the included script:
```bash
python3 test.py
```
This creates `test_cases.txt` in the repo root (10,000 cases by default).

### Run
Run either executable with the path to the test file:

The program will:
- Parse each line as a test case: `k,pattern,text,expected_positions`.
- Compute all k-mismatch occurrences using the PILLAR Section 4 pipeline.
- Report PASS/FAIL per test and a final score summary.

### Input format example
```
1,aacc,aaaccc,0;1;2
1,ababaxabab,abacababababababab,4;6;8
1,abcrstdefxyz,012abcrstXefxyz01abcrstdefXyz,3;17
```

### Optional LCP experiments
The `Makefile` also provides three LCP-related targets:
- `lcp_acl`: builds from `lcp_acl.cpp` using AtCoder Library headers under `ac-library/`.
- `lcp_libsais`: builds from `lcp_libsais.cpp` and `libsais/src/libsais.c`.
- `lcp_sdsl`: builds from `lcp_sdsl.cpp` linking against SDSL and divsufsort (paths configurable at the top of the `Makefile`).

### Troubleshooting
- SDSL not found: adjust `SDSL_INC` and `SDSL_LIB` at the top of the `Makefile` to point to your installation, then rebuild `lcp_sdsl`.
- `libsais` missing: ensure `libsais/include` and `libsais/src/libsais.c` exist (e.g., add the submodule or copy sources) for building `main` and `lcp_libsais`.
- `ac-library` missing: clone or copy AtCoder Library headers into `ac-library/` to build `lcp_acl`.
- Build flags: you can tweak optimization or standards via `CXXFLAGS` in the `Makefile`.

### Notes
- The demos use a `MockPillar` to expose PILLAR operations over in-memory strings while faithfully following the paper’s Section 4 structure (Analyze, Periodic/Break/Repetitive cases, and verification).
- `main` and `main2` print both expected and found positions for every test, aiding debugging.

### License
