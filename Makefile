# Simple Makefile to build lcp_acl, lcp_libsais, and lcp_sdsl

CXX := g++
CXXFLAGS := -std=c++17 -O2

# Paths
ACL_INC := ac-library
SDSL_INC := /u20/hasibih/include
SDSL_LIB := /u20/hasibih/lib
LIBSAIS_INC := libsais/include
LIBSAIS_SRC := libsais/src/libsais.c

.PHONY: all clean

all: lcp_acl lcp_libsais lcp_sdsl main

# ACL-based build (AtCoder Library)
lcp_acl: lcp_acl.cpp
	$(CXX) $(CXXFLAGS) -I $(ACL_INC) $< -o $@

# LibSAIS-based build (compile libsais.c directly)
lcp_libsais: lcp_libsais.cpp $(LIBSAIS_SRC)
	$(CXX) $(CXXFLAGS) -I $(LIBSAIS_INC) $^ -o $@

# SDSL-based build
lcp_sdsl: lcp_sdsl.cpp
	$(CXX) $(CXXFLAGS) -I $(SDSL_INC) $< -L $(SDSL_LIB) -lsdsl -ldivsufsort -ldivsufsort64 -o $@

# PILLAR Hamming Distance with LibSAIS
main: main.cpp $(LIBSAIS_SRC)
	$(CXX) $(CXXFLAGS) -I $(LIBSAIS_INC) $^ -o $@

clean:
	rm -f lcp_acl lcp_libsais lcp_sdsl main

