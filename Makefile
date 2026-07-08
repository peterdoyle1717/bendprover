# mpfr/gmp paths: /opt/homebrew for arm64 macs; plain -lmpfr -lgmp on linux
MPFRFLAGS ?= $(shell [ -d /opt/homebrew/include ] && echo "-I/opt/homebrew/include -L/opt/homebrew/lib")
CLERS_BIN ?= clers

csrc/euclid_lm_mp: csrc/euclid_lm_mp.c
	cc -O2 -Wall $(MPFRFLAGS) -o $@ $< -lmpfr -lgmp

test: csrc/euclid_lm_mp
	python3 tests/selftest_format.py

test-oracle:
	CLERS_BIN=$(CLERS_BIN) python3 tests/run_tests.py

.PHONY: test test-oracle
