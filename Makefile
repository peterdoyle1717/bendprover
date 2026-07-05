# override on the command line as needed, e.g.
#   make test CLERS_BIN=~/Dropbox/neo/clers/bin/clers
CLERS_BIN ?= clers
DUMPS ?= /Users/doyle/Dropbox/neo/_parking/bendq_sandbox/dumps

test:
	CLERS_BIN=$(CLERS_BIN) python3 tests/run_tests.py

census:
	ls $(DUMPS)/*.dump > /tmp/bendprover_census_list.txt
	CLERS_BIN=$(CLERS_BIN) python3 bendprover.py --batch /tmp/bendprover_census_list.txt census322

.PHONY: test census
