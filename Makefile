.PHONY: test
test:
	pytest --doctest-modules -vvv


# Generate expected output for integration tests
# Double check diffs to make sure they are correct!
INTEGRATION_SRC := $(wildcard tests/test_*_integration.py)
INTEGRATION_MOD := $(subst .py,,$(INTEGRATION_SRC))
INTEGRATION_MOD := $(subst /,.,$(INTEGRATION_MOD))
.PHONY: test_generate
test_generate: $(INTEGRATION_SRC)
	for mod in $(INTEGRATION_MOD); do python -m "$$mod"; done


.PHONY: README.md
README.md:
	$(MAKE) -C doc README.md
	cp doc/README.md .
	$(MAKE) -C doc clean

.PHONY: clean
clean:
	find . -depth -name "__pycache__" -type d -exec rm -r {} +
	find . -name "*.pyc" -delete
