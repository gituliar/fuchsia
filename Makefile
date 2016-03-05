.PHONY: clean install test

DTESTS := ${wildcard fuchsia/*.py}
UTESTS := ${wildcard test/test_*.py}
.PHONY: $(DTESTS) $(FTESTS) $(UTESTS)

install:
	@sudo sage -python setup.py develop

test: $(DTESTS) $(FTESTS) $(UTESTS)

$(DTESTS): fuchsia/%.py:
	env SAGE_PATH=$(CURDIR) sage -python -mdoctest fuchsia/$*.py

$(UTESTS): test/%.py:
	env SAGE_PATH=$(CURDIR) sage -python test/$*.py

clean:
	@rm -fr fuchsia/*.pyc fuchsia/__pycache__
