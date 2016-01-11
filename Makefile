.PHONY: clean install test

DTESTS := ${wildcard delirium/*.py}
UTESTS := ${wildcard test/test_*.py}
.PHONY: $(DTESTS) $(UTESTS)

install:
	@date +"%y.%m.%d" > delirium/VERSION
	@sudo sage-python setup.py develop

test: $(DTESTS) $(UTESTS)

$(DTESTS): delirium/%.py:
	sage-python -mdoctest delirium/$*.py

$(UTESTS): test/%.py:
	sage-python test/$*.py

clean:
	@rm -fr delirium/*.pyc delirium/__pycache__
