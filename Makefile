.PHONY: clean dist test

TESTS := ${wildcard test/test_*.py}
.PHONY: $(TESTS)

dist:
	sage -python setup.py sdist

clean:
	@rm -fr build dist fuchsia.egg-info


test:
	sage -python setup.py test

$(TESTS): test/%.py:
	env SAGE_PATH="$(CURDIR)" sage -python test/$*.py
