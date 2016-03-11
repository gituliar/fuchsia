.PHONY: build clean install test

TESTS := ${wildcard test/test_*.py}
.PHONY: fuchsia.py $(TESTS)

install:
	sage -python setup.py develop

build:
	sage -python setup.py build_py

clean:
	@rm -fr build dist fuchsia/*.pyc fuchsia/__pycache__ fuchsia.egg-info


test: fuchsia.py $(TESTS)

fuchsia.py:
	sage -python -mdoctest fuchsia.py

$(TESTS): test/%.py:
	sage -python test/$*.py
