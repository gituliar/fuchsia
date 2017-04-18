.PHONY: dist test

# Test actions

test-maxima:
	env SAGE_PATH="$(CURDIR)" \
		sage -python -munittest -fv test.test_suite_maxima \
            2>&1 | tee test/test_suite_maxima.log

test-maple:
	env SAGE_PATH="$(CURDIR)" \
		sage -python -munittest -fv test.test_suite_maple \
            2>&1 | tee test/test_suite_maple.log

fuchsia.py: test/__init__.py

test/*.py::
	env SAGE_PATH="$(CURDIR)" \
		sage -python -munittest -fv test.$(basename $(notdir $@))
