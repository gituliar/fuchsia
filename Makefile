.PHONY: dist

# Test actions

test:
	# This is the same as `test-maxima`
	sage -python setup.py test

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

# PyPI actions

dist:
	sage -python setup.py sdist

register:
	sage -python setup.py register

upload-pypi:
	sage -python setup.py sdist upload

upload-http:
	@echo "Not implemented"

clean:
	@rm -fr build dist fuchsia.egg-info
	@rm -fr *.pyc */*.pyc

