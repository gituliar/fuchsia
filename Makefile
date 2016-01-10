.PHONY: install test-all test-unit

install:
	@date +"%y.%m.%d" > delirium/VERSION
	@sudo sage-python setup.py develop

test-all: test-func test-unit

test-func: install
	@echo -e "\nRun Functional Tests"
	@sage-python -m pytest ./test/

test-unit: install
	@echo -e "\nRun Unit Tests"
	@sage-python -m pytest --doctest-modules ./delirium/

test:
	sage -python -mdoctest delirium.py
	sage -python test_delirium.py

clean:
	rm -f *.pyc
