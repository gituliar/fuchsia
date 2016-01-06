test:
	sage -python -mdoctest delirium.py
	sage -python test_delirium.py

clean:
	rm -f *.pyc
