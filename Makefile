install:
	python3 -m pip install -r requirements.txt

test:
	python3 -m pytest -v --flake8 --pylint --pylint-rcfile=./.pylintrc --mypy dif.py tests/test_dif.py