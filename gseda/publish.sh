conda activate py38
rm -rf dist
python3 -m build
python3 -m twine upload --repository pypi dist/*