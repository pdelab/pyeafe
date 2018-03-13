#!/bin/sh
# Unstall PyEAFE

echo "Uninstalling PyEAFE..."

if [ "$1" = "user" ]; then
	echo "this one"
	python setup.py install --user --record installed_files_tmp.txt
else
	echo "this one 2"
	python setup.py install --record installed_files_tmp.txt
fi

cat installed_files_tmp.txt | xargs sudo rm -rf
rm -rf installed_files_tmp.txt

echo "Done."


