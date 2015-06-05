#!/bin/bash
echo Cleaning python code in \'`pwd`\' and subdirectories
# split output of find at newlines.
IFS=$'\n'
# send all relevant files to the code cleaner
find *.py *.sh *.cpp *.h *.pxd *.pyx | xargs ./codecleaner.py
