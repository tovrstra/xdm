#!/bin/bash
echo Cleaning python code in \'`pwd`\' and subdirectories
# send all relevant files to the code cleaner
for ext in .py .sh .cpp .h .pxd .pyx; do
    find . | grep "${ext}$" | xargs ./tools/codecleaner.py
done
