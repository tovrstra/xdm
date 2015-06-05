#!/bin/bash
rm -r build
for i in $(find xdm scripts | egrep "\.pyc$|\.py~$|\.pyc~$|\.bak$|\.so$") ; do rm -v ${i}; done
