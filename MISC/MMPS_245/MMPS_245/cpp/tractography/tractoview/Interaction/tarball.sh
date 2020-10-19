#!/bin/sh
MAKE=gmake
if [ $# -lt 1 ]; then
	echo "syntax:   $0 somefilepath"
	echo "optional: $0 somefilepath make.exe"
        exit
fi
if [ $# -eq 2 ]; then
	MAKE=$2
fi
echo "Making the tarball target with TARBALLPATH=\"$1\""
$MAKE -j tarball TARBALLPATH="$1"
