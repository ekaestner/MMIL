#!/bin/sh
echo "Making the tarball target with TARBALLPATH=\"$1\""
gmake -j tarball TARBALLPATH="$1"
