#!/bin/bash
# Usage: add_comment.sh "comment"
# Description: Add comment to the README in current directory
# Date: 2023-08-22
# contact: siva.selvanayagam[at]wur.nl

set -u
set -e

print_log(){
    echo "[`date +'%Y-%m-%d %H:%M:%S'`]: $*"
}

if [ $# -eq 0 ]; then
    echo "No comments supplied to add to README"
    exit 1
fi

COMMENT=""
for i in "$@"
do
    if [ -f $i ]; then
        COMMENT=${COMMENT}" FILE:"$i
    elif [ -d $i ]; then
        COMMENT=${COMMENT}" DIR:"$i
    else
        COMMENT=${COMMENT}" "$i
    fi
done
print_log $COMMENT >> README
