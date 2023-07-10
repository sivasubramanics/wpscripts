#!/bin/bash

FILE=$1

usage(){
    echo "Usage: $0 <file>"
    exit 1
}

if [ -z "$FILE" ]; then
    usage
fi

if [ -f "$FILE" ]; then
    echo "$FILE:yes"
else 
    echo "$FILE:no"
fi