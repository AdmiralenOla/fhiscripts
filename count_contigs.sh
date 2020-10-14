#!/bin/bash

programname=$0

function usage {
  echo "usage: $programname fastafile"
  exit 1
}

if [ $# == 0 ]; then
  usage
fi
for file in "$@"
do
  number=$(cat "$file" | grep -c ">")
  echo -e "$file\t$number"
done

