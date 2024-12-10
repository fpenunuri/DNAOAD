#!/bin/bash

#Check if at least one argument (the source file) is provided
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <source_file.f90> [output_executable_name]"
  exit 1
fi

#Use the provided source file and output executable name or default to "a.out"
SOURCE_FILE=$1
OUTPUT_NAME=${2:-a.out}

gfortran -I../LibDualzn128 -o "$OUTPUT_NAME" "$SOURCE_FILE" -L../LibDualzn128 -ldualzn
