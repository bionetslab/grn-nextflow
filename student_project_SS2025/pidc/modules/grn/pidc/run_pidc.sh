#!/bin/bash
set -e
INPUT_FILE=$1
OUTPUT_FILE=$2
julia runPIDC.jl "$INPUT_FILE" "$OUTPUT_FILE"


