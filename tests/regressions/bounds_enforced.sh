#!/usr/bin/env bash

output=$("${1}" 2>&1)
status=$?

if [[ "${status}" -ne 1 ]]; then
  echo "No input, the program should return 1, it returned ${status} instead".
  exit 1
fi

echo "${output}" | grep "Usage: "
if [[ "${?}" -ne 0 ]]; then
  echo "The program should have given a help message on stderr and quit with status 1."
  exit 1
fi
