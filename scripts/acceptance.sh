#!/bin/bash

set -e -x

current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

python3 -B -m unittest ${current_dir}/../tests/acceptance/*_test.py