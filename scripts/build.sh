
#!/bin/bash

export PYTHONDONTWRITEBYTECODE=true

set -e -x

current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

coverage run -m unittest ${current_dir}/../tests/unit/*_test.py

coverage report -m

unset PYTHONDONTWRITEBYTECODE