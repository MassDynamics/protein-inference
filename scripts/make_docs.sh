#!/bin/bash

cd docs
make html # make docs
make gh-pages # push docs to gh-pages branch

cd ..
git checkout main # switch back to main branch