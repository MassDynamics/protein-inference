#!/bin/bash

docker run -id --name protein_inference \
	--mount type=bind,source=$HOME/data,destination=/data \
	--mount type=bind,source="$(pwd)",destination=/protein_inference \
	proteininference:latest 
