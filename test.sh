#!/usr/local/bash

rm -rf datasets
mkdir datasets
./createdataset.py

motif/findmotifs.py datasets/*

