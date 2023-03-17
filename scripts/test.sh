#!/bin/bash

dos2unix test_objects/*
./build/bin/test
unix2dos test_objects/*
./build/bin/test
