#!/bin/bash

rm -rf build/*

cmake  -DCMAKE_BUILD_TYPE=Release -B build
cmake --build build
