#!/bin/bash

./build/TP.out 8 10 2 5 5 -0.01000 0 1 0.1 0 1 "config.json"
./build/TP.out 8 10 3 5 5 5 -0.01000 0 0 1 0.1 0 0 1 0 0 0 1 "config.json"
./build/TP.out 8 10 4 5 5 5 5 -0.01000 0 0 5 1 0.1 0 0 1 0 0 0 1 0 0 0 0 0 0 0 "config.json"
# ./build/TP.out 8 1000 2 500 500  0 0 1 1 1 1 "data/"
# ./build/TP.out 8 1000 2 250 750  0 0 1 1 1 1 "data/"
# ./build/TP.out 8 1000 2 750 250  0 0 1 1 1 1 "data/"