#!/bin/bash

g++ -fPIC -shared FEFSG.cpp dllmain.cpp -o FEFSG.o -std=c++11 -I/home/eschwarz/FEBio/ -L/home/eschwarz/FEBio-build/lib -lfebiomech -lfecore
