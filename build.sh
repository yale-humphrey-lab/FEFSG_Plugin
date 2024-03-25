#!/bin/bash

g++ -fPIC -shared FEFSG.cpp dllmain.cpp -o FEFSG.o -std=c++11 -I/home/eschwarz/FEBio-FSG/ -L/home/eschwarz/FEBio-FSG-build/lib -lfebiomech -lfecore
