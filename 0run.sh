#!/bin/bash

gfortran -o user.exe test.f90
./user.exe > output.dat
