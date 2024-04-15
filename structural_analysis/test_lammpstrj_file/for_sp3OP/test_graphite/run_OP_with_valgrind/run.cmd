#!/bin/bash

code='/usr/workspace/wsa/lyu1/utils/structural_analysis/OP'

valgrind --tool=memcheck --log-file=val.log $code graphite_288.lammpstrj 1.7 12 12 1.0 1 > test.out
