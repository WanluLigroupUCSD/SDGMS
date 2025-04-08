#!/bin/bash
gcc -lstdc++ ModifiedBasinHopping.cpp chemistry.hpp DistanceHeuristics.hpp calculators.hpp Chemistry.cpp DistanceHeuristics.cpp calculators.cpp -lm -o SDGMS
