cmake_minimum_required(VERSION 3.13)
project(SwarmOptimization)

set(CMAKE_CXX_STANDARD 14)

add_executable(SwarmOptimization main.cpp BenchmarkFunctions.cpp BenchmarkFunctions.h utilities.cpp utilities.h ParticleSwarm.cpp ParticleSwarm.h SwarmUtilities.cpp SwarmUtilities.h FireflyAlgorithm.cpp FireflyAlgorithm.h HarmonySearch.cpp HarmonySearch.h)