/**
 * @file SwarmUtilities.h
 * @author  Al Timofeyev
 * @date    May 10, 2019
 * @brief   Utilities library for swarm algorithms.
 */

#ifndef SWARMOPTIMIZATION_SWARMUTILITIES_H
#define SWARMOPTIMIZATION_SWARMUTILITIES_H


#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "BenchmarkFunctions.h"

using namespace std;


/** Prints all the possible Function IDs to the screen. */
void printAllFunctionIDs();

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------
/** Creates a matrix with the given min/max bound for the given number of rows/columns.*/
vector<vector<double>> createMatrix(int rows, int columns, double minBound, double maxBound);
/** Creates a matrix with the given min/max bound for the given number of rows/columns.*/
vector<vector<double>> createMatrixMT(int rows, int columns, double minBound, double maxBound, mt19937 &randGenerator);

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------
/** Calculates the fitness of a single vector.*/
double calculateFitnessOfVector(vector<double> &vect, int functionID);
/** Calculates the fitness of all vectors in matrix.*/
vector<double> calculateFitnessOfMatrix(vector<vector<double>> matrix, int functionID);

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------
/** Calculates the average value of a vector of doubles.*/
double calculateAverage(vector<double> vect);

/** Calculates the standard deviation value of a vector of doubles.*/
double calculateStandardDeviation(vector<double> vect);

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------
/** Special Quicksort implementation for fitness/matrices.*/
void quicksort(vector<double> &fitnessList, vector<vector<double>> &matrix, int L, int R);
/** Swap function for the Quicksort.*/
void swap(vector<double> &fitnessList, vector<vector<double>> &matrix, int x, int y);

// ------------------------------------------------------------------------------------------------------
/** Normal Quicksort implementation for vector arrays.*/
void quicksort(vector<double> &vec, int L, int R);
void swap(vector<double> &v, int x, int y);

#endif //SWARMOPTIMIZATION_SWARMUTILITIES_H
