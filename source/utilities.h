/**
 * @file utilities.h
 * @author  Al Timofeyev
 * @date    April 15, 2019
 * @brief   This utilities file is used as a helper file for ProcessFunctions.h
 *          and SearchAlgorithms.h, and to create matricies using the Mersenne Twister.
 */

#ifndef BENCHMARKFUNCTIONS_UTILITIES_H
#define BENCHMARKFUNCTIONS_UTILITIES_H

#define _USE_MATH_DEFINES // Uncomment if cmath constants are desirable, like M_PI.

#include <iostream>
#include <string>
#include <string.h>
#include <vector>
#include <cmath>

using namespace std;

/** Parses a string of numbers into a vector of doubles.*/
vector<double> parseStringDbl(string str, string delimiter);

/** Parses a string of numbers into a vector of integers.*/
vector<int> parseStringInt(string str, string delimiter);

/** Parses a string of characters into a vector of strings.*/
vector<string> parseStringStr(string str, string delimiter);

/** Preps the setup vector for the matrix of a function by resizing to size 3.*/
void prepForFunctionMatrix(vector<double> &setup);

#endif //BENCHMARKFUNCTIONS_UTILITIES_H
