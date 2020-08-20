/**
 * @file utilities.cpp
 * @author  Al Timofeyev
 * @date    April 15, 2019
 * @brief   This utilities file is used as a helper file for ProcessFunctions.h
 *          and SearchAlgorithms.h, and to create matricies using the Mersenne Twister.
 */

#include "utilities.h"

// Utility Functions.
// --------------------------------------------------------------------------------------------------------------------
/**
 * @brief Parses a string of numbers into a vector of doubles.
 *
 * Constructs and returns a vector of doubles, given a string list of
 * numbers and a delimiter.
 *
 * @note The input string str MUST be a list of doubles!
 *
 * @param str  A string list of numbers.
 * @param delimiter A string of character(s) used to separate the numbers in the string list.
 *
 * @return Returns a vector filled with doubles that were extracted from the string list.
 */
vector<double> parseStringDbl(string str, string delimiter)
{
    // Copy the string str over to a char[] array.
    char list[str.size()+1];
    strcpy(list, str.c_str());

    // Copy the delimiter into a char[] array.
    char delims[delimiter.size()+1];
    strcpy(delims, delimiter.c_str());

    // Set up a vector where to store the numbers.
    vector<double> numList;

    // Set up a tokenizing variable.
    char * token;
    token = strtok (list,delims);

    // While there are still tokens
    while(token != 0)
    {
        if(strcmp(token, "pi") == 0)
            numList.push_back(M_PI);
        else
            numList.push_back(stod(token));
        token = strtok(NULL, delims);
    }

    // Return the vector of doubles.
    return numList;
}

/**
 * @brief Parses a string of numbers into a vector of integers.
 *
 * Constructs and returns a vector of integers, given a string list of
 * numbers and a delimiter.
 *
 * @note The input string list MUST be a list of integers!
 *
 * @param str  A string list of numbers.
 * @param delimiter A string of character(s) used to separate the numbers in the string list.
 *
 * @return Returns a vector filled with integers that were extracted from the string list.
 */
vector<int> parseStringInt(string str, string delimiter)
{
    // Copy the string str over to a char[] array.
    char list[str.size()+1];
    strcpy(list, str.c_str());

    // Copy the delimiter into a char[] array.
    char delims[delimiter.size()+1];
    strcpy(delims, delimiter.c_str());

    // Set up a vector where to store the numbers.
    vector<int> numList;

    // Set up a tokenizing variable.
    char * token;
    token = strtok (list,delims);

    // While there are still tokens
    while(token != 0)
    {
        if(strcmp(token, "pi") == 0)
            numList.push_back((int)M_PI);
        else
            numList.push_back(stoi(token));
        token = strtok(NULL, delims);
    }

    // Return the vector of integers.
    return numList;
}

/**
 * @brief Parses a string of elements into a vector of strings.
 *
 * Constructs and returns a vector of strings, given a string list of
 * elements and a delimiter.
 *
 * @param str  A string list of characters.
 * @param delimiter A string of character(s) used to separate the numbers in the string list.
 *
 * @return Returns a vector filled with integers that were extracted from the string list.
 */
vector<string> parseStringStr(string str, string delimiter)
{
    // Copy the string str over to a char[] array.
    char list[str.size()+1];
    strcpy(list, str.c_str());

    // Copy the delimiter into a char[] array.
    char delims[delimiter.size()+1];
    strcpy(delims, delimiter.c_str());

    // Set up a vector where to store the string elements.
    vector<string> strList;

    // Set up a tokenizing variable.
    char * token;
    token = strtok (list,delims);

    // While there are still tokens
    while(token != 0)
    {
        strList.push_back(token);
        token = strtok(NULL, delims);
    }

    // Return the vector of strings.
    return strList;
}


/**
 * @brief Resizes the vector to size 3.
 *
 * Resizes the given vector to size three in order to prep it for the
 * matrix of a function. Because to generate a matrix, you only need
 * 3 values: function ID, minimum bound, maximum bound.
 *
 * @param setup The vector that's going to be resized for the matrix setup.
 */
void prepForFunctionMatrix(vector<double> &setup)
{
    // Resize the vector to size 3.
    setup.resize(3);
}