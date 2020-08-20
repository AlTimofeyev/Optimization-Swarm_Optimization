/**
 * @file HarmonySearch.h
 * @author  Al Timofeyev
 * @date    May 16, 2019
 * @brief   This is an implementation of a Harmony Search.
 */

#ifndef SWARMOPTIMIZATION_HARMONYSEARCH_H
#define SWARMOPTIMIZATION_HARMONYSEARCH_H


#include <fstream>
#include <chrono>
#include "SwarmUtilities.h"

/**
 * @brief Holds all the user defined variables.
 * Harmony Search Configuration Structure, where all user defined
 * variables that are used to configure the Harmony Search are stored.
 */
struct HS_Config
{
    int dimensions;     /**< Number of dimensions per individual in population. */
    int popSize;        /**< Population size. */
    int iterations;     /**< Maximum number of iterations. */
    double HMCR;        /**< Harmony Memory Consideration Rate - range [0,1]. */
    double PAR;         /**< Pitch Adjustment Rate - range [0,1]. */
    double bandwidth;   /**< The Bandwidth. */
};

/**
 * @brief Holds all the population information.
 * Harmony Search Population Structure, holds all the data related
 * to the population of the Harmony Search.
 */
struct HS_Population
{
    int functionID;                     /**< The ID determines which benchmark function to call. */
    vector<double> bounds;              /**< Holds the (min,max) bounds of the values for each individual in the population. */
    int functionCounter = 0;            /**< The function counter keeps track of how many times the benchmark function was called. */
    vector<vector<double>> pop;         /**< The population matrix. */
    vector<double> fitness;             /**< The fitness for each vector in the population matrix. */
    vector<double> bestGlobFit;         /**< A list of best global fitness values from each iteration. */
    vector<vector<double>> worstSol;    /**< Population of the worst solutions from the actual population. */
    vector<double> worstFitness;        /**< A list of worst fitness values from the worstSol population. */
    double executionTime = -1.0;        /**< Time(ms) it took to run the Harmony Search on this population. */
};

/**
 * @brief Harmony Search Analysis
 * Harmony Search Analysis Structure, to keep track of the analysis
 * performed on each population in the population list.
 */
struct HS_Analysis
{
    string header = "Function ID,Average Fitness,Standard Deviation,Range(min),Range(max),Median,Time(ms),Function Calls\n"; /**< Header used when saving the data.*/
    vector<int> functionIDs;                /**< List of function IDs.*/
    vector<double> avgFunctionFitness;      /**< List of the average fitness from function.*/
    vector<double> standardDeviation;       /**< List of standard fitness deviations.*/
    vector<vector<double>> ranges;          /**< List of ranges for each fitness function.*/
    vector<double> medianFunctionFitness;   /**< List of the Median fitness for each function.*/
    vector<double> executionTimes;          /**< List of execution times in ms for all functions.*/
    vector<int> functionCalls;              /**< List of the amount of times a function was called. */
};

/**
 * @brief Harmony Search Analysis Worst
 * Harmony Search Analysis Worst Structure, to keep track of the analysis
 * performed on the list of worst solutions of each population in the population list.
 */
struct HS_Analysis_Worst
{
    string header = "Function ID,Average Fitness,Standard Deviation,Range(min),Range(max),Median,Time(ms),Function Calls\n"; /**< Header used when saving the data.*/
    vector<int> functionIDs;                /**< List of function IDs.*/
    vector<double> avgFunctionFitness;      /**< List of the average fitness from function.*/
    vector<double> standardDeviation;       /**< List of standard fitness deviations.*/
    vector<vector<double>> ranges;          /**< List of ranges for each fitness function.*/
    vector<double> medianFunctionFitness;   /**< List of the Median fitness for each function.*/
    vector<double> executionTimes;          /**< List of execution times in ms for all functions.*/
    vector<int> functionCalls;              /**< List of the amount of times a function was called. */
};


class HarmonySearch
{
public:
    // ---------------------- CONSTRUCTORS ----------------------
    HarmonySearch(int dimensions, int populationSize, int maxIterations, double HMCR, double PAR, double bandwidth);    /**< The Harmony Search constructor. */

    // ------------------------- METHODS ------------------------
    double runHarmonySearch(int functionID, double minBound, double maxBound);  /**< Runs the Harmony Search with set parameters. */
    void analyzeHSResults();                                                    /**< Analyzes the results of the Harmony Search. */
    void analyzeHSWorstResults();                                               /**< Analyses the worst results of the Harmony Search. */

    void printHSResults();                                                      /**< Prints the Results of the Harmony Search. */
    void printHSAnalysis();                                                     /**< Prints the Analysis of the Harmony Search Results. */
    void printHSWorstAnalysis();                                                /**< Prints the Analysis of the worst Harmony Search Results. */

    void saveHSResults();                                                       /**< Saves all Harmony Search Results to file. */
    void saveHSAnalysis();                                                      /**< Saves the Analysis of the Harmony Search to file. */
    void saveHSWorstAnalysis();                                                 /**< Saves the Analysis of the worst Harmony Search Results to file. */

private:
    // ------------------------ VARIABLES -----------------------
    HS_Config hsConfig;
    vector<HS_Population> popList;
    HS_Analysis hsAnalysis;
    HS_Analysis_Worst hsWorstAnalysis;

    // ------------------------- METHODS ------------------------
    void generateHSPopulation(HS_Population &population, mt19937 &randGenerator);       /**< Generates the initial population. */

    void evaluatePopulation(int functionID, vector<vector<double>> &pop, vector<double> &fitness, int &functionCounter);    /**< Calculates fitness of all solutions in population. */
    void evaluateIndividual(const int &functionID, vector<double> &indiv, double &fitness, int &functionCounter);           /**< Calculate the fitness of an individual solution of the population. */

    double chooseRandomHarmonic(const vector<vector<double>> &pop, mt19937 &randGenerator);                                                             /**< Choose a random harmonic (dimension) from the population. */
    double adjustHarmonicPitch(const vector<vector<double>> &pop, const int &dim, const double &PAR, const double &bandwidth, mt19937 &randGenerator);  /**< Adjust the pitch of a random solution's harmonic (dimension). */
    double generateNewRandHarmonic(const double &minBound, const double &maxBound, mt19937 &randGenerator);                                             /**< Generate a new harmonic (dimension) within the bounds. */

    vector<double> generateNewRandHarmony(HS_Population &population, const double &HMCR, const double &PAR, const double &bandwidth, mt19937 &randGenerator);   /**< Generates a new random harmony from existing population. */
    void iterateHarmony(HS_Population &population, const double &HMCR, const double &PAR, const double &bandwidth, mt19937 &randGenerator);                     /**< Iterate the harmonic population. */
};


#endif //SWARMOPTIMIZATION_HARMONYSEARCH_H
