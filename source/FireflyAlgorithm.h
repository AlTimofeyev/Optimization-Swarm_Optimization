/**
 * @file FireflyAlgorithm.h
 * @author  Al Timofeyev
 * @date    May 14, 2019
 * @brief   This is an implementation of a Firefly Algorithm.
 */

#ifndef SWARMOPTIMIZATION_FIREFLYALGORITHM_H
#define SWARMOPTIMIZATION_FIREFLYALGORITHM_H

#include <fstream>
#include <chrono>
#include "SwarmUtilities.h"

/**
 * @brief Holds all the user defined variables.
 * Firefly Algorithm Configuration Structure, where all user defined
 * variables that are used to configure the Firefly Algorithm are stored.
 */
struct FA_Config
{
    int dimensions;     /**< Number of dimensions per individual in population. */
    int popSize;        /**< Population size. */
    int iterations;     /**< Maximum number of iterations. */
    double alpha;       /**< Alpha scaling factor - range [0,1]. */
    double betaMin;     /**< Beta scaling factor. */
    double gamma;       /**< Gamma scaling factor. */
};

/**
 * @brief Holds all the population information.
 * Firefly Algorithm Population Structure, holds all the data related
 * to the population of the Firefly Algorithm.
 */
struct FA_Population
{
    int functionID;                     /**< The ID determines which benchmark function to call. */
    vector<double> bounds;              /**< Holds the (min,max) bounds of the values for each individual in the population. */
    int functionCounter = 0;            /**< The function counter keeps track of how many times the benchmark function was called. */
    vector<vector<double>> pop;         /**< The population matrix. */
    vector<double> fitness;             /**< The fitness for each vector in the population matrix. */
    vector<double> bestGlobFit;         /**< A list of best global fitness values from each iteration. */
    double executionTime = -1.0;        /**< Time(ms) it took to run the Firefly Algorithm on this population. */
};

/**
 * @brief Firefly Algorithm Analysis
 * Firefly Algorithm Analysis Structure, to keep track of the analysis
 * performed on each population in the population list.
 */
struct FA_Analysis
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


class FireflyAlgorithm
{
public:
    // ---------------------- CONSTRUCTORS ----------------------
    FireflyAlgorithm(int dimensions, int populationSize, int maxIterations, double alpha, double beta, double gamma);   /**< The Firefly Algorithm constructor. */

    // ------------------------- METHODS ------------------------
    double runFireflyAlgorithm(int functionID, double minBound, double maxBound);   /**< Runs the Firefly Algorithm with set parameters. */
    void analyzeFAResults();                                                        /**< Analyzes the results of the Firefly Algorithm. */

    void printFAResults();                                                          /**< Prints the Results of the Firefly Algorithm. */
    void printFAAnalysis();                                                         /**< Prints the Analysis of the Firefly Algorithm Results. */

    void saveFAResults();                                                           /**< Saves all Firefly Algorithm Results to file. */
    void saveFAAnalysis();                                                          /**< Saves the Analysis of the Firefly Algorithm to file. */
    void saveEndingPopulation();                                                    /**< Saves the ending population solutions to file. */

private:
    // ------------------------ VARIABLES -----------------------
    FA_Config faConfig;
    vector<FA_Population> popList;
    FA_Analysis faAnalysis;

    // ------------------------- METHODS ------------------------
    void generateFAPopulation(FA_Population &population, mt19937 &randGenerator);       /**< Generates the initial population. */

    void evaluatePopulation(int functionID, vector<vector<double>> &pop, vector<double> &fitness, int &functionCounter);    /**< Calculates fitness of all solutions in population. */
    void evaluateIndividual(const int &functionID, vector<double> &indiv, double &fitness, int &functionCounter);           /**< Calculate the fitness of an individual solution of the population. */

    double calculateDistanceBetweenFireflies(const vector<double> &firefly1, const vector<double> &firefly2);                                                   /**< Calculates the distance between two fireflies. */
    void calculateLightIntensity(const double &ffFit1, const double &ffFit2, const double &r, const double &gamma, double &lIntensity1, double &lIntensity2);   /**< Calculates the light intensity of two fireflies. */
    double calculateAttractiveness(const double &betaMin, const double &gamma, const double &r);                                                                /**< Calculates the attractiveness between two fireflies. */

    void moveFirefly(const vector<double> &brightFF, vector<double> &lessBrightFF, const double &alpha, const double &betaMin, const double &gamma, const double &r, mt19937 &randGenerator);   /**< Moves the less brighter firefly towards the brighter firefly. */
    void moveLessBrightFireflies(const int &currFFIndex, FA_Population &population, const double &alpha, const double &betaMin, const double &gamma, mt19937 &randGenerator);                   /**< Moves all less brighter fireflies towards the current firefly. */
    void iterateFireflies(FA_Population &population, const double &alpha, const double &betaMin, const double &gamma, mt19937 &randGenerator);                                                  /**< Iterates the firefly population to the next generation. */
};


#endif //SWARMOPTIMIZATION_FIREFLYALGORITHM_H
