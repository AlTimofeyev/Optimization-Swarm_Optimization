/**
 * @file ParticleSwarm.h
 * @author  Al Timofeyev
 * @date    May 10, 2019
 * @brief   This is an implementation of a Particle Swarm.
 */

#ifndef SWARMOPTIMIZATION_PARTICLESWARM_H
#define SWARMOPTIMIZATION_PARTICLESWARM_H

#include <fstream>
#include <chrono>
#include "SwarmUtilities.h"

/**
 * @brief Holds all the user defined variables.
 * Particle Swarm Configuration Structure, where all user defined
 * variables that are used to configure the Particle Swarm are stored.
 */
struct PS_Config
{
    int dimensions;     /**< Number of dimensions per individual in population. */
    int popSize;        /**< Population size. */
    int iterations;     /**< Maximum number of iterations. */
    double k;           /**< Dampening factor for the velocity. */
    double c1;          /**< Scaling factor to bring velocity closer to personal best. */
    double c2;          /**< Scaling factor to bring velocity closer to global best. */
};

/**
 * @brief Holds all the population information.
 * Particle Swarm Population Structure, holds all the data related
 * to the population of the Particle Swarm.
 */
struct PS_Population
{
    int functionID;                     /**< The ID determines which benchmark function to call. */
    vector<double> bounds;              /**< Holds the (min,max) bounds of the values for each individual in the population. */
    int functionCounter = 0;            /**< The function counter keeps track of how many times the benchmark function was called. */
    vector<vector<double>> pop;         /**< The population matrix. */
    vector<vector<double>> velocity;    /**< Velocity matrix holds velocity vectors of each individual in the population. */
    vector<double> fitness;             /**< The fitness for each vector in the population matrix. */
    vector<vector<double>> pBestInd;    /**< The personal best solution for each individual in the population. */
    vector<double> pBestFitness;        /**< The personal best fitness for each individual in the population. */
    vector<double> gBestIndividual;     /**< The global best solution of the population. */
    double gBestFitness;                /**< The global best fitness (of gBestIndividual) of the population. */
    vector<double> gBestFitnessList;    /**< List of the best global fitness values from each iteration. */
    double executionTime = -1.0;        /**< Time(ms) it took to run the Particle Swarm on this population. */
};

/**
 * @brief Particle Swarm Analysis
 * Particle Swarm Analysis Structure, to keep track of the analysis
 * performed on each population in the population list.
 */
struct PS_Analysis
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

class ParticleSwarm
{
public:
    // ---------------------- CONSTRUCTORS ----------------------
    ParticleSwarm(int dimensions, int populationSize, int maxIterations, double kDampeningFactor, double c1, double c2);    /**< The Particle Swarm constructor. */

    // ------------------------- METHODS ------------------------
    double runParticleSwarm(int functionID, double minBound, double maxBound);  /**< Runs the Particle Swarm with set parameters. */
    void analyzePSResults();                                                    /**< Analyzes the results of the Particle Swarm. */

    void printPSResults();                                                      /**< Prints the Results of the Particle Swarm. */
    void printPSAnalysis();                                                     /**< Prints the Analysis of the Particle Swarm Results. */

    void savePSResults();                                                       /**< Saves all Particle Swarm Results to file. */
    void savePSAnalysis();                                                      /**< Saves the Analysis of the Particle Swarm to file. */
    void saveEndingPopulation();                                                /**< Saves the ending population solutions to file. */

private:
    // ------------------------ VARIABLES -----------------------
    PS_Config psConfig;
    vector<PS_Population> popList;
    PS_Analysis psAnalysis;

    // ------------------------- METHODS ------------------------
    void generatePSPopulation(PS_Population &population, mt19937 &randGenerator);                                           /**< Generates the initial population for Particle Swarm. */

    void calculateParticleVelocity(vector<double> &pVelocity, vector<double> particle, vector<double> pBest, vector<double> gBest, const double &k, const double &c1, const double &c2, mt19937 &randGenerator);    /**< Calculates the velocity of a particle in the Particle Swarm. */
    void evaluatePopulation(int functionID, vector<vector<double>> &pop, vector<double> &fitness, int &functionCounter);    /**< Calculates fitness of all solutions in population. */
    void evaluateIndividual(const int &functionID, vector<double> &indiv, double &fitness, int &functionCounter);           /**< Calculate the fitness of an individual solution of the population. */

    void updateParticle(PS_Population &population, int particleIndex, mt19937 &randGenerator);                              /**< Updates a single particle in the population. */
    void iteratePopulation(PS_Population &population, mt19937 &randGenerator);                                              /**< Iterates the Particle Swarm population to the next generation. */
};


#endif //SWARMOPTIMIZATION_PARTICLESWARM_H
