#include <iostream>
#include <fstream>
#include "utilities.h"
#include "ParticleSwarm.h"
#include "FireflyAlgorithm.h"
#include "HarmonySearch.h"

using namespace std;

void test_PSO(const string &configFilename);
void test_FA(const string &configFilename);
void test_HS(const string &configFilename);

int main(int argc, char **argv)
{
    // Default configuration filenames.
    string psConfigFilename = "PSconfig.txt";
    string faConfigFilename = "FAconfig.txt";
    string hsConfigFilename = "HSconfig.txt";

    string configFilename = hsConfigFilename; // Set the name of the current config file.

    // If a configuration filename was provided, reassign the variable to user input.
    if (argc == 2)
        configFilename = argv[1];


    // Test.
    test_PSO(configFilename);
    //test_FA(configFilename);
    //test_HS(configFilename);

    return 0;
}

// ******************************************************************************************************
// ******************************************************************************************************
void test_PSO(const string &configFilename)
{
    // Open the config file.
    ifstream configFile;
    configFile.open(configFilename);
    if(configFile.fail())
    {
        cout << "Failed to open file: " << configFilename << endl;
        cout << "-----------------------------------------------\n";
        cout << "File is either not in the right directory\n";
        cout << "or does not exist.\n";
        cout << "-----------------------------------------------\n";
        cout << "Accepted File Formats: .txt" << endl;
        cout << "-----------------------------------------------\n";
        cout << "******** TERMINATING PROGRAM EXECUTION ********\n\n";
        exit(1);
    }

    // Get the configuration values from file.
    string line;
    vector<double> configParams;
    while(configFile.good())
    {
        getline(configFile, line);
        vector<double> tokens = parseStringDbl(line, ",\r\n\t");

        if(tokens.size() == 1)
            configParams.push_back(tokens[0]);
        else
            break;
    }

    // Create Particle Swarm Object.
    ParticleSwarm particleSwarm(configParams[0], configParams[1], configParams[2], configParams[3], configParams[4], configParams[5]);

    // Reset the file pointer to the beginning of file.
    configFile.clear();
    configFile.seekg(0, ios::beg);

    // Search from beginning of file.
    while(configFile.good())
    {
        getline(configFile, line);
        vector<double> tokens = parseStringDbl(line, ",\r\n\t");

        // If there aren't 3 parameters, skip that line.
        if(tokens.size() != 3)
            continue;

        // Do Particle Swarm.
        particleSwarm.runParticleSwarm(tokens[0], tokens[1], tokens[2]);
    }

    //particleSwarm.printPSResults();

    particleSwarm.analyzePSResults();
    particleSwarm.printPSAnalysis();

    //particleSwarm.savePSResults();
    particleSwarm.savePSAnalysis();
    //particleSwarm.saveEndingPopulation();

    // Close the configuration file.
    configFile.close();
}


// ******************************************************************************************************
// ******************************************************************************************************

void test_FA(const string &configFilename)
{
    // Open the config file.
    ifstream configFile;
    configFile.open(configFilename);
    if(configFile.fail())
    {
        cout << "Failed to open file: " << configFilename << endl;
        cout << "-----------------------------------------------\n";
        cout << "File is either not in the right directory\n";
        cout << "or does not exist.\n";
        cout << "-----------------------------------------------\n";
        cout << "Accepted File Formats: .txt" << endl;
        cout << "-----------------------------------------------\n";
        cout << "******** TERMINATING PROGRAM EXECUTION ********\n\n";
        exit(1);
    }

    // Get the configuration values from file.
    string line;
    vector<double> configParams;
    while(configFile.good())
    {
        getline(configFile, line);
        vector<double> tokens = parseStringDbl(line, ",\r\n\t");

        if(tokens.size() == 1)
            configParams.push_back(tokens[0]);
        else
            break;
    }

    // Create Firefly Algorithm Object.
    FireflyAlgorithm fireflyAlgorithm(configParams[0], configParams[1], configParams[2], configParams[3], configParams[4], configParams[5]);

    // Reset the file pointer to the beginning of file.
    configFile.clear();
    configFile.seekg(0, ios::beg);

    // Search from beginning of file.
    while(configFile.good())
    {
        getline(configFile, line);
        vector<double> tokens = parseStringDbl(line, ",\r\n\t");

        // If there aren't 3 parameters, skip that line.
        if(tokens.size() != 3)
            continue;

        // Do Firefly Algorithm.
        fireflyAlgorithm.runFireflyAlgorithm(tokens[0], tokens[1], tokens[2]);
    }

    //fireflyAlgorithm.printFAResults();

    fireflyAlgorithm.analyzeFAResults();
    fireflyAlgorithm.printFAAnalysis();

    //fireflyAlgorithm.saveFAResults();
    fireflyAlgorithm.saveFAAnalysis();
    //fireflyAlgorithm.saveEndingPopulation();

    // Close the configuration file.
    configFile.close();
}

// ******************************************************************************************************
// ******************************************************************************************************

void test_HS(const string &configFilename)
{
    // Open the config file.
    ifstream configFile;
    configFile.open(configFilename);
    if(configFile.fail())
    {
        cout << "Failed to open file: " << configFilename << endl;
        cout << "-----------------------------------------------\n";
        cout << "File is either not in the right directory\n";
        cout << "or does not exist.\n";
        cout << "-----------------------------------------------\n";
        cout << "Accepted File Formats: .txt" << endl;
        cout << "-----------------------------------------------\n";
        cout << "******** TERMINATING PROGRAM EXECUTION ********\n\n";
        exit(1);
    }

    // Get the configuration values from file.
    string line;
    vector<double> configParams;
    while(configFile.good())
    {
        getline(configFile, line);
        vector<double> tokens = parseStringDbl(line, ",\r\n\t");

        if(tokens.size() == 1)
            configParams.push_back(tokens[0]);
        else
            break;
    }

    // Create Harmony Object.
    HarmonySearch harmonySearch(configParams[0], configParams[1], configParams[2], configParams[3], configParams[4], configParams[5]);

    // Reset the file pointer to the beginning of file.
    configFile.clear();
    configFile.seekg(0, ios::beg);

    // Search from beginning of file.
    while(configFile.good())
    {
        getline(configFile, line);
        vector<double> tokens = parseStringDbl(line, ",\r\n\t");

        // If there aren't 3 parameters, skip that line.
        if(tokens.size() != 3)
            continue;

        // Do Harmony Search.
        harmonySearch.runHarmonySearch(tokens[0], tokens[1], tokens[2]);
    }

    //harmonySearch.printHSResults();

    harmonySearch.analyzeHSResults();
    harmonySearch.printHSAnalysis();

    //harmonySearch.saveHSResults();
    harmonySearch.saveHSAnalysis();

    // Close the configuration file.
    configFile.close();
}