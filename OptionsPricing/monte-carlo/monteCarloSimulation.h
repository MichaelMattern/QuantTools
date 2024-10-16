#ifndef MONTE_CARLO_SIMULATION_H
#define MONTE_CARLO_SIMULATION_H
#include "AmericanOptionPricing.h"
#include <vector>
#include <random>
#include <vector>
#include <utility> // For std::pair

class MonteCarloSimulation {
public:
    // Constructor
    MonteCarloSimulation(unsigned int numSimulations, unsigned int steps);

    // Method to calculate option price
    std::pair<double, double> calculateOptionPrice(AmericanOption& option);

    // Getter for number of simulations
    unsigned int getNumSimulations() const;

    // Setter for number of simulations
    void setNumSimulations(unsigned int numSimulations);

private:
    unsigned int numSimulations_; // Number of simulations
    unsigned int steps_;          // Number of time steps

    // Method to simulate a price path using Geometric Brownian Motion
    std::vector<double> simulatePricePath(double spotPrice, double riskFreeRate, double volatility, double expiry);

    // Random number generator for simulating paths
    std::mt19937 generator_;
};

#endif // MONTE_CARLO_SIMULATION_H
