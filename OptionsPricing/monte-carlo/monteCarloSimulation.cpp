#include "MonteCarloSimulation.h"
#include <cmath>
#include <algorithm>
#include <iostream>


// Constructor
MonteCarloSimulation::MonteCarloSimulation(unsigned int numSimulations, unsigned int steps)
    : numSimulations_(numSimulations), steps_(steps), generator_(std::random_device()()) {}

// Simulate a single price path using Geometric Brownian Motion (GBM)
std::vector<double> MonteCarloSimulation::simulatePricePath(double spotPrice, double riskFreeRate, double volatility, double expiry) {
    std::vector<double> prices(steps_ + 1);
    prices[0] = spotPrice;
    double dt = expiry / steps_; // Time increment for each step
    double drift = (riskFreeRate - 0.5 * std::pow(volatility, 2)) * dt;
    double diffusion = volatility * std::sqrt(dt);

    std::normal_distribution<double> normal(0.0, 1.0); // Standard normal distribution

    // Generate the price path
    for (unsigned int i = 1; i <= steps_; ++i) {
        double randomShock = normal(generator_);
        prices[i] = prices[i - 1] * std::exp(drift + diffusion * randomShock);
    }

    return prices;
}

// Monte Carlo simulation to calculate the option price
std::pair<double, double> MonteCarloSimulation::calculateOptionPrice(AmericanOption& option) {
    double totalPayoff = 0.0;
    unsigned int inTheMoneyCount = 0;

    for (unsigned int i = 0; i < numSimulations_; ++i) {
        // Simulate a price path
        std::vector<double> pricePath = simulatePricePath(option.getSpotPrice(), option.getRiskFreeRate(), option.getVolatility(), option.getExpiry());

        // Get the final price (at maturity)
        double finalPrice = pricePath.back();

        // Calculate payoff depending on whether it's a call or put option
        double payoff = 0.0;
        if (option.isCallOption()) {
            payoff = std::max(0.0, finalPrice - option.getStrike()); // Call option payoff
        }
        else {
            payoff = std::max(0.0, option.getStrike() - finalPrice);  // Put option payoff
        }

        // Count the option if it finishes in the money
        if (payoff > 0) {
            inTheMoneyCount++;
        }

        // Accumulate total payoff
        totalPayoff += payoff;
    }

    // Discount the total payoff to present value using the risk-free rate
    double discountFactor = std::exp(-option.getRiskFreeRate() * option.getExpiry());
    double expectedPayoff = (totalPayoff / numSimulations_) * discountFactor;

    // Calculate the probability of the option being in-the-money
    double inTheMoneyProbability = static_cast<double>(inTheMoneyCount) / numSimulations_;

    // Output the results
    std::cout << "Expected Option Payoff:\t" << expectedPayoff << std::endl;
    std::cout << "ITM Probability:\t" << inTheMoneyProbability * 100 << "%" << std::endl;

    return std::make_pair(expectedPayoff, inTheMoneyProbability);
}

// Getter for the number of simulations
unsigned int MonteCarloSimulation::getNumSimulations() const {
    return numSimulations_;
}

// Setter for the number of simulations
void MonteCarloSimulation::setNumSimulations(unsigned int numSimulations) {
    numSimulations_ = numSimulations;
}
