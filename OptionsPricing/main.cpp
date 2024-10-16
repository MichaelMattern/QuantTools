#include "option-models/AmericanOptionPricing.h"
#include "monte-carlo/MonteCarloSimulation.h"
#include <numeric> 
#include <iostream>

int main() {
    std::cout << "Strategy One: binomial tree" << "\n" << std::endl;

// American Option Pricing
    double strike = 581.00;                                         // Strike price
    double expiry = 0.032;                                         // Time to maturity (in years) | 1 week => 1/52 | 2 weeks => 2/52 | etc...
    double spotPrice = 580.53;                                      // Current stock price
    double riskFreeRate = 0.046;                                    // Risk-free interest rate (%)
    double volatility = 0.1160;                                     // Volatility (%)
    unsigned int steps = 1000;                                      // Number of steps in the binomial tree (American Option)
    bool isCallOption = true;                                       // Call option

    // Create an AmericanOption instance
    AmericanOption option(strike, expiry, spotPrice, riskFreeRate, volatility, steps, isCallOption);

    std::cout << "Strike Price:\t" << strike << std::endl;
    std::cout << "Spot Price:\t" << spotPrice << std::endl;

    // Calculate the option price using the binomial tree method
    double optionPrice = option.calculateOptionPrice();

    // Print Result
    std::cout << "Option Price:\t" << optionPrice << std::endl;

    // Simulate final stock price and calculate payoff
    double finalPrice = 585.0;  // Final stock price at option expiry
    double payoff = option.calculatePayoff(finalPrice);

    std::cout << "Option Payoff:\t" << payoff << "\n" << std::endl;


    std::cout << "----------------------------------------------" << "\n" << std::endl;


// Monte Carlo Simulation for a call option
    std::cout << "Strategy Two: Monte Carlo Simulation" << "\n" << std::endl;

    std::vector<std::pair<double, double>> regressionArray;          // Stores the expected payoff

//   double strike = 581.00;                                         // Strike price
//   double expiry = (2.0 / 52.0);                                   // Time to maturity (in years) | 1 week => 1/52 | 2 weeks => 2/52 | etc...
//   double spotPrice = 580.32;                                      // Current stock price
//   double riskFreeRate = 0.046;                                    // Risk-free interest rate (%)
//   double volatility = 0.1157;                                     // Volatility (%)
//   unsigned int steps = 1;                                         // Number of step the random path takes from start
//   bool isCallOption = true;                                       // Call option

    steps = 1;                                         // Number of step the random path takes from start


    // Create an AmericanOption object
//    AmericanOption option(strike, expiry, spotPrice, riskFreeRate, volatility, steps, isCallOption);



// Create a MonteCarloSimulation object
MonteCarloSimulation mcSim(10000, steps);

    std::vector<double> expectedPayoffs; // Stores the expected payoff
    std::vector<double> inTheMoneyProbabilities; // Stores the probability of ending ITM
    for (int i = 0; i < 390; i++) {
        std::cout << "Simulation " << i + 1 << std::endl;

        // Calculate option price using Monte Carlo simulation
        auto result = mcSim.calculateOptionPrice(option); // This returns a pair
        expectedPayoffs.push_back(result.first); // Store expected payoff
        inTheMoneyProbabilities.push_back(result.second); // Store ITM probability
        std::cout << "\n";

        // Update the number of simulations for the next iteration if needed
        unsigned int newSteps = steps * 7; // Update the number of steps as needed
        mcSim.setNumSimulations(newSteps);
        std::cout << "\n";

        //Update # of steps
        newSteps = steps * 7;
        mcSim.setNumSimulations(newSteps);
    }

    // Calculate the mean of expected payoffs
    double expectedPayoffSum = std::accumulate(expectedPayoffs.begin(), expectedPayoffs.end(), 0.0);
    double expectedPayoffMean = expectedPayoffSum / expectedPayoffs.size(); // Use expectedPayoffs.size() for accuracy

    // Calculate the mean of in-the-money probabilities
    double inTheMoneyProbabilitySum = std::accumulate(inTheMoneyProbabilities.begin(), inTheMoneyProbabilities.end(), 0.0);
    double inTheMoneyProbabilityMean = inTheMoneyProbabilitySum / inTheMoneyProbabilities.size(); // Use inTheMoneyProbabilities.size() for accuracy

    std::cout << "\nMean of the expected payoff: " << expectedPayoffMean << std::endl;
    std::cout << "Mean of the probability of ending ITM: " << inTheMoneyProbabilityMean * 100 << "%" << std::endl; // Convert to percentage

    return 0;



}