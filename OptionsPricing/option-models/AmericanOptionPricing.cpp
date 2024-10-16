#include "AmericanOptionPricing.h"
#include <cmath>
#include <algorithm>

// Constructor
AmericanOption::AmericanOption(double strike, double expiry, double spotPrice, double riskFreeRate,
    double volatility, unsigned int steps, bool isCallOption)
    : strike_(strike), expiry_(expiry), spotPrice_(spotPrice), riskFreeRate_(riskFreeRate),
    volatility_(volatility), steps_(steps), isCallOption_(isCallOption) {}

// Destructor
AmericanOption::~AmericanOption() {}

// Calculate the up factor in the binomial tree
double AmericanOption::calculateUpFactor() const {
    return std::exp(volatility_ * std::sqrt(expiry_ / steps_));
}

// Calculate the down factor in the binomial tree
double AmericanOption::calculateDownFactor() const {
    return 1.0 / calculateUpFactor();
}

// Calculate the risk-neutral probability
double AmericanOption::calculateRiskNeutralProbability() const {
    double up = calculateUpFactor();
    double down = calculateDownFactor();
    return (std::exp(riskFreeRate_ * (expiry_ / steps_)) - down) / (up - down);
}

// Generate the stock price tree
std::vector<std::vector<double>> AmericanOption::generateStockPriceTree() const {
    std::vector<std::vector<double>> stockTree(steps_ + 1, std::vector<double>(steps_ + 1));

    double up = calculateUpFactor();
    double down = calculateDownFactor();

    for (unsigned int i = 0; i <= steps_; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            stockTree[i][j] = spotPrice_ * std::pow(up, static_cast<int>(j)) * std::pow(down, static_cast<int>(i - j));
        }
    }

    return stockTree;
}

// Generate the option price tree
std::vector<std::vector<double>> AmericanOption::generateOptionPriceTree(const std::vector<std::vector<double>>& stockTree) const {
    std::vector<std::vector<double>> optionTree(steps_ + 1, std::vector<double>(steps_ + 1));

    double p = calculateRiskNeutralProbability();
    double discountFactor = std::exp(-riskFreeRate_ * (expiry_ / steps_));

    // Calculate option value at the final nodes (expiry time)
    for (unsigned int j = 0; j <= steps_; ++j) {
        double payoff = 0.0;
        if (isCallOption_) {
            payoff = std::max(0.0, stockTree[steps_][j] - strike_);
        }
        else {
            payoff = std::max(0.0, strike_ - stockTree[steps_][j]);
        }
        optionTree[steps_][j] = payoff;
    }

    // Backward induction to calculate the option price at earlier nodes
    for (int i = steps_ - 1; i >= 0; --i) {
        for (unsigned int j = 0; j <= i; ++j) {
            double earlyExerciseValue = 0.0;
            if (isCallOption_) {
                earlyExerciseValue = std::max(0.0, stockTree[i][j] - strike_);
            }
            else {
                earlyExerciseValue = std::max(0.0, strike_ - stockTree[i][j]);
            }

            double holdValue = discountFactor * (p * optionTree[i + 1][j + 1] + (1.0 - p) * optionTree[i + 1][j]);
            optionTree[i][j] = std::max(earlyExerciseValue, holdValue);
        }
    }

    return optionTree;
}

// Main function to calculate option price
double AmericanOption::calculateOptionPrice() {
    std::vector<std::vector<double>> stockTree = generateStockPriceTree();
    std::vector<std::vector<double>> optionTree = generateOptionPriceTree(stockTree);

    // Return the option price at the root node
    return optionTree[0][0];
}
