#include <algorithm>
#include <vector>
#ifndef AMERICAN_OPTION_PRICING_H
#define AMERICAN_OPTION_PRICING_H
class AmericanOption {
public:
    // Constructor
    AmericanOption(double strike, double expiry, double spotPrice, double riskFreeRate,
        double volatility, unsigned int steps, bool isCallOption);

    // Destructor
    ~AmericanOption();

    // Method to calculate option price using the binomial tree
    double calculateOptionPrice();

    // Method to calculate the payoff based on final price (call or put)
    double calculatePayoff(double finalPrice) const {
        if (this->isCallOption_) {
            return std::max(0.0, finalPrice - this->strike_);  // Call option payoff
        }
        else {
            return std::max(0.0, this->strike_ - finalPrice);  // Put option payoff
        }
    }

    // Getters
    double getStrike() const {
        return strike_;
    }

    double getExpiry() const {
        return expiry_;
    }

    double getSpotPrice() const {
        return spotPrice_;
    }

    double getRiskFreeRate() const {
        return riskFreeRate_;
    }

    double getVolatility() const {
        return volatility_;
    }

    unsigned int getSteps() const {
        return steps_;
    }

    bool isCallOption() const {
        return isCallOption_;
    }

    // Setters
    void setStrike(double strike) {
        strike_ = strike;
    }

    void setExpiry(double expiry) {
        expiry_ = expiry;
    }

    void setSpotPrice(double spotPrice) {
        spotPrice_ = spotPrice;
    }

    void setRiskFreeRate(double riskFreeRate) {
        riskFreeRate_ = riskFreeRate;
    }

    void setVolatility(double volatility) {
        volatility_ = volatility;
    }

    void setSteps(unsigned int steps) {
        steps_ = steps;
    }

    void setCallOption(bool isCallOption) {
        isCallOption_ = isCallOption;
    }

private:
    // Member variables
    double strike_;         // Strike price
    double expiry_;         // Time to maturity (in years)
    double spotPrice_;      // Current stock price
    double riskFreeRate_;   // Risk-free interest rate
    double volatility_;     // Volatility of the underlying asset
    unsigned int steps_;    // Number of steps in the binomial tree
    bool isCallOption_;     // Call option if true, otherwise put option

    // Helper methods for binomial tree construction
    double calculateUpFactor() const;
    double calculateDownFactor() const;
    double calculateRiskNeutralProbability() const;

    // Function to create the binomial tree and calculate option price
    std::vector<std::vector<double>> generateStockPriceTree() const;
    std::vector<std::vector<double>> generateOptionPriceTree(const std::vector<std::vector<double>>& stockTree) const;

};

#endif // AMERICAN_OPTION_PRICING_H
