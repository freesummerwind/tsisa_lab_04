#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <vector>

double fitnessFunction(const double x, const double y) {
    return sin(x) * exp(-1.*pow(y, 2)) / (1. + pow(x, 2) + pow(y, 2));
}

double randomInRange(const double lower, const double upper) {
    return lower + rand() * 1./RAND_MAX * (upper - lower);
}

struct Point {
    double _x;
    double _y;
};

bool operator <(const Point& lhs, const Point& rhs) {
    return fitnessFunction(lhs._x, lhs._y) < fitnessFunction(rhs._x, rhs._y);
}

Point randomPoint(const double xLower, const double xUpper,
                  const double yLower, const double yUpper) {
    return {randomInRange(xLower, xUpper), randomInRange(yLower, yUpper)};
}

using Population = std::vector<Point>;

void printLine() {
    std::cout << '+' << std::string(5, '-') << '+' << std::string(12, '-')
              << '+' << std::string(12, '-') << '+' << std::string(12, '-')
              << '+' << std::string(12, '-') << '+' << std::string(12, '-')
              << '+' << std::endl;
}

void printTableHead() {
    printLine();
    std::cout << '|' << "  N  " << '|' << "     X      "
              << '|' << "     Y      " << '|' << "    FIT     "
              << '|' << "    MAX     " << '|' << "  AVERAGE   "
              << '|' << std::endl;
    printLine();
}

void printPopulation(const int generation, const Population& population) {
    for (size_t i = 0; i < population.size(); ++i) {
        if (i == 0) {
            std::cout << '|' << std::setw(4) << generation << " | ";
        } else {
            std::cout << '|' << std::setw(7) << " | ";
        }
        std::cout << std::setw(10) << population[i]._x << " | "
                  << std::setw(10) << population[i]._y << " | "
                  << std::setw(10) << fitnessFunction(population[i]._x, population[i]._y) << " | ";
        if (i == 0) {
            double max = fitnessFunction(population[0]._x, population[0]._y), average = 0;
            for (auto i : population) {
                if (fitnessFunction(i._x, i._y) > max) max = fitnessFunction(i._x, i._y);
                average += fitnessFunction(i._x, i._y);
            }
            average /= population.size();
            std::cout << std::setw(10) << max << " | "
                      << std::setw(10) << average << " |\n";
        } else {
            std::cout << std::setw(12) << '|' << std::setw(14) << "|\n";
        }
    }
    printLine();
}

void geneticAlgorithm(const double xLower, const double xUpper,
                      const double yLower, const double yUpper,
                      const int pointsNumber, const int generationsNumber,
                      const double mutationProbability) {
    Population currentPopulation;
    for (size_t i = 0; i < pointsNumber; ++i) {
        currentPopulation.push_back(randomPoint(xLower, xUpper, yLower, yUpper));
    }
    std::sort(currentPopulation.rbegin(), currentPopulation.rend());
    const double epsilon = .5;

    printTableHead();
    for(int i = 1; i <= generationsNumber; ++i) {
        printPopulation(i - 1, currentPopulation);
        //generation of new population
        Population newPopulation;
        newPopulation.push_back(Point{currentPopulation[0]._x, currentPopulation[1]._y});
        newPopulation.push_back(Point{currentPopulation[1]._x, currentPopulation[0]._y});
        newPopulation.push_back(Point{currentPopulation[0]._x, currentPopulation[2]._y});
        newPopulation.push_back(Point{currentPopulation[2]._x, currentPopulation[0]._y});
        newPopulation.push_back(Point{currentPopulation[0]._x, currentPopulation[3]._y});
        newPopulation.push_back(Point{currentPopulation[3]._x, currentPopulation[0]._y});
        //mutation
        for (size_t j = 0; j < pointsNumber; ++j) {
            double probability = randomInRange(0., 1.);
            if (probability < mutationProbability) {
                double deltaX = randomInRange(-epsilon/2, epsilon/2);
                double deltaY = randomInRange(-epsilon/2, epsilon/2);
                newPopulation[j]._x += deltaX;
                if (newPopulation[j]._x < xLower) newPopulation[j]._x = xLower;
                if (newPopulation[j]._x > xUpper) newPopulation[j]._x = xUpper;
                newPopulation[j]._y += deltaY;
                if (newPopulation[j]._y < yLower) newPopulation[j]._y = yLower;
                if (newPopulation[j]._y > yUpper) newPopulation[j]._y = yUpper;
            }
        }
        std::sort(newPopulation.rbegin(), newPopulation.rend());
        currentPopulation = newPopulation;
    }
    printPopulation(generationsNumber, currentPopulation);
}

const double X_LOWER = 0.;
const double X_UPPER = 2.;
const double Y_LOWER = -2.;
const double Y_UPPER = 2.;
const int POINTS_NUMBER = 6;
const double MUTATION_PROBABILITY = .25;

int main() {
    std::cout << "Variant 15.\nFunction: sin(x) * exp(-y^2) / (1 + x^2 + y^2), D: ("
         << X_LOWER << ", " << X_UPPER << ") x (" << Y_LOWER << ", " << Y_UPPER << ")\n";
    srand(time(nullptr));
    std::cout << "For 10 generations:\n";
    geneticAlgorithm(X_LOWER, X_UPPER, Y_LOWER, Y_UPPER, POINTS_NUMBER,
            10, MUTATION_PROBABILITY);

    std::cout << "For 100 generations:\n";
    geneticAlgorithm(X_LOWER, X_UPPER, Y_LOWER, Y_UPPER, POINTS_NUMBER,
                     100, MUTATION_PROBABILITY);

    return 0;
}
