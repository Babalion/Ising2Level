//
// Created by chris on 14.06.21.
//
#include "SpinLattice2level.h"

#include <cmath>


SpinLattice2level::SpinLattice2level(unsigned int sights)
        : J(1), performedSweeps(0), mt(rd()), u_int_dist(0, 1), u_float_dist(0, 1), sights(sights),
          spins(sights * sights), h(0) {
    initRandom();
}

SpinLattice2level::SpinLattice2level(unsigned int sights, short h)
        : J(1), performedSweeps(0), mt(rd()), u_int_dist(0, 1), u_float_dist(0, 1), sights(sights),
          spins(sights * sights), h(h) {
    initRandom();
}

SpinLattice2level::SpinLattice2level(const SpinLattice2level &sl)
        : J(sl.J), performedSweeps(sl.performedSweeps), mt(rd()), u_int_dist(0, 1), u_float_dist(0, 1),
          sights(sl.sights), h(sl.h) {
    spins = sl.spins;
}

void SpinLattice2level::printSpins() const {
    for (unsigned int i = 0; i < sights * sights; i++) {
        std::cout << spins[i];
        if ((i + 1) % sights == 0) {//right boarder
            std::cout << std::endl;
        } else {
            std::cout << "\t";
        }
    }
    std::cout << std::endl;
}

void SpinLattice2level::initRandom() {
    for (auto &s : spins) {
        s = u_int_dist(mt) == 0 ? -1 : 1;
    }
}

float SpinLattice2level::calcEnergy() const {
    int energyIt = 0;
    for (size_t i = 0; i < sights; ++i) {
        for (size_t j = 0; j < sights; ++j) {
            energyIt += calcEnergy(i, j);
        }
    }
    float energy= static_cast<float>(energyIt) / static_cast<float>(2*4 * sights * sights) + 0.5f;//scale to [0,1]
    return energy;
}

int SpinLattice2level::calcEnergy(unsigned int x, unsigned int y) const {
    return calcEnergy(x, y, spins[x + y * sights]);
}

// we choose cyclic boundary conditions
int SpinLattice2level::calcEnergy(unsigned int x, unsigned int y, int newSpin) const {
#ifdef DEBUG
    if(1!=std::abs(newSpin)){
        std::cerr<<"this is no valid spin ("<<newSpin<<")\n";
        exit(12);
    }
#endif
    const int J_val = J;
    const unsigned int i = x + y * sights;
    int energy = 0;

    // not at left boarder
    if (x > 0) {
        energy += (spins[i - 1] + h);
    } else {
        // give right boarder
        energy += (spins[i + (sights - 1)] + h);
    }

    // not at right boarder
    if (x < sights - 1) {
        energy += (spins[i + 1] + h);
    } else {
        //give right boarder
        energy += (spins[i - (sights - 1)] + h);
    }

    // not at top boarder
    if (y > 0) {
        energy += (spins[i - sights] + h);
    } else {
        //give bottom boarder
        energy += (spins[i + sights * (sights - 1)] + h);
    }

    // not at bottom boarder
    if (y < sights - 1) {
        energy += (spins[i + sights] + h);
    } else {
        //give top boarder
        energy += (spins[i - sights * (sights - 1)] + h);
    }

    energy=-1 * J_val * energy * newSpin;
    return energy;
}

float SpinLattice2level::calcMagnetization() const {
    float magnet = 0;
    for (const auto &s : spins) {
        magnet += s;
#ifdef DEBUG
        if (std::abs(magnet) > static_cast<float>(spins.size())) {
            std::cerr << "magnetization is " << magnet << " this is higher than possible\n";
            printSpins();
            exit(13);
        }
#endif
    }
#ifdef DEBUG
    std::cout << magnet << std::endl;
#endif
    return magnet / static_cast<float>(spins.size());
}

float SpinLattice2level::calcSusceptibility() {
    return 0;
}

int SpinLattice2level::calcHeatCapacity() {
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Markov Algorithms
////////////////////////////////////////////////////////////////////////////////


void metropolisSweep(SpinLattice2level &spinLattice, float temp) {
    for (size_t i = 0; i < spinLattice.getSights(); ++i) {
        for (size_t j = 0; j < spinLattice.getSights(); ++j) {
            const short newSpin = spinLattice.u_int_dist(spinLattice.mt) == 0 ? -1 : 1;
#ifdef DEBUG
            if (std::abs(newSpin))!=1){
                std::cerr << "new calculated Spin " << newSpin << " is not valid\n";
                exit(15);
            }
#endif
            if (newSpin == spinLattice(i, j)) {// spin has not changed, so skip all the work
            } else {
                const int newEnergy = spinLattice.calcEnergy(i, j, newSpin);
                const int deltaE = 2 * newEnergy; //if spin flips, energy changes from either -4 to 4 or -2 to 2

                if (deltaE < 0) {// energy decreases, so accept
                    spinLattice(i, j) *= -1;
                } else {
                    const float rand = spinLattice.u_float_dist(spinLattice.mt);
                    if (rand < std::exp(static_cast<float>(-1 * deltaE) / temp)) {
                        spinLattice(i, j) *= -1;
                    }
                }
            }
        }
    }
}

void metropolisSweep(SpinLattice2level &spinLattice, float temp, unsigned int iterations) {
    for (size_t i = 0; i < iterations; ++i) {
        metropolisSweep(spinLattice, temp);
    }
    spinLattice.performedSweeps++;
}

int heatBathSumOfNeighbours(SpinLattice2level &sl, unsigned int x, unsigned int y) {
    const auto sights = sl.getSights();
    const auto &spins = sl.getSpins();
    const unsigned int i = x + y * sights;
    int sum = 0;
    if (x > 0) {// not at left boarder
        sum += spins[i - 1];
    }
    if (x < sl.getSights() - 1) {// not at right boarder
        sum += spins[i + 1];
    }
    if (y > 0) {// not at top boarder
        sum += spins[i - sights];
    }
    if (y < sl.getSights() - 1) {// not at bottom boarder
        sum += spins[i + sights];
    }

    return sum;
}

void heatBathSweep(SpinLattice2level &spinLattice, float temp) {
    const auto J_val = spinLattice.J;
    auto mt = std::mt19937(spinLattice.rd());
    for (size_t i = 0; i < spinLattice.getSights(); i += 1) {
        for (size_t j = 0; j < spinLattice.getSights(); j += 1) {
            const int delta = heatBathSumOfNeighbours(spinLattice, i, j);
            const float k = static_cast<float>(-1 * J_val * delta) / temp;
            const float q = std::exp(-1.0f * k) / 2.0f / std::cosh(k);
            const float r = spinLattice.u_float_dist(mt);
            if (r < q) {
                spinLattice(i, j) = 1;
            } else {
                spinLattice(i, j) = -1;
            }
        }
    }
    spinLattice.performedSweeps++;
}

void heatBathSweepRandChoice(SpinLattice2level &spinLattice, float temp) {
    const auto J_val = spinLattice.J;
    static std::uniform_int_distribution<unsigned int> loc(0, spinLattice.getSights() - 1);

    for (size_t i = 0; i < spinLattice.getSights() * spinLattice.getSights(); i++) {
        const unsigned int x = loc(spinLattice.mt);
        const unsigned int y = loc(spinLattice.mt);

        const int delta = heatBathSumOfNeighbours(spinLattice, x, y);
        const float k = static_cast<float>(-1 * J_val * delta) / temp;
        const float q = std::exp(-1.0f * k) / 2.0f / std::cosh(k);
        const float r = spinLattice.u_float_dist(spinLattice.mt);
        if (r < q) {
            spinLattice(x, y) = 1;
        } else {
            spinLattice(x, y) = -1;
        }
    }
    spinLattice.performedSweeps++;
}
