//
// Created by chris on 14.06.21.
//
#include "SpinLattice2level.h"


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

void SpinLattice2level::initCold() {
    std::fill(spins.begin(), spins.end(), -1);
}

inline int SpinLattice2level::calcEnergy(const SpinLattice2level::Loc2d &loc, int newSpinVal) const {
#ifdef DEBUG
    if(1!=std::abs(newSpin)){
        std::cerr<<"this is no valid spin ("<<newSpin<<")\n";
        exit(12);
    }
#endif
    const auto neighbours = getNeighbours(loc);
    int energy = 0;

    for (const auto &n:neighbours) {
        energy += operator()(n);
    }
    energy = -1 * J * energy * newSpinVal;
    return energy;

}

float SpinLattice2level::calcEnergy() const {
    int energyIt = 0;
    for (size_t i = 0; i < sights; ++i) {
        for (size_t j = 0; j < sights; ++j) {
            energyIt += calcEnergy(SpinLattice2level::Loc2d(i, j));
        }
    }
    float energy = static_cast<float>(energyIt) / static_cast<float>(2 * 4 * sights * sights) + 0.5f;//scale to [0,1]
    return energy;
}


float SpinLattice2level::calcMagnetization() const {
    float magnet = 0;
    for (const auto &s : spins) {
        magnet += static_cast<float>(s);
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

////////////////////////////////////////////////////////////////////////////////
/// Markov Algorithms
////////////////////////////////////////////////////////////////////////////////


void metropolisSweep(SpinLattice2level &spinLattice, const float &temp) {
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
                const auto loc = SpinLattice2level::Loc2d(i, j);
                const int newEnergy = spinLattice.calcEnergy(loc, newSpin);
                const int deltaE = 2 * newEnergy; //if spin flips, energy changes from either -4 to 4 or -2 to 2

                if (deltaE < 0) {// energy decreases, so accept
                    spinLattice(loc) *= -1;
                } else {
                    const float rand = spinLattice.u_float_dist(spinLattice.mt);
                    if (rand < std::exp(static_cast<float>(-1 * deltaE) / temp)) {
                        spinLattice(loc) *= -1;
                    }
                }
            }
        }
    }
    spinLattice.performedSweeps++;
}

void metropolisSweep(SpinLattice2level &spinLattice, const float &temp, const unsigned int &iterations) {
    for (size_t i = 0; i < iterations; ++i) {
        metropolisSweep(spinLattice, temp);
    }
}

inline int heatBathSumOfNeighbours(const SpinLattice2level &sl, const SpinLattice2level::Loc2d &loc) {
    int sum = 0;
    auto neighbours = sl.getNeighbours(loc);
    for (const auto &n:neighbours) {
        sum += sl(n);
    }
    return sum;
}

void heatBathSweep(SpinLattice2level &spinLattice, const float &temp) {
    const auto J_val = spinLattice.J;
    auto mt = std::mt19937(spinLattice.rd());
    for (size_t i = 0; i < spinLattice.getSights(); i += 1) {
        for (size_t j = 0; j < spinLattice.getSights(); j += 1) {
            const auto loc = SpinLattice2level::Loc2d(i, j);

            const int delta = heatBathSumOfNeighbours(spinLattice, loc);
            const float k = static_cast<float>(-1 * J_val * delta) / temp;
            const float q = std::exp(-1.0f * k) / 2.0f / std::cosh(k);
            const float r = spinLattice.u_float_dist(mt);
            if (r < q) {
                spinLattice(loc) = 1;
            } else {
                spinLattice(loc) = -1;
            }
        }
    }
    spinLattice.performedSweeps++;
}

void heatBathSweep(SpinLattice2level &spinLattice, const float &temp, const unsigned int &iterations) {
    for (size_t i = 0; i < iterations; ++i) {
        heatBathSweep(spinLattice, temp);
    }
}

void heatBathSweepRandChoice(SpinLattice2level &spinLattice, const float &temp) {
    const auto J_val = spinLattice.J;
    static std::uniform_int_distribution<unsigned int> u(0, spinLattice.getSights() - 1);

    for (size_t i = 0; i < spinLattice.getSights() * spinLattice.getSights(); i++) {
        const unsigned int x = u(spinLattice.mt);
        const unsigned int y = u(spinLattice.mt);

        const auto loc = SpinLattice2level::Loc2d(x, y);

        const int delta = heatBathSumOfNeighbours(spinLattice, loc);
        const float k = static_cast<float>(-1 * J_val * delta) / temp;
        const float q = std::exp(-1.0f * k) / 2.0f / std::cosh(k);
        const float r = spinLattice.u_float_dist(spinLattice.mt);
        if (r < q) {
            spinLattice(loc) = 1;
        } else {
            spinLattice(loc) = -1;
        }
    }
    spinLattice.performedSweeps++;
}

void wolffSweep(SpinLattice2level &sl, const float &temp) {
    std::uniform_int_distribution<unsigned int> u(0, sl.getSights() - 1);

    // queue to save all locations of cluster
    // initialize with random location
    const SpinLattice2level::Loc2d startLoc{u(sl.mt), u(sl.mt)};
    sl(startLoc) *= -1;
    std::deque<SpinLattice2level::Loc2d> queue{startLoc};

    while (!queue.empty()) {
        // process first loc
        auto loc = queue.front();
        queue.pop_front();

        // calculate neighbours locations
        auto neighbours = sl.getNeighbours(loc);

        std::uniform_real_distribution<float> r(0, 1);
        for (auto n:neighbours) {
            if (sl(loc) * -1 == sl(n) &&
                1 - std::exp(-2.0f * static_cast<float>(sl.J) / temp) > r(sl.mt)) {
                sl(n) *= -1;
                queue.push_back(n);
            }
        }
    }
    sl.performedSweeps++;
}

void wolffSweep(SpinLattice2level &spinLattice, const float &temp, const unsigned int &iterations) {
    for (unsigned int i = 0; i < iterations; ++i) {
        wolffSweep(spinLattice, temp);
    }
}
