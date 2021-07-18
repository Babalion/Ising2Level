//
// Created by chris on 14.06.21.
//
#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>

// a quadratic 2-level ising-lattice
class SpinLattice2level {
public:
    /**
     * initialize an ising-lattice with sights² spins and no external magnetic field
     * @param sights length of the quadratic lattice
     */
    explicit SpinLattice2level(unsigned int sights);

    /**
     * initialize an ising-lattice with sights² spins and extern magnetic field h
     * @param sights length of the quadratic lattice
     * @param h strength of extern magnetic field
     */
    SpinLattice2level(unsigned int sights, short h);

    /**
     * Copy-constructor: Doesn't initialize random.
     * @param sl
     */
    SpinLattice2level(const SpinLattice2level &sl);

    ~SpinLattice2level() = default;


    // prints a matrix-scheme to the console
    void printSpins() const;

    // Reinitialize all spins with random values
    void initRandom();


    short operator()(unsigned int x, unsigned int y) const {
        assert(x < sights && y < sights);
        return spins[x + y * sights];
    }

    short &operator()(unsigned int x, unsigned int y) {
        assert(x < sights && y < sights);
        return spins[x + y * sights];
    }

    /**
     * calculates normalized energy of system: sum over all spins, divided by 4N²
     * @return energy between 0 and 1
    */
    [[nodiscard]] float calcEnergy() const;

    [[nodiscard]] int calcEnergy(unsigned int x, unsigned int y) const;

    [[nodiscard]] int calcEnergy(unsigned int x, unsigned int y, int newSpinVal) const;

    /**
     *calculates normalized magnetization: sum over all spins, divided by N²
     * @return magnetization between -1 and 1
     */
    [[nodiscard]] float calcMagnetization() const;

    [[nodiscard]] static float calcSusceptibility();

    [[nodiscard]] static int calcHeatCapacity();


    [[nodiscard]] unsigned int getSights() const {
        return sights;
    }

    [[nodiscard]] const std::vector<short> &getSpins() const {
        return spins;
    }

    int J;

    unsigned int performedSweeps;

    std::random_device rd;
    std::mt19937 mt;
    /**
     * returns int either 0 or 1
     */
    std::uniform_int_distribution<short> u_int_dist;
    /**
     * returns a float between 0 and 1
     */
    std::uniform_real_distribution<float> u_float_dist;
private:
    unsigned int sights;
    std::vector<short> spins;
    int h;
};

void metropolisSweep(SpinLattice2level &spinLattice, const float &temp);

void metropolisSweep(SpinLattice2level &spinLattice, const float &temp, const unsigned int &iterations);

void heatBathSweep(SpinLattice2level &spinLattice, const float &temp);
void heatBathSweep(SpinLattice2level &spinLattice, const float &temp,const unsigned int &iterations);

void heatBathSweepRandChoice(SpinLattice2level &spinLattice, const float &temp);

void wolffSweep(SpinLattice2level &sl, const float &temp);

void wolffSweep(SpinLattice2level &spinLattice, const float &temp, const unsigned int &iterations);
