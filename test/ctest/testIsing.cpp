//
// Created by chris on 06.07.21.
//
#include <algorithm>
#include "assert_macro.h"
#include <chrono>

#include "../../Simulation.h"

int test_SpinLattice2level() {
    int err_code = 0;

    // check initialization...
    SpinLattice2level sl0(3);
    sl0.initRandom();
    assertEqual (sl0.getSights() == 3);
    for (auto &s:sl0.getSpins()) {
        assertEqual (std::abs(s) == 1);
    }
    metropolisSweep(sl0, 2.8, 1);
    for (auto &s:sl0.getSpins()) {
        assertEqual (std::abs(s) == 1);
    }
    SpinLattice2level sl1(137);
    assertEqual (sl1.getSights() == 137);
    sl1.initRandom();
    for (auto &s:sl1.getSpins()) {
        assertEqual (std::abs(s) == 1);
    }
    metropolisSweep(sl1, 1.4, 30);
    for (auto &s:sl1.getSpins()) {
        assertEqual (std::abs(s) == 1);
    }
    for (int i = 0; i < 1E3; ++i) {
        metropolisSweep(sl1, 2.4);
        assertEqual (std::abs(sl1.calcMagnetization()) <= 1);
        assertEqual (sl1.calcMagnetization() >= -1);
        assertEqual (sl1.calcEnergy() <= 1);
        assertEqual (sl1.calcEnergy() >= 0);
    }

    return err_code;
}


int test_Simulation_seq() {
    std::cout << std::endl << "Testing sequential mode" << std::endl << std::endl;
    int err_code = 0;

    // check initialization...
    const int numOfTemps = 16;
    const int numOfIterations = 1000;
    Simulation Sim(128, numOfTemps, 1, 8, numOfIterations, UINT32_MAX);
    Sim.sweepsPerIteration = 2;
    Sim.thermalizeSweeps = 50;
    assertEqual(Sim.getSights() == 128);
    assertEqual(Sim.getNumOfTemps() == numOfTemps);
    assertEqual(Sim.getTemps().size() == numOfTemps * numOfIterations);
    Sim.simulate_seq();
    assertEqual(Sim.getTemps().back() == 8);
    std::cout << "Sim.getTemps().back()=" << Sim.getTemps().back() << std::endl;
    assertEqual(Sim.getTemps()[0] == 1);
    assertEqual(Sim.getTemps().size() == numOfIterations * numOfTemps);
    assertEqual(std::is_sorted(Sim.getTemps().begin(), Sim.getTemps().end()) == true);

    assertEqual(Sim.getEnergies().size() == numOfIterations * numOfTemps);
    assertEqual(Sim.getMagnetization().size() == numOfIterations * numOfTemps);
    std::cout << "Sim.getTemps().size()= " << Sim.getTemps().size() << std::endl;
    for (int i = 0; i < numOfIterations * numOfTemps; ++i) {
        std::cout << Sim.getTemps()[i] << std::endl;
        assertEqual(Sim.getEnergies()[i] >= 0);
        assertEqual(Sim.getEnergies()[i] <= 1);
        assertEqual(Sim.getMagnetization()[i] <= 1);
        assertEqual(Sim.getMagnetization()[i] >= -1);
    }
    //Sim.simulate_seq();

    return err_code;
}

int test_Simulation_par() {
    std::cout << std::endl << "Testing parallel mode" << std::endl << std::endl;
    int err_code = 0;

    // check initialization...
    const float startTemp = 2;
    const int numOfTemps = 16;
    const int numOfIterations = 50;
    Simulation Sim(16, numOfTemps, startTemp, 8, numOfIterations, 10);
    assertEqual(Sim.getSights() == 16);
    assertEqual(Sim.getNumOfTemps() == numOfTemps);
    assertEqual(Sim.getTemps().size() == numOfTemps * numOfIterations);

    Sim.simulate_par();

    assertEqual(Sim.getTemps().back() == 8);
    std::cout << "Sim.getTemps().back()=" << Sim.getTemps().back() << std::endl;
    assertEqual(Sim.getTemps()[0] == startTemp);
    assertEqual(Sim.getTemps().size() == numOfIterations * numOfTemps);
    assertEqual(std::is_sorted(Sim.getTemps().begin(), Sim.getTemps().end()) == true);

    assertEqual(Sim.getEnergies().size() == numOfIterations * numOfTemps);
    assertEqual(Sim.getMagnetization().size() == numOfIterations * numOfTemps);
    std::cout << "Sim.getTemps().size()= " << Sim.getTemps().size() << std::endl;
    for (unsigned int i = 0; i < Sim.getTemps().size(); ++i) {
        //std::cout << Sim.getTemps()[i] << std::endl;
        assertEqual(Sim.getEnergies()[i] >= 0);
        assertEqual(Sim.getEnergies()[i] <= 1);
        assertEqual(Sim.getMagnetization()[i] <= 1);
        assertEqual(Sim.getMagnetization()[i] >= -1);
    }

    return err_code;
}

int main() {
    int err_code = 0;
    auto begin = std::chrono::steady_clock::now();
    assertEqual (test_SpinLattice2level() == 0);
    assertEqual (test_Simulation_seq() == 0);
    assertEqual (test_Simulation_par() == 0);

    auto end = std::chrono::steady_clock::now();
    std::cout << "Time needed = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]"
              << std::endl;
    return err_code;
}
