//
// Created by chris on 03.07.21.
//
#pragma once

#include "SpinLattice2level.h"

#include <chrono>
#include <ctime>
#include <fstream>
#include <string>
#include <thread>

class Simulation {
public:
    /**
     * Creates a set for a Ising-Simulation with given parameters
     * @param sights
     * @param numOfTemps
     * @param tempStart
     * @param tempEnd
     * @param numIterations
     * @param shuffleAgainAfter
     */
    Simulation(unsigned int sights, int numOfTemps, float tempStart, float tempEnd, int numIterations,
               unsigned int shuffleAgainAfter)
            : thermalizeSweeps(1000),sweepsPerIteration(1E4), sights(sights), tempStart(tempStart), tempEnd(tempEnd), numOfTemps(numOfTemps),
              numOfIterations(numIterations),
              shuffleAgainAfter(shuffleAgainAfter), sl(sights), isSimulated(false) {
        init();
    }

    void simulate_seq() {
        if (isSimulated) {
            std::cerr << "This simulation is already finished.\n";
        } else {
            for (unsigned int i = 0; i < temps.size(); i++) {
                // shuffle sl to obtain maybe a different equilibrate state
                if (i % shuffleAgainAfter == 0) {
                    std::cout<<"Shuffle!!<<\n";
                    sl.initRandom();
                    metropolisSweep(sl, temps[i], thermalizeSweeps);
                }
                if(i%(temps.size()/20)==0){
                    auto time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                    std::cout.precision(3);
                    std::cout << std::ctime(&time) << "N="<<sights<<"\tprogress:" << 100.0f*i/temps.size() <<"%"<< std::endl;
                }
                metropolisSweep(sl, temps[i], sweepsPerIteration);
                energies.push_back(sl.calcEnergy());
                magnetization.push_back(sl.calcMagnetization());
                susceptibility.push_back(sl.calcSusceptibility());
                heatCapacity.push_back(sl.calcHeatCapacity());
                if (i > 0 && i % numOfIterations == 0) {
                    auto time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                    std::cout << std::ctime(&time) << "done temp " << temps[i] << ". Remaining: "
                              << (temps.size() - i) / numOfIterations << std::endl;
                }

            }
            isSimulated = true;
        }
    }


    /**
 * sumlate the simulation parallelized
 */
    void simulate_par() {
        // -------------- Divide amount of work -------------------------------------------------
        static const unsigned int hardwareCon = std::thread::hardware_concurrency();
        static const unsigned int supportedThreads = hardwareCon == 0 ? 2 : hardwareCon;

        static const unsigned int workPerThread = numOfTemps / supportedThreads;
        static const unsigned int workRemaining = numOfTemps % supportedThreads;

        static unsigned int amountOfThreads;
        if (workPerThread == 0 && workRemaining > 0) //we have less work than threads
            amountOfThreads = numOfTemps;
        else //we have enough work --> use all cores
            amountOfThreads = supportedThreads;
        //TODO split up remaining work to all cores, not only the last one
        std::cout << amountOfThreads << " threads will be used for calculation." << std::endl;
        // --------------------------------------------------------------------------------------


        /// Create now data structure to divide and store work: initialize different ensembles
        std::vector<Simulation> Sims;
        Sims.reserve(amountOfThreads);
        float tempStep = (tempEnd - tempStart) / static_cast<float>(numOfTemps - 1);
        for (unsigned int i = 0; i < amountOfThreads; ++i) {
            float startTemp = tempStart + static_cast<float>(i * workPerThread) * tempStep;
            float endTemp;

            if (i == amountOfThreads - 1) {//this is the last case, include the rest now
                Sims.emplace_back(sights, workPerThread + workRemaining, startTemp, tempEnd, numOfIterations,
                                  shuffleAgainAfter);
            } else {
                endTemp = tempStart + static_cast<float>((i + 1) * workPerThread - 1) * tempStep;
                Sims.emplace_back(sights, workPerThread, startTemp, endTemp, numOfIterations, shuffleAgainAfter);
            }

        }
#ifdef DEBUG
        for (const auto &s:Sims) {
            std::cout.precision(3);
            std::cout << "Simulate from tS=" << s.getTempStart() << "\tuntil tE=" << s.getTempEnd() << "\twith\t"
                      << s.getNumOfTemps() << " temps, temps.size()=\t" << s.getTemps().size() << std::endl;
        }
#endif

        //start threads
        std::vector<std::thread> threads(amountOfThreads - 1);
        for (unsigned int i = 0; i < threads.size(); i++) {
            threads[i] = std::thread([&Sims, i]() { Sims[i].simulate_seq(); });
        }

        // do the remaining work
        Sims[amountOfThreads - 1].simulate_seq();

        for (auto &i : threads) {
            i.join();
        }

        for (const auto &Simulations:Sims) {
            //TODO check if measurements correspond to temps
#ifdef DEBUG
            std::cout << "Simulations.getTemps().size()=" << Simulations.getTemps().size() << std::endl;
#endif
            energies.insert(end(energies), begin(Simulations.getEnergies()), end(Simulations.getEnergies()));
            magnetization.insert(end(magnetization), begin(Simulations.getMagnetization()),
                                 end(Simulations.getMagnetization()));
            susceptibility.insert(end(susceptibility), begin(Simulations.getSusceptibility()),
                                  end(Simulations.getSusceptibility()));
            heatCapacity.insert(end(heatCapacity), begin(Simulations.getHeatCapacity()),
                                end(Simulations.getHeatCapacity()));
        }
        isSimulated = true;
#ifdef DEBUG
        if(temps.size()!=numOfTemps*numOfIterations){
            std::cerr<<"results size not matched!\n";
            exit(15);
        }
#endif

    }


    [[nodiscard]] unsigned int getSights() const {
        return sights;
    }

    [[nodiscard]] unsigned int getNumOfTemps() const {
        return numOfTemps;
    }

    [[nodiscard]] float getTempStart() const {
        return tempStart;
    }

    [[nodiscard]] float getTempEnd() const {
        return tempEnd;
    }

    [[nodiscard]] unsigned int getNumOfIterations() const {
        return numOfIterations;
    }

    [[nodiscard]] unsigned int getShuffleAgainAfter() const {
        return shuffleAgainAfter;
    }

    [[nodiscard]] const std::vector<float> &getTemps() const {
        return temps;
    }

    [[nodiscard]] const std::vector<float> &getEnergies() const {
        return energies;
    }

    [[nodiscard]] const std::vector<float> &getMagnetization() const {
        return magnetization;
    }

    [[nodiscard]] const std::vector<float> &getSusceptibility() const {
        return susceptibility;
    }

    [[nodiscard]] const std::vector<float> &getHeatCapacity() const {
        return heatCapacity;
    }

public:
    unsigned int thermalizeSweeps;
    unsigned int sweepsPerIteration;
private:
    /// Parameters for simulation
    unsigned int sights;
    float tempStart;
    float tempEnd;
    unsigned int numOfTemps;
    unsigned int numOfIterations;
    unsigned int shuffleAgainAfter;

    /// Results of simulation
    std::vector<float> temps;
    std::vector<float> energies;
    std::vector<float> magnetization;
    std::vector<float> susceptibility;
    std::vector<float> heatCapacity;

    // SpinLattice for simulation
    SpinLattice2level sl;
    bool isSimulated;
private:
    void init() {
        // reserve memory for results
        temps.reserve(numOfTemps * numOfIterations);
        energies.reserve(numOfTemps * numOfIterations);
        magnetization.reserve(numOfTemps * numOfIterations);
        susceptibility.reserve(numOfTemps * numOfIterations);
        heatCapacity.reserve(numOfTemps * numOfIterations);

        // calculate temps
        for (unsigned int i = 0; i < numOfTemps; ++i) {
            float temp = tempStart + static_cast<float>(i) * (tempEnd - tempStart) / static_cast<float>(numOfTemps - 1);
            for (unsigned int j = 0; j < numOfIterations; ++j) {
                temps.push_back(temp);
            }
        }

    }
};
