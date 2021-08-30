//
// Created by chris on 03.07.21.
//
#pragma once

#include "SpinLattice2level.h"

#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <string>
#include <thread>
#include <iomanip>

class Simulation {
public:
    /**
     * Creates a set for a Ising-Simulation with given parameters
     * @param sights
     * @param numOfTemps
     * @param tempStart
     * @param tempEnd
     * @param numIterations
     * @param shuffleAgainAfter this leads to reinitialize the spins after given number. Set to UINT32_MAX if you would like to use always the same ensemble
     */
    Simulation(unsigned int sights, unsigned int numOfTemps, float tempStart, float tempEnd, unsigned int numIterations,
               unsigned int shuffleAgainAfter)
            : thermalizeSweeps(10), sweepsPerIteration(1), sights(sights), tempStart(tempStart), tempEnd(tempEnd),
              numOfTemps(numOfTemps), numOfIterations(numIterations), shuffleAgainAfter(shuffleAgainAfter),
              tempIndexATM(0), amountOfThreads(1), amountOfWorkingThreads(0), printStat(true), sl(sights),
              isSimulated(false) {
        // reserve memory for results
        temps.reserve(numOfTemps * numOfIterations);
        energies.reserve(numOfTemps * numOfIterations);
        magnetization.reserve(numOfTemps * numOfIterations);

        // calculate temps
        for (unsigned int i = 0; i < numOfTemps; ++i) {
            float temp = tempStart + static_cast<float>(i) * (tempEnd - tempStart) / static_cast<float>(numOfTemps - 1);
            for (unsigned int j = 0; j < numOfIterations; ++j) {
                temps.push_back(temp);
            }
        }
    }

    void simulate_seq() {
        if (isSimulated) {
            std::cerr << "This simulation is already finished.\n";
        } else {
            amountOfWorkingThreads = 1;
            for (unsigned int i = 0; i < temps.size(); i++) {
                // shuffle sl to obtain maybe a different equilibrate state
                if (i % shuffleAgainAfter == 0) {
                    sl.initRandom();
                    wolffSweep(sl, temps[i], thermalizeSweeps);
                }
                if (printStat && i % 10000 == 0) {
                    printStatus();
                }

                wolffSweep(sl, temps[i], sweepsPerIteration);
                energies.push_back(sl.calcEnergy());
                magnetization.push_back(sl.calcMagnetization());
                tempIndexATM++;

            }
            isSimulated = true;
            amountOfWorkingThreads = 0;
        }
    }


    /**
 * simulate the simulation parallelized
 */
    void simulate_par() {
        // -------------- Divide amount of work -------------------------------------------------
        static const unsigned int hardwareCon = std::thread::hardware_concurrency();
        static const unsigned int supportedThreads = hardwareCon == 0 ? 2 : hardwareCon;

        const unsigned int workPerThread = numOfTemps / supportedThreads;
        const unsigned int workRemaining = numOfTemps % supportedThreads;

        if (workPerThread == 0 && workRemaining > 0) //we have less work than threads
            amountOfThreads = numOfTemps;
        else //we have enough work --> use all cores
            amountOfThreads = supportedThreads;
        //TODO split up remaining work to all cores, not only the last one
        std::cout << amountOfThreads << " threads will be used for calculation." << std::endl;
        // --------------------------------------------------------------------------------------


        /// Create now data structure to divide and store work: initialize different ensembles
        // TODO shuffle temperatures bec. different temps may take different simulation time (see Wolff for low vs high T)
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
            // deactivate std::cout of those sims
            Sims.back().printStat = false;

        }
#ifdef DEBUG
        for (const auto &s:Sims) {
            std::cout.precision(3);
            std::cout << "Simulate from tS=" << s.getTempStart() << "\tuntil tE=" << s.getTempEnd() << "\twith\t"
                      << s.getNumOfTemps() << " temps, temps.size()=\t" << s.getTemps().size() << std::endl;
        }
#endif

        //start threads
        std::vector<std::thread> threads(amountOfThreads);
        for (unsigned int i = 0; i < threads.size(); i++) {
            threads[i] = std::thread([&Sims, i]() { Sims[i].simulate_seq(); });
        }

        amountOfWorkingThreads = amountOfThreads;
        // Let threads initialize first
        printStatus();
        std::this_thread::sleep_for(std::chrono::seconds(20));

        while (amountOfWorkingThreads > 0) {
            tempIndexATM = 0;
            amountOfWorkingThreads = 0;
            for (auto &s : Sims) {
                amountOfWorkingThreads += s.amountOfWorkingThreads;
                tempIndexATM += s.tempIndexATM;
            }
            tempIndexATM = std::min(tempIndexATM, temps.size() - 1);
            printStatus();
            std::this_thread::sleep_for(std::chrono::seconds(20));
        }

        for (auto &i : threads) {
            i.join();
        }

        for (const auto &Simulations:Sims) {
            //TODO check if measurements correspond to temps
#ifdef DEBUG
            std::cout << "Simulations.getTemps().size()=" << Simulations.getTemps().size() << std::endl;
#endif
            energies.insert(energies.end(), Simulations.getEnergies().begin(), Simulations.getEnergies().end());
            magnetization.insert(end(magnetization), Simulations.getMagnetization().begin(),
                                 Simulations.getMagnetization().end());
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

    /**
     * prints status of simulation to console
     */
    void printStatus() const {
        //TODO add shuffle for new ensemble

        // TODO add ETA

        const std::string sep = " | ";
        const std::string tempSize = std::string(std::to_string(temps.size()));

        const auto time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        const std::string timeString = std::string(std::ctime(&time));


        std::cout << sep << "N=" << std::left << std::setw(5) << sights
                  << sep << "run:" << std::right << std::setw(static_cast<int>(tempSize.size()))
                  << tempIndexATM + 1 << "/" << std::left << temps.size()

                  << sep << "T=" << std::setprecision(3) << std::setw(6) << temps[tempIndexATM]
                  << sep << std::setprecision(4) << std::setw(6)
                  << static_cast<float>(tempIndexATM * 100) / static_cast<float>(temps.size() - 1) << "%"

                  << sep << "threads: " << std::right << std::setw(3)
                  << amountOfWorkingThreads << "/" << std::left << amountOfThreads

                  << sep << timeString.substr(0, timeString.size() - 1)
                  << sep << '\n';
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

    /// Monitoring simulation parameters for std::cout
    unsigned long tempIndexATM;
    unsigned int amountOfThreads;
    unsigned int amountOfWorkingThreads;

public:
    bool printStat;

private:
    /// SpinLattice for simulation
    SpinLattice2level sl;
    bool isSimulated;
};
