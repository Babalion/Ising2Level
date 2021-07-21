//
// Created by chris on 14.06.21.
//
#include "Simulation.h"
#include <sys/resource.h>

/** TASK 1:
 *
 * Simulate Systems with sight-size 128,256,512,1024
 * Reconstruct for each size dependent on temperature
 *  -energy E(T)
 *  -magnetization m(T)
 *
 *  We need at least 1E4 measurements per temperature for less autocorrelations
 */
void simulateAndPlot() {
    const float TCritical = 2.0f / std::log(1.0f + std::sqrt(2.0f));
    const float minTemp = 2;
    const float maxTemp = 2.6;
    const int numOfTemps = 16;
    const int numIterations = 1E5;
    const unsigned int shuffleAgainAfter = UINT32_MAX;

    std::vector<Simulation> Sims = {
            //Simulation(16, numOfTemps, 1.5, 3.5, numIterations, shuffleAgainAfter),
            //Simulation(32, numOfTemps, 1.6, 3, numIterations, shuffleAgainAfter),
            //Simulation(64, numOfTemps, 1.6, 3, numIterations, shuffleAgainAfter),
            //Simulation(128, numOfTemps, 2, 2.6, numIterations, shuffleAgainAfter),
            //Simulation(256, numOfTemps, 2, 2.4, numIterations, shuffleAgainAfter),
            //Simulation(512, numOfTemps, 2.1, 2.4, 2 * numIterations, shuffleAgainAfter),
            Simulation(1024, numOfTemps, 2.1, 2.4, 2 * numIterations, shuffleAgainAfter)};

    for (auto &S:Sims) {
        S.sweepsPerIteration = 5;
        S.thermalizeSweeps = 200;
    }

    std::ofstream file("IsingResultsWolff1024.tsv");
    file << "numOfTemps:\t" << numOfTemps << std::endl;
    file << "numOfIterations:\t" << numIterations << std::endl << std::endl;
    file << "N\ttemp\tmagnetization\tenergy\tsusceptibility\theatCapacity\n";
    file << std::fixed;
    file.precision(10);
    for (auto &S:Sims) {
        std::cout << "Simulating now N=" << S.getSights() << "\n";
        S.simulate_par();
        std::cout << "Simulation finished. Save results now...\n";

        /// Save measurements to file
        for (size_t i = 0; i < S.getNumOfTemps() * S.getNumOfIterations(); ++i) {
            file << S.getSights() << "\t" << S.getTemps()[i] << "\t"
                 << S.getMagnetization()[i] << "\t" << S.getEnergies()[i] << "\n";
        }
    }

    file.close();


}

int main() {
    simulateAndPlot();
}