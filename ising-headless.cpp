//
// Created by chris on 14.06.21.
//
#include "Simulation.h"

/** TASK 1:
 *
 * Simulate Systems with sight-size 128,256,512,1024
 * Reconstruct for each size dependent on temperature
 *  -energy E(T)
 *  -magnetization m(T)
 *  -heat capacity C(T)
 *  -susceptibility chi(T)
 *
 *  We need at least 1E4 measurements per temperature
 *  Calculate the autocorrelations
 */
void simulateAndPlot() {
    const float TCritical = 2.0f / std::log(1.0f + std::sqrt(2.0f));
    const float minTemp = 2;
    const float maxTemp = 2.6;
    const int numOfTemps = 8;
    const int numIterations = 1000;
    const unsigned int shuffleAgainAfter = UINT32_MAX;

    std::vector<Simulation> Sims = {
            Simulation(16, 8, minTemp, maxTemp, numIterations, shuffleAgainAfter),
            Simulation(32, 8, TCritical - 0.8f, TCritical + 0.8f, numIterations, shuffleAgainAfter),
            Simulation(64, 8, TCritical - 0.5f, TCritical + 0.5f, numIterations, shuffleAgainAfter),
            Simulation(128, 8, TCritical - 0.1f, TCritical + 0.1f, numIterations, shuffleAgainAfter),
            Simulation(256, 8, TCritical - 0.1f, TCritical + 0.1f, numIterations, shuffleAgainAfter)};

    for (auto &S:Sims) {
        S.sweepsPerIteration = 5;
        S.thermalizeSweeps = 50;
    }

    std::ofstream file("IsingResultsWolff.tsv");
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