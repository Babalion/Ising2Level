//
// Created by chris on 14.06.21.
//

#include "ising.h"

#include "SpinLattice2level.h"
#include <CvPlot/cvplot.h>

void simulateAndPlot() {
    const int noOfTemps = 8;
    const float tempStart = 2;
    const float tempEnd = 3;
    const int numIterations = 1000;
    const unsigned int shuffleAgainAfter = 100;

    std::vector<Simulation> Sim = {Simulation(2, noOfTemps, tempStart, tempEnd, numIterations, shuffleAgainAfter),
                                   Simulation(4, noOfTemps, tempStart, tempEnd, numIterations, shuffleAgainAfter),
                                   Simulation(8, noOfTemps, tempStart, tempEnd, numIterations, shuffleAgainAfter),
                                   Simulation(16, noOfTemps, tempStart, tempEnd, numIterations, shuffleAgainAfter)};

    /// Start simulation
    for (auto &s:Sim) {
        s.thermalizeSweeps = 1000;
        s.sweepsPerIteration = 2;
        std::cout << "Simulating now N=" << s.getSights() << "\n";
        s.simulate_par();
    }

    /// ---------------------------------------------------------------------------------------------------------------
    /// Plot and calculate measured parameters
    /// ---------------------------------------------------------------------------------------------------------------

    std::vector<std::string> colors = {"-r", "-g", "-b", "-c"};

    auto axesMagnetization = CvPlot::makePlotAxes();
    for (size_t i = 0; i < Sim.size(); i++) {
        axesMagnetization.create<CvPlot::Series>(Sim[i].getTemps(), Sim[i].getMagnetization(), colors[i])
                .setLineSpec("o").setMarkerSize(2)
                .setName(std::string("N^2=").append(std::to_string(Sim[i].getSights())));
    }
    axesMagnetization.xLabel("temperature T").yLabel("magnetization m");
    axesMagnetization.title("magnetization vs temperature");
    CvPlot::show("magnetization", axesMagnetization);

    auto axesEnergy = CvPlot::makePlotAxes();
    for (size_t i = 0; i < Sim.size(); i++) {
        axesMagnetization.create<CvPlot::Series>(Sim[i].getTemps(), Sim[i].getEnergies(), colors[i])
                .setLineSpec("o").setMarkerSize(2)
                .setName(std::string("N^2=").append(std::to_string(Sim[i].getSights())));
    }
    axesEnergy.xLabel("temperature T").yLabel("energy E");
    axesEnergy.title("energy vs temperature");
    CvPlot::show("energy", axesEnergy);
    /// ---------------------------------------------------------------------------------------------------------------
    /// Save measurements to file
    /// ---------------------------------------------------------------------------------------------------------------

    std::cout << "Simulation finished. Save results now...\n";
    std::ofstream file("IsingResults.tsv");
    file << "numOfTemps:\t" << noOfTemps << std::endl;
    file << "numOfIterations:\t" << numIterations << std::endl << std::endl;
    file << "N\ttemp\tmagnetization\tenergy\tsusceptibility\theatCapacity\n";
    file << std::fixed;
    file.precision(5);
    for (auto &s:Sim) {
        for (size_t i = 0; i < s.getTemps().size(); ++i) {
            file << s.getSights() << "\t" << s.getTemps()[i] << "\t"
                 << s.getMagnetization()[i] << "\t" << s.getEnergies()[i]
                 << "\t" << s.getSusceptibility()[i] << "\t" << s.getHeatCapacity()[i] << "\n";
        }
    }
    file.close();


}

int main() {
    simulateAndPlot();
}