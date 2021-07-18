//
// Created by chris on 14.06.21.
//

#include "ising.h"
#include "utils.h"

#include "SpinLattice2level.h"
#include <CvPlot/cvplot.h>

void simulateAndPlot() {
    const int noOfTemps = 8;
    const float tempStart = 2;
    const float tempEnd = 3;
    const int numIterations = 1000;
    const unsigned int shuffleAgainAfter = UINT32_MAX;

    std::vector<Simulation> Sim = {Simulation(16, noOfTemps, tempStart, tempEnd, numIterations, shuffleAgainAfter),
                                   Simulation(32, noOfTemps, tempStart, tempEnd, numIterations, shuffleAgainAfter),
                                   Simulation(64, noOfTemps, tempStart, tempEnd, numIterations, shuffleAgainAfter),
                                   Simulation(128, noOfTemps, tempStart, tempEnd, numIterations, shuffleAgainAfter)};

    /// Start simulation
    for (auto &s:Sim) {
        s.thermalizeSweeps = 50;
        s.sweepsPerIteration = 2;
        std::cout << "Simulating now N=" << s.getSights() << "\n";
        s.simulate_seq();
    }

    for (const auto &i : Sim[1].getEnergies()) {
        std::cout<<i<<std::endl;
    }

    /// ---------------------------------------------------------------------------------------------------------------
    /// Plot and calculate measured parameters
    /// ---------------------------------------------------------------------------------------------------------------

    std::vector<std::string> colors = {"-r", "-g", "-b", "-c"};

    auto axesMagnetization = CvPlot::makePlotAxes();
    for (size_t i = 0; i < Sim.size(); i++) {
        axesMagnetization.create<CvPlot::Series>(Sim[i].getTemps(), Sim[i].getMagnetization(), colors[i])
                .setLineSpec("o").setMarkerSize(2)
                .setName(std::string("N=").append(std::to_string(Sim[i].getSights())).append("^2"));
    }
    axesMagnetization.xLabel("temperature T").yLabel("magnetization m");
    axesMagnetization.create<Legend>()._parentAxes = &axesMagnetization;
    axesMagnetization.title("magnetization vs temperature");
    CvPlot::show("magnetization", axesMagnetization);

    auto axesEnergy = CvPlot::makePlotAxes();
    for (size_t i = 0; i < Sim.size(); i++) {
        axesEnergy.create<CvPlot::Series>(Sim[i].getTemps(), Sim[i].getEnergies(), colors[i])
                .setLineSpec("o").setMarkerSize(2)
                .setName(std::string("N=").append(std::to_string(Sim[i].getSights())).append("^2"));
    }
    axesEnergy.xLabel("temperature T").yLabel("energy E");
    axesEnergy.create<Legend>()._parentAxes = &axesEnergy;
    axesEnergy.title("energy vs temperature");
    CvPlot::show("energy", axesEnergy);
    /// ---------------------------------------------------------------------------------------------------------------
    /// Save measurements to file
    /// ---------------------------------------------------------------------------------------------------------------

    std::cout << "Simulation finished. Save results now...\n";
    std::ofstream file("IsingResultsTemp.tsv");
    file << "numOfTemps:\t" << noOfTemps << std::endl;
    file << "numOfIterations:\t" << numIterations << std::endl << std::endl;
    file << "N\ttemp\tmagnetization\tenergy\n";
    file << std::fixed;
    file.precision(5);
    for (auto &s:Sim) {
        for (size_t i = 0; i < s.getTemps().size(); ++i) {
            file << s.getSights() << "\t" << s.getTemps()[i] << "\t"
                 << s.getMagnetization()[i] << "\t" << s.getEnergies()[i]<< "\n";
        }
    }
    file.close();


}

int main() {
    simulateAndPlot();
}