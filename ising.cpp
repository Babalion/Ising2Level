//
// Created by chris on 14.06.21.
//

#include "ising.h"

#include "SpinLattice2level.h"
#include <CvPlot/cvplot.h>
#include "utils.h"


#ifdef DEBUG
void runTest() {
    SpinLattice2level sl(75);
    // create arrays to save measurements to
    std::vector<std::vector<float>> lowTemp;
    std::vector<std::vector<float>> lowTempHB;

    // Register the GUI
    // You can specify the dimensions of the window
    RuntimeGUI gui(1000, 1200);
    gui.waitTime = 1;
    gui.notify(sl);

    const int numEnsembles = 100;
    const int numSteps = 100;
    // reserve memory
    lowTemp.reserve(numEnsembles);
    lowTempHB.reserve(numEnsembles);/*
    for (int i = 0; i < numEnsembles; ++i) {
        lowTemp[i].reserve(numSteps);
        lowTempHB[i].reserve(numSteps);
    }*/


    /// Start simulation
    for (int i = 0; i < numEnsembles; i++) {
        lowTemp.push_back(std::vector<float>{(float) sl.calcEnergy()});
        for (int j = 0; j < numSteps; ++j) {
            metropolisSweep(sl, 2, 1);
            lowTemp[i].push_back(sl.calcEnergy());
            if (j == 0)
                lowTemp[i].reserve(numSteps);
            if (j % 1000 == 0)
                gui.notify(sl);
        }
        sl = SpinLattice2level(75);
    }

    sl = SpinLattice2level(75);
    for (int i = 0; i < numEnsembles; i++) {
        lowTempHB.push_back(std::vector<float>{(float) sl.calcEnergy()});
        for (int j = 0; j < numSteps; ++j) {
            heatBathSweep(sl, 2);
            lowTempHB[i].push_back(sl.calcEnergy());
            //if (i % 100 == 0)
            //gui.notify(sl);
        }
        sl = SpinLattice2level(75);
    }


    /// we can save a screenshot of the gui to disk (so you also can generate videos with ffmpeg)
    //gui.saveImageToDisk("2levelIsing.png", "./results/");

    /// ---------------------------------------------------------------------------------------------------------------
    /// Plot and calculate measured energy
    /// ---------------------------------------------------------------------------------------------------------------
    /*
    std::vector<float> x(lowTemp.size());
    for (size_t i = 0; i < x.size(); i++) {
        x[i] = static_cast<float>(i);
    }
    auto axes = CvPlot::makePlotAxes();
    axes.create<CvPlot::Series>(x, lowTemp, "-g").setMarkerType(CvPlot::MarkerType::Circle).setMarkerSize(
            10).setLineType(CvPlot::LineType::None);
    axes.create<CvPlot::Series>(x, lowTempHB, "-r").setMarkerType(CvPlot::MarkerType::Circle).setMarkerSize(
            10).setLineType(CvPlot::LineType::None);
*/
    std::cout << "mean(lowT)= " << mean(lowTemp) << std::endl;
    std::cout << "mean(lowT_HB)= " << mean(lowTempHB) << std::endl;

    //CvPlot::show("plotPhysics", axes);

    /// ---------------------------------------------------------------------------------------------------------------
    /// Calculate and plot autocorrelation
    /// ---------------------------------------------------------------------------------------------------------------

    auto axesAC = CvPlot::makePlotAxes();
    for (int i = 0; i < numEnsembles; ++i) {
        std::vector<float> ac_le = autoCorr(lowTemp[i]);
        normalize(ac_le);
        std::vector<float> acHB_le = autoCorr(lowTempHB[i]);
        normalize(acHB_le);

        axesAC.create<CvPlot::Series>(ac_le, "-b").setName("Metropolis (T=0.1)");
        axesAC.create<CvPlot::Series>(acHB_le, "-g").setName("Heatbath (T=0.1)");
    }
    axesAC.xLabel("simulation steps n").yLabel("normalized autocorrelation");
    axesAC.title("Autocorrelation Metropolis vs. Heatbath in 2D Ising");
    axesAC.create<Legend>()._parentAxes = &axesAC;
    CvPlot::show("autocorr", axesAC);

    cv::waitKey(0);
}
#endif

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