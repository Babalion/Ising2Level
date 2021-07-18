//
// Created by chris on 16.07.21.
//
#include "ising.h"
#include "utils.h"


void runIsingLive() {
    SpinLattice2level sl(100);
    // create arrays to save measurements to
    std::vector<std::vector<float>> energyMetro;
    std::vector<std::vector<float>> energyHB;
    std::vector<std::vector<float>> energyWolff;

    // Register the GUI
    // You can specify the dimensions of the window
    RuntimeGUI gui(1000, 1200);
    gui.waitTime = 1;
    gui.notify(sl);

    const int numEnsembles = 1;
    const int numSteps = 1000;
    const float temp = 2.0f / std::log(1.0f + std::sqrt(2.0f));

    // reserve memory
    energyMetro.reserve(numEnsembles);
    energyHB.reserve(numEnsembles);
    energyWolff.reserve(numEnsembles);


    /// Start simulation
    for (int i = 0; i < numEnsembles; i++) {
        sl.initCold();
        energyMetro.push_back(std::vector<float>{(float) sl.calcEnergy()});
        for (int j = 0; j < numSteps; ++j) {
            metropolisSweep(sl, temp, 10);
            if (j == 0)
                energyMetro[i].reserve(numSteps);
            energyMetro[i].push_back(sl.calcEnergy());
            if (j % 20 == 0)
                gui.notify(sl);
        }
    }

    for (int i = 0; i < numEnsembles; i++) {
        sl.initCold();
        energyHB.push_back(std::vector<float>{(float) sl.calcEnergy()});
        for (int j = 0; j < numSteps; ++j) {
            heatBathSweep(sl, temp, 10);
            if (j == 0)
                energyHB[i].reserve(numSteps);
            energyHB[i].push_back(sl.calcEnergy());
            if (j % 20 == 0)
                gui.notify(sl);
        }
    }

    for (int i = 0; i < numEnsembles; i++) {
        sl.initCold();
        energyWolff.push_back(std::vector<float>{(float) sl.calcEnergy()});
        for (int j = 0; j < numSteps; ++j) {
            wolffSweep(sl, temp, 10);
            if (j == 0)
                energyWolff[i].reserve(numSteps);
            energyWolff[i].push_back(sl.calcEnergy());
            if (j % 20 == 0)
                gui.notify(sl);
        }
    }


    /// we can save a screenshot of the gui to disk (so you also can generate videos with ffmpeg)
    //gui.saveImageToDisk("2levelIsing.png", "./results/");

    /// ---------------------------------------------------------------------------------------------------------------
    /// Plot and calculate measured energy
    /// ---------------------------------------------------------------------------------------------------------------

    auto axesEnergy = CvPlot::makePlotAxes();
    for (size_t i = 0; i < energyMetro.size(); i++) {
        if (i == 0) {
            axesEnergy.create<CvPlot::Series>(energyMetro[i], "-g").setMarkerType(
                    CvPlot::MarkerType::Circle).setMarkerSize(
                    10).setLineType(CvPlot::LineType::None).setName("Metropolis (T=T_c)");

            axesEnergy.create<CvPlot::Series>(energyHB[i], "-r").setMarkerType(
                    CvPlot::MarkerType::Circle).setMarkerSize(
                    10).setLineType(CvPlot::LineType::None).setName("HeatBath (T=T_c)");

            axesEnergy.create<CvPlot::Series>(energyWolff[i], "-b").setMarkerType(
                    CvPlot::MarkerType::Circle).setMarkerSize(
                    10).setLineType(CvPlot::LineType::None).setName("Wolff (T=T_c)");
        } else {
            axesEnergy.create<CvPlot::Series>(energyMetro[i], "-g").setMarkerType(
                    CvPlot::MarkerType::Circle).setMarkerSize(
                    10).setLineType(CvPlot::LineType::None);

            axesEnergy.create<CvPlot::Series>(energyHB[i], "-r").setMarkerType(
                    CvPlot::MarkerType::Circle).setMarkerSize(
                    10).setLineType(CvPlot::LineType::None);

            axesEnergy.create<CvPlot::Series>(energyWolff[i], "-b").setMarkerType(
                    CvPlot::MarkerType::Circle).setMarkerSize(
                    10).setLineType(CvPlot::LineType::None);
        }

    }
    axesEnergy.xLabel("simulation steps n").yLabel("normalized energies");
    axesEnergy.title("Energies algorithms in 2D Ising");
    axesEnergy.create<Legend>()._parentAxes = &axesEnergy;
    CvPlot::show("plotPhysics", axesEnergy);


    //
    std::cout << "mean(metropolis)= " << mean(energyMetro) << std::endl;
    std::cout << "mean(heatBath)= " << mean(energyHB) << std::endl;
    std::cout << "mean(wolff)= " << mean(energyWolff) << std::endl;

    /// ---------------------------------------------------------------------------------------------------------------
    /// Calculate and plot autocorrelation
    /// ---------------------------------------------------------------------------------------------------------------

    auto axesAC = CvPlot::makePlotAxes();
    for (int i = 0; i < numEnsembles; ++i) {
        std::vector<float> ac_metro = autoCorr(energyMetro[i]);
        normalize(ac_metro);
        std::vector<float> ac_HB = autoCorr(energyHB[i]);
        normalize(ac_HB);
        std::vector<float> ac_wolff = autoCorr(energyWolff[i]);
        normalize(ac_wolff);
        if (i == 0) {
            axesAC.create<CvPlot::Series>(ac_metro, "-g").setName("Metropolis (T=T_c)");
            axesAC.create<CvPlot::Series>(ac_HB, "-r").setName("HeatBath (T=T_c)");
            axesAC.create<CvPlot::Series>(ac_wolff, "-b").setName("Wolff (T=T_c)");
        }
    }
    axesAC.xLabel("simulation steps n").yLabel("normalized autocorrelation");
    axesAC.title("Autocorrelations algorithms in 2D Ising");
    axesAC.create<Legend>()._parentAxes = &axesAC;
    CvPlot::show("autocorr", axesAC);

    cv::waitKey(0);
}

int main() {
    runIsingLive();
}