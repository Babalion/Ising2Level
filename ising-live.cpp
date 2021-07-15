//
// Created by chris on 16.07.21.
//
#include "ising.h"
#include "utils.h"


void runIsingLive() {
    SpinLattice2level sl(75);
    // create arrays to save measurements to
    std::vector<std::vector<float>> lowTemp;
    std::vector<std::vector<float>> lowTempHB;

    // Register the GUI
    // You can specify the dimensions of the window
    RuntimeGUI gui(1000, 1200);
    gui.waitTime = 1;
    gui.notify(sl);

    const int numEnsembles = 10;
    const int numSteps = 1000;
    // reserve memory
    lowTemp.reserve(numEnsembles);
    lowTempHB.reserve(numEnsembles);


    /// Start simulation
    for (int i = 0; i < numEnsembles; i++) {
        lowTemp.push_back(std::vector<float>{(float) sl.calcEnergy()});
        for (int j = 0; j < numSteps; ++j) {
            metropolisSweep(sl, 2, 1);
            lowTemp[i].push_back(sl.calcEnergy());
            if (j == 0)
                lowTemp[i].reserve(numSteps);
            if (j % 20 == 0)
                gui.notify(sl);
        }
        sl.initRandom();
    }

    sl.initRandom();
    for (int i = 0; i < numEnsembles; i++) {
        lowTempHB.push_back(std::vector<float>{(float) sl.calcEnergy()});
        for (int j = 0; j < numSteps; ++j) {
            heatBathSweep(sl, 2);
            lowTempHB[i].push_back(sl.calcEnergy());
            if (i % 20 == 0)
                gui.notify(sl);
        }
        sl.initRandom();
    }


    /// we can save a screenshot of the gui to disk (so you also can generate videos with ffmpeg)
    //gui.saveImageToDisk("2levelIsing.png", "./results/");

    /// ---------------------------------------------------------------------------------------------------------------
    /// Plot and calculate measured energy
    /// ---------------------------------------------------------------------------------------------------------------

    auto axes = CvPlot::makePlotAxes();
    for(size_t i=0;i<lowTemp.size();i++) {
        axes.create<CvPlot::Series>(lowTemp[i], "-g").setMarkerType(CvPlot::MarkerType::Circle).setMarkerSize(
                10).setLineType(CvPlot::LineType::None);
        axes.create<CvPlot::Series>(lowTempHB[i], "-r").setMarkerType(CvPlot::MarkerType::Circle).setMarkerSize(
                10).setLineType(CvPlot::LineType::None);
    }
    std::cout << "mean(lowT)= " << mean(lowTemp) << std::endl;
    std::cout << "mean(lowT_HB)= " << mean(lowTempHB) << std::endl;

    CvPlot::show("plotPhysics", axes);

    /// ---------------------------------------------------------------------------------------------------------------
    /// Calculate and plot autocorrelation
    /// ---------------------------------------------------------------------------------------------------------------

    auto axesAC = CvPlot::makePlotAxes();
    for (int i = 0; i < numEnsembles; ++i) {
        std::vector<float> ac_le = autoCorr(lowTemp[i]);
        normalize(ac_le);
        std::vector<float> acHB_le = autoCorr(lowTempHB[i]);
        normalize(acHB_le);
        if(i==0) {
            axesAC.create<CvPlot::Series>(ac_le, "-b").setName("Metropolis (T=0.1)");
            axesAC.create<CvPlot::Series>(acHB_le, "-g").setName("Heatbath (T=0.1)");
        }
    }
    axesAC.xLabel("simulation steps n").yLabel("normalized autocorrelation");
    axesAC.title("Autocorrelation Metropolis vs. Heatbath in 2D Ising");
    axesAC.create<Legend>()._parentAxes = &axesAC;
    CvPlot::show("autocorr", axesAC);

    cv::waitKey(0);
}

int main() {
    runIsingLive();
}