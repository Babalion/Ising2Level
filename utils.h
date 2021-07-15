//
// Created by chris on 15.06.21.
//
#pragma once

#include <filesystem>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>

#include "SpinLattice2level.h"

/**
 * This is the runtime GUI that let's you watch what happens during the
 * optimization procedure
 */
class RuntimeGUI {
public:
    /**
     * Constructor
     */
    RuntimeGUI(int rows, int cols) : waitTime(25), gui(rows, cols, CV_8UC3) {
        // Open the window
        cv::namedWindow("GUI", 1);
    }

    /**
     * Destructor
     */
    virtual ~RuntimeGUI() {
        cv::destroyWindow("GUI");
    }

    /**
     * Paint the gui
     */
    void notify(const SpinLattice2level &instance);

    /**
     * The time the GUI pauses after each update. Set to 0 to let
     * it wait for a keypress
     */
    int waitTime;

    void saveImageToDisk(const std::string &imageName, const std::string &path);

private:
    /**
     * The GUI matrix
     */
    cv::Mat gui;
};


////////////////////////////////////////////////////////////////////////////////
/// RuntimeGUI
////////////////////////////////////////////////////////////////////////////////

void RuntimeGUI::notify(const SpinLattice2level &instance) {
    // The screen is split as follows:
    // 75% points
    // 25% status

    // Clear the gui
    gui = cv::Scalar(40, 40, 40);

    // Get the status marker
    int statusCol = 0.75 * gui.cols;

    // Write the status
    std::stringstream ss;
    ss << "energy = " << instance.calcEnergy();
    cv::putText(gui,
                ss.str(),
                cv::Point(statusCol, 15),
                cv::FONT_HERSHEY_PLAIN,
                0.9,
                cv::Scalar(255, 255, 255));
    ss.str("");
    ss << "sweeps = " << instance.performedSweeps;
    cv::putText(gui,
                ss.str(),
                cv::Point(statusCol, 30),
                cv::FONT_HERSHEY_PLAIN,
                0.9,
                cv::Scalar(255, 255, 255));
    ss.str("");
    /*ss << "inner = " << config.inner;
    cv::putText(gui,
                ss.str(),
                cv::Point(statusCol, 45),
                cv::FONT_HERSHEY_PLAIN,
                0.9,
                cv::Scalar(255, 255, 255));
    ss.str("");
    ss << "energy = " << config.energy;
    cv::putText(gui,
                ss.str(),
                cv::Point(statusCol, 60),
                cv::FONT_HERSHEY_PLAIN,
                0.9,
                cv::Scalar(255, 255, 255));
    ss.str("");
    ss << "best energy = " << config.bestEnergy;
    cv::putText(gui,
                ss.str(),
                cv::Point(statusCol, 75),
                cv::FONT_HERSHEY_PLAIN,
                0.9,
                cv::Scalar(255, 255, 255));

    ss.str("");
    ss << "epsilon =  = " << config.bestEnergy;
    cv::putText(gui,
                ss.str(),
                cv::Point(statusCol, 90),
                cv::FONT_HERSHEY_PLAIN,
                0.9,
                cv::Scalar(255, 255, 255));
*/

    // Plot the charts
    // [...]

    // Plot the cities
    // Determine the minimum and maximum X/Y
    float minX = 0;
    float minY = 0;
    auto maxX = static_cast<float>(instance.getSights());
    auto maxY = static_cast<float>(instance.getSights());


    // Calculate the compression factor
    float width = maxX - minX;
    float height = maxY - minY;
    float compression = (statusCol - 10) / width;
    if (height * compression > gui.rows - 10) {
        compression = (gui.rows - 10) / height;
    }

    // Paint the spins
    for (size_t i = 0; i < instance.getSights() * instance.getSights(); i++) {
        cv::Point p1;
        p1.x = static_cast<float>((i % instance.getSights())) * compression + 10;
        p1.y = static_cast<float>((i / instance.getSights())) * compression + 10;

        if (instance.getSpins()[i] == -1) {
            //Coler Ordering BGR (blue,green,red)
            cv::circle(gui, p1, 2, cv::Scalar(212, 188, 0), 5);
        } else {
            cv::circle(gui, p1, 2, cv::Scalar(255, 0, 255), 5);
        }
    }

    cv::imshow("GUI", gui);
    cv::waitKey(waitTime);
    /*
    if (config.terminated) {
        cv::waitKey(0);
    } else {
        cv::waitKey(waitTime);
    }*/
}

/**
 * Save the gui-frame as image to disk
 * @param imageName name to save the image to
 * @param pathname must end with a /
 */
void RuntimeGUI::saveImageToDisk(const std::string &imageName, const std::string &pathname) {
    struct stat info{};

    /// Try to create a directory to save images to
    if (stat(pathname.c_str(), &info) != 0) {
        std::cout << "cannot access, try to create now. " << pathname << std::endl;
        try {
            std::filesystem::create_directories(pathname);
        } catch (const std::exception &e) { // caught by reference to base
            std::cout << " a standard exception was caught, with message '"
                      << e.what() << "'\n";
        }
    } else if (info.st_mode & S_IFDIR) {  // S_ISDIR() doesn't exist on my windows
#ifdef DEBUG
        std::cout<<pathname<<" is a directory"<<std::endl;
#endif
    } else {
        std::cerr << pathname << "is no directory" << std::endl;
    }

    // compression params from https://docs.opencv.org/3.4/d4/da8/group__imgcodecs.html#gabbc7ef1aa2edfaa87772f1202d67e0ce
    static const std::vector<int> compression_params{cv::IMWRITE_PNG_COMPRESSION, 9};

    std::string filePath = pathname;
    filePath = filePath.append(imageName);

    cv::imwrite(filePath, gui, compression_params);
}

/**
 * @tparam T
 * @param vec
 * @return maximum element of given vector
 */
template<typename T>
T max(const std::vector<T> &vec) {
    T maxVal = vec[0];
    for (int i = 1; i < vec.size(); ++i) {
        maxVal = std::max(maxVal, vec[i]);
    }
    return maxVal;
}

/**
 * @tparam T
 * @param vec
 * @return minimum element of given vector
 */
template<typename T>
T min(const std::vector<T> &vec) {
    T minVal = vec[0];
    for (int i = 1; i < vec.size(); ++i) {
        minVal = std::min(minVal, vec[i]);
    }
    return minVal;
}

/**
 * Calculates the arithmetic mean of given vector-elements
 * @tparam T
 * @param vec
 * @return the arithmetic mean
 */
template<typename T>
T mean(const std::vector<T> &vec) {
    return std::accumulate(vec.begin(), vec.end(), 0.0)/vec.size();
}

/**
 * Calculates the arithmetic mean of given vector-of-vector-of-elements
 * @tparam T
 * @param vec vector of vector of elements (tensor of rank 2)
 * @return the arithmetic mean of all values
 */
template<typename T>
T mean(const std::vector<std::vector<T>> &vec) {
    std::vector<T> firstMeans;
    firstMeans.reserve(vec.size());
    for(auto &i:vec){
        firstMeans.push_back(mean(i));
    }
    return mean(firstMeans);
}

/**
 * Calculates the standard deviation of given vector-elements
 * @tparam T
 * @param vec
 * @return standard deviation
 */
template<typename T>
T stdev(const std::vector<T> &vec) {
    T summation=0;
    const T vec_mean=mean(vec);
    for (size_t i = 0; i < vec.size(); ++i) {
        summation+=pow(vec[i]-vec_mean,2);
    }
    return sqrt(summation*1.0/(vec.size()-1));
}

/**
 * Calculates the autocorrelation-function of given data-sample
 * @tparam T
 * @param vec of data-sample with length n
 * @return autocorrelation-function in vector with length n/2
 */
template<typename T>
std::vector<T> autoCorr(const std::vector<T> &vec) {
    float mean = std::accumulate(vec.begin(), vec.end(), 0.0f) / vec.size();

    std::vector<float> autocorrelation(vec.size() / 2);
    for (size_t t = 0; t < autocorrelation.size(); t++) {
        float n = 0; // Numerator
        float d = 0; // Denominator
        for (size_t i = 0; i < vec.size() - t; i++) {
            float xim = vec[i] - mean;
            n += xim * (vec[(i + t) % vec.size()] - mean);
            d += xim * xim;
        }
        autocorrelation[t] = n / d;
    }
    return autocorrelation;
}

/**
 * normalizes given vector call by reference
 * @tparam T
 * @param vec vector which should be normalized
 */
template<typename T>
void normalize(std::vector<T> &vec) {
    const T vec_max = max(vec);
    for (T &i : vec) {
        i /= 1.0 * vec_max;
    }
}


struct LegendLabel : public CvPlot::Drawable {
    std::string _text;
    cv::Point2f _position;
    const int _fontFace = cv::FONT_HERSHEY_SIMPLEX;
    const double _fontScale = .4;
    const int _fontThickness = 1;
    cv::Scalar _color = cv::Scalar(0, 0, 0);

    void render(CvPlot::RenderTarget &renderTarget) override {
        int baseline;
        cv::Size size = cv::getTextSize(_text, _fontFace, _fontScale, _fontThickness, &baseline);
        auto pos = renderTarget.innerToOuter(renderTarget.project(_position)) +
                   cv::Point2d(size.height * 2, size.height / 2);
        cv::putText(renderTarget.outerMat(), _text, pos, _fontFace, _fontScale, _color, _fontThickness, cv::LINE_AA);
    }
};

struct Legend : public CvPlot::Drawable {
    CvPlot::Axes *_parentAxes;
    int _width = 200;
    int _height = 70;
    int _margin = 20;

    void render(CvPlot::RenderTarget &renderTarget) override {
        std::vector<CvPlot::Series *> seriesVec;
        for (const auto &drawable : _parentAxes->drawables()) {
            auto series = dynamic_cast<CvPlot::Series *>(drawable.get());
            if (series) {
                seriesVec.push_back(series);
            }
        }
        CvPlot::Axes axes;
        axes.setMargins(5, _width - 2 * _margin - 60, 5, 5)
                .setXLim({-.2, 1.2})
                .setYLim({-.2, seriesVec.size() - 1 + .2})
                .setYReverse();
        for (size_t i = 0; i < seriesVec.size(); i++) {
            auto &series = *seriesVec[i];
            axes.create<CvPlot::Series>(std::vector<cv::Point2f>{{0,   (float) i},
                                                                 {0.2, (float) i}})
                    .setLineType(series.getLineType())
                    .setLineWidth(series.getLineWidth())
                    .setColor(series.getColor())
                    .setMarkerType(series.getMarkerType())
                    .setMarkerSize(series.getMarkerSize());
            auto &label = axes.create<LegendLabel>();
            label._position = {0.2, (float) i};
            label._text = series.getName();
        }
        cv::Rect rect(renderTarget.innerMat().cols - _width - _margin, _margin, _width, _height);
        if (rect.x >= 0 && rect.x + rect.width < renderTarget.innerMat().cols && rect.y >= 0 &&
            rect.y + rect.height < renderTarget.innerMat().rows) {
            axes.render(renderTarget.innerMat()(rect));
            cv::rectangle(renderTarget.innerMat(), rect, cv::Scalar::all(0));
        }
    }
};