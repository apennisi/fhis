/*
 *  FHIS Image Segmentation Library
 *
 *
 *  Written by Andrea Pennisi
 *
 *  Please, report suggestions/comments/bugs to
 *  andrea.pennisi@gmail.com
 *
 */

/** 
* \file edgedetector.h
* 
* \class EdgeDetectir
* 
* \brief Class for detecting edges
* 
**/ 

#ifndef EDGEDETECTOR_H
#define EDGEDETECTOR_H

#include <vector>
#include <queue>
#include <iostream>
#include <math.h>
#include <opencv2/opencv.hpp>

#include "utils.h"


class EdgeDetector {

    typedef struct {
        int i, j;
    } Position;

public:
    /**
    * \brief Create a new object EdgeDetector
    *
    * \param src source image
    * \param edgeImg the image where the edges are stored
    * \param edge_thresh filter threshold
    * \param edge_sigma filter threshold
    *
    */
    EdgeDetector( const cv::Mat &src, cv::Mat &edgeImg, const double &edge_thresh, const double &edge_sigma );
    virtual ~EdgeDetector();
    /**
    * \brief Return the edge points
    *
    * \return return edge epoints
    *
    */
    std::vector<cv::Point2i> getEdgesPoints() { return points; }
private:
    void computGradient(const double &edge_sigma);
    bool im2single(const cv::Mat &src, std::vector<double> &mat );
    bool imfilter(const std::vector<double> &src, const std::vector<double> &filter,
                  std::vector<double> &_mat, const bool &isTranspose);
    bool smoothGradient(const std::vector<double> &src, std::vector<double> &dx,
                        std::vector<double> &dy);
    void hypot(const std::vector<double> &dx, const std::vector<double> &dy,
               std::vector<double> &magGrad);
    void selectThresholds(double& lowThresh, double& highThresh, const double &thresh);
    bool cannyFindLocalMaxima(const std::vector<double> &src, const std::vector<double> &dx,
                              const std::vector<double> &dy, const cv::Point2i &pos);
    void thinAndThreshold(const std::vector<double> &src, const std::vector<double> &dx, const std::vector<double> &dy,
                           const double& lowThresh, const double& highThresh, std::vector<uchar> &dst );
    void bwselect(std::vector<uchar> &idxWeak, std::vector<uchar> &idxWeakNeg,
                  const std::vector<cv::Point2i> &idxStrongPts, const int &size);
    void dtMat2Img(const std::vector<uchar> &mat, cv::Mat &img);
    double _eps;
    double PercentOfPixelsNotEdges;
    double ThresholdRatio;
    int dimKernel;
    std::vector<double> gaussKernel;
    std::vector<double> derivGaussKernel;
    int rows, cols, numPixels;
    std::vector<cv::Point2i> points;
    void gradient();
};

#endif
