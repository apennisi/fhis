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

#ifndef RASTERIZETRIANGLES_H
#define RASTERIZETRIANGLES_H

#include "utils.h"
#include <opencv2/opencv.hpp>

typedef struct {
    cv::Point2i P1, P2;
} edge;


class RasterizeTriangles
{
public:
    RasterizeTriangles(const int &rows, const int &cols, const std::vector<cv::Point2i> &XYPoints,
                       const std::vector<int> &indices, const int &size);
    inline std::vector<double> getRasterizedImage() { return rasterizedImage; }
    virtual ~RasterizeTriangles();
private:
    edge makeEdge (const cv::Point2i &P, const cv::Point2i &Q);
    int edgeHeight (edge *theEdge);
    std::vector<double> rasterizedImage;
    std::vector<int> min_col_buffer;
    std::vector<int> max_col_buffer;
    void updateLineExtents (const int &row, const int &col);
    void DrawSpans (const int &row1, const int &row2, const int &label);
    void rasterizeEdge(edge *theEdge);
    void rasterizeTriangle (edge *edge1, edge *edge2, edge *edge3, const int &label);
    int nrows, ncols;
};

#endif // RASTERIZETRIANGLES_H
