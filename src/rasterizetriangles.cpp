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


#include "rasterizetriangles.h"

RasterizeTriangles::RasterizeTriangles(const int &rows, const int &cols, const std::vector<cv::Point2i> &XYPoints,
                                       const std::vector<int> &indices, const int &size)
{
    this->nrows = rows;
    this->ncols = cols;

    rasterizedImage.resize(rows * cols, 0);

    min_col_buffer.resize(nrows);
    max_col_buffer.resize(nrows);

    int i;

    for(i = 0; i < size; ++i) {
        int a = indices[i*3];
        int b = indices[i*3 + 1];
        int c = indices[i*3 + 2];

        // Extract the vertices
        cv::Point2i P1 = cv::Point2i( XYPoints[a - 1].x - 1,
                                     XYPoints[a - 1].y - 1);

        cv::Point2i P2 = cv::Point2i( XYPoints[b - 1].x - 1,
                                     XYPoints[b - 1].y - 1);

        cv::Point2i P3 = cv::Point2i( XYPoints[c - 1].x - 1,
                                     XYPoints[c - 1].y - 1);


        // Make the edges
        edge edge1 = makeEdge (P1, P2);
        edge edge2 = makeEdge (P2, P3);
        edge edge3 = makeEdge (P3, P1);

        // Draw the triangle
        rasterizeTriangle (&edge1, &edge2, &edge3, i + 1);
   }

}

edge RasterizeTriangles::makeEdge(const cv::Point2i &P, const cv::Point2i &Q) {

    edge theEdge;

    if (P.y < Q.y) {
        theEdge.P1 = P;
        theEdge.P2 = Q;
    } else {
        theEdge.P1 = Q;
        theEdge.P2 = P;
    }

    return theEdge;
}


RasterizeTriangles::~RasterizeTriangles() {

}


int RasterizeTriangles::edgeHeight (edge *theEdge) {
    return (theEdge->P2.y - theEdge->P1.y);
}

void RasterizeTriangles::updateLineExtents (const int &row, const int &col)
{
    if ( (row >= 0) && (row < nrows) ) {
        if (col < min_col_buffer[row]) min_col_buffer[row] = col;
        if (col > max_col_buffer[row]) max_col_buffer[row] = col;
    }

}

void RasterizeTriangles::DrawSpans (const int &row1, const int &row2, const int &label)
{


    int row, col, min_col, max_col;

    for (row = row1; row <= row2; ++row) {
        min_col = min_col_buffer[row];
        max_col = max_col_buffer[row];



        if (min_col < 0) min_col = 0;
        if (max_col >= ncols) max_col = ncols-1;



        for (col = min_col; col <= max_col; ++col) {
            rasterizedImage[row*ncols + col] = (double)label;
        }
    }
}

void RasterizeTriangles::rasterizeEdge(edge *theEdge) {
    int x1 = theEdge->P1.x;
    int y1 = theEdge->P1.y;

    int x2 = theEdge->P2.x;
    int y2 = theEdge->P2.y;

    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int length, i, x, y;

    length = (dx > dy) ? dx : dy;

    if (length == 0) {

        updateLineExtents (y1, x1);

    } else {

        int doubleLength = 2*length;

        for (i = 0; i <= length; ++i) {

            x = ( 2*(x1*i + x2*(length-i)) + length ) / (doubleLength);
            y = ( 2*(y1*i + y2*(length-i)) + length ) / (doubleLength);

            updateLineExtents (y, x);
        }

    }
}


void RasterizeTriangles::rasterizeTriangle(edge *edge1, edge *edge2, edge *edge3, const int &label)
{
    int height1 = edgeHeight(edge1);
    int height2 = edgeHeight(edge2);
    int height3 = edgeHeight(edge3);
    edge *tall_edge, *short_edge1, *short_edge2;
    int i, row1, row2;

    // Find the tallest edge
    if ( (height1 >= height2) && (height1 >= height3) ) {
        tall_edge   = edge1;
        short_edge1 = edge2;
        short_edge2 = edge3;
    } else if ( (height2 >= height1) && (height2 >= height3) ) {
        tall_edge   = edge2;
        short_edge1 = edge1;
        short_edge2 = edge3;
    } else {
        tall_edge   = edge3;
        short_edge1 = edge2;
        short_edge2 = edge1;

    }

    row1 = tall_edge->P1.y;
    row2 = tall_edge->P2.y;


    if (row1 < 0) row1 = 0;
    if (row2 >= nrows) row2 = nrows-1;


    // Clear the min and max idx buffers
    for (i = row1; i <= row2; ++i) {
        min_col_buffer[i] = ncols + 1;
        max_col_buffer[i] = -1;
    }

    // render the edges to get the triangle extents on each row
    rasterizeEdge (edge1);
    rasterizeEdge (edge2);
    rasterizeEdge (edge3);

    // Go through and fill in the spans
    DrawSpans (row1, row2, label);
}
