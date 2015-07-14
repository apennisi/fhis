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
* \file delaunaytri.h
* 
* \class DelaunayTri
* 
* \brief Class for calculating Delaunay Triangulation
* 
**/ 

#ifndef DELAUNAYTRI_H
#define DELAUNAYTRI_H

#include <opencv2/opencv.hpp>
#include <map>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <math.h>
#include "utils.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                       Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
typedef K::Point_2 PointCGAL;

class DelaunayTri
{
public:
    /**
    * \brief Create a new object DelaunaTri
    *
    * \param XYpoints edge points
    * \param width image width
    * \param height image height
    * \param size XYPoints size
    *
    */
    DelaunayTri(const std::vector<cv::Point2i> &XYpoints, const int &width,
                    const int &height, const int &size);
    virtual ~DelaunayTri();
    /**
    * \brief return the indices of the triangulation
    *
    * \return return the indices of the triangulation
    */
    inline std::vector<int> getIndices(){ return indices; }
    /**
    * \brief return the radii of the circles that circumscribe the triangles
    *
    * \return return the radii of the circles that circumscribe the triangles
    */
    inline std::vector<double> getRadii() { return radii; }
    /**
    * \brief return the centers of the circles that circumscribe the triangles
    *
    * \return return the centers of the circles that circumscribe the triangles
    */
    inline std::vector<cv::Point2d> getCenters() { return centers; }
    /**
    * \brief return the triangle neighbor list
    *
    * \return return the triangle neighbor list
    */
    inline std::vector<double> getNeighbors () { return neighbors; }
    /**
    * \brief return the size of the triangle list
    *
    * \return return the size of the triangle list
    */
    inline int getSize() { return size_; }
private:
    Delaunay dt;
    void draw_subdiv();
    cv::Mat img;
    int size_;
    std::vector<cv::Vec6d> triangleList;
    std::vector<int> indices;
    std::vector<cv::Point2d> centers;
    std::vector<double> radii;
    int width, height;
    std::vector< std::pair<PointCGAL,unsigned> > points;
    /**
    * \brief insert a new index in the map
    *
    * \param ind current index
    * \param second previous index
    * 
    */
    void insertMap(const int *ind, const int &second);
    /**
    * \brief compute the circumcenter of the circle that circumscribes a triangle
    *
    * \param a triangle vertex
    * \param b triangle vertex
    * \param c traingle vertex
    * 
    * \return return the circumcenter
    */
    cv::Point2d GetCircumcenter(const cv::Point2d &a, const cv::Point2d &b,
                                       const cv::Point2d &c);
    std::map<std::vector<int>, std::vector<int> > mapTriangles;
    /**
    * \brief compute the neighbors of each triangle
    * 
    */
    void computeNeighbors();
    std::vector<double> neighbors;
};

#endif // DELAUNAY_H
