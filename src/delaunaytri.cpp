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

#include "delaunaytri.h"


DelaunayTri::DelaunayTri(const std::vector<cv::Point2i> &XYpoints,
                         const int &width, const int &height, const int &size)
{

    this->width = width;
    this->height = height;

    img = cv::Mat::zeros(cv::Size(width, height), CV_8UC1);


    int i;
    points.resize(size);

    for(i = 0; i < size; ++i)
    {
        points[i] = std::make_pair(PointCGAL(XYpoints[i].x, XYpoints[i].y), i + 1);
    }

    dt.insert(points.begin(),points.end());

    int idxSize = 0;
    size_ = dt.number_of_faces();

    indices.resize(size_*3, -1);
    centers.resize(size_, cv::Point2d(0, 0));
    neighbors.resize(size_*3, 0.0);
    radii.resize(size_, 0.0);

    for(Delaunay::Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++) {

        cv::Point2d t[3];
        int idx[3];

        Delaunay::Face_handle face = it;

        int offset = idxSize * 3;
        for(i = 0; i < 3; i++) {
            t[i] = cv::Point2d(dt.triangle(face)[i].x(), dt.triangle(face)[i].y());
            *(idx + i) = face->vertex(i)->info();
            indices [offset + i] = face->vertex(i)->info();
        }

        this->insertMap(idx, idxSize);

        cv::Point2d circumCenter = this->GetCircumcenter(t[0], t[1], t[2]);

        centers[idxSize] = circumCenter;

        radii[idxSize] = distPoint2Point<double>(circumCenter.x, circumCenter.y,
                                                     t[1].x,  t[1].y);

        idxSize++;
    }

    //this->draw_subdiv();

    //cv::imwrite("delaunay.png", img);

    this->computeNeighbors();

}

DelaunayTri::~DelaunayTri() {

    /*if(!deallocate<int>(indices)) {
        std::cerr << "Error in deallocating in indices" << std::endl;
        exit(-1);
    }

    if(!deallocate<double>(neighbors)) {
        std::cerr << "Error in deallocating in neighbors" << std::endl;
        exit(-1);
    }

    if(!deallocate<Point<double> >(centers)) {
        std::cerr << "Error in deallocating in centers" << std::endl;
        exit(-1);
    }

    if(!deallocate<double>(radii)) {
        std::cerr << "Error in deallocating in radii" << std::endl;
        exit(-1);
    }*/

}


void DelaunayTri::insertMap(const int *ind, const int &secondElem) {

    int first = *ind, second = *(ind + 1), third = *(ind + 2);

    if(mapTriangles.size() == 0) {
        std::vector<int> edge;
        std::vector<int> idx;
        idx.push_back(secondElem);
        edge.push_back(first);
        edge.push_back(second);
        mapTriangles.insert(std::pair<std::vector<int>, std::vector<int> >(edge, idx));
        edge.clear();
        edge.push_back(second);
        edge.push_back(third);
        mapTriangles.insert(std::pair<std::vector<int>, std::vector<int> >(edge, idx));
        edge.clear();
        edge.push_back(first);
        edge.push_back(third);
        mapTriangles.insert(std::pair<std::vector<int>, std::vector<int> >(edge, idx));
    } else {
        std::vector<int> edge_1, edge_2;
        std::vector<int> idx;
        std::map<std::vector<int>, std::vector<int> >::iterator iter_1, iter_2;
        edge_1.push_back(first);
        edge_1.push_back(second);
        edge_2.push_back(second);
        edge_2.push_back(first);
        if((iter_1 = mapTriangles.find(edge_1)) == mapTriangles.end() &&
                (iter_2 = mapTriangles.find(edge_2)) == mapTriangles.end()) {
           idx.clear();
           idx.push_back(secondElem);
           mapTriangles.insert(std::pair<std::vector<int>, std::vector<int> >(edge_1, idx));
        } else {
           idx.clear();
           (iter_1 != mapTriangles.end()) ? idx = (*iter_1).second : idx = (*iter_2).second;
           idx.push_back(secondElem);
           (iter_1 != mapTriangles.end()) ? mapTriangles[edge_1] = idx : mapTriangles[edge_2] = idx;
        }

        edge_1.clear();
        edge_2.clear();
        edge_1.push_back(second);
        edge_1.push_back(third);
        edge_2.push_back(third);
        edge_2.push_back(second);
        if((iter_1 = mapTriangles.find(edge_1)) == mapTriangles.end() &&
                (iter_2 = mapTriangles.find(edge_2)) == mapTriangles.end()) {
           idx.clear();
           idx.push_back(secondElem);
           mapTriangles.insert(std::pair<std::vector<int>, std::vector<int> >(edge_1, idx));
        } else {
            idx.clear();
            (iter_1 != mapTriangles.end()) ? idx = (*iter_1).second : idx = (*iter_2).second;
            idx.push_back(secondElem);
            (iter_1 != mapTriangles.end()) ? mapTriangles[edge_1] = idx : mapTriangles[edge_2] = idx;
        }

        edge_1.clear();
        edge_2.clear();
        edge_1.push_back(first);
        edge_1.push_back(third);
        edge_2.push_back(third);
        edge_2.push_back(first);
        if((iter_1 = mapTriangles.find(edge_1)) == mapTriangles.end() &&
                (iter_2 = mapTriangles.find(edge_2)) == mapTriangles.end()) {
           idx.clear();
           idx.push_back(secondElem);
           mapTriangles.insert(std::pair<std::vector<int>, std::vector<int> >(edge_1, idx));
        } else {
            idx.clear();
            (iter_1 != mapTriangles.end()) ? idx = (*iter_1).second : idx = (*iter_2).second;
            idx.push_back(secondElem);
            (iter_1 != mapTriangles.end()) ? mapTriangles[edge_1] = idx : mapTriangles[edge_2] = idx;
        }


    }

}


void DelaunayTri::draw_subdiv() {
    std::vector<cv::Point> pt(3);
    for( size_t i = 0; i < triangleList.size(); i++ ) {
        cv::Vec6d t = triangleList[i];
        pt[0] = cv::Point(cvRound(t[0]) - 1, cvRound(t[1]) - 1);
        pt[1] = cv::Point(cvRound(t[2]) - 1, cvRound(t[3]) - 1);
        pt[2] = cv::Point(cvRound(t[4]) - 1, cvRound(t[5]) - 1);
        cv::line(img, pt[0], pt[1], cv::Scalar(255), 1, 8, 0);
        cv::line(img, pt[1], pt[2], cv::Scalar(255), 1, 8, 0);
        cv::line(img, pt[2], pt[0], cv::Scalar(255), 1, 8, 0);
    }
}


cv::Point2d DelaunayTri::GetCircumcenter(const cv::Point2d &a, const cv::Point2d &b, const cv::Point2d &c) {
    PointCGAL p = dt.circumcenter(PointCGAL(a.x, a.y), PointCGAL(b.x, b.y), PointCGAL(c.x, c.y));
    return cv::Point2d(p.x(), p.y());
}


void DelaunayTri::computeNeighbors() {
    int i, offset = 0;
    for(i = 0; i < size_; ++i) {
        int triangle[3];
        *triangle = indices[offset];
        *(triangle + 1) = indices[offset + 1];
        *(triangle + 2) = indices[offset + 2];
        std::vector<double> tempNeighbors;
        std::map<std::vector<int>, std::vector<int> >::iterator iter_1, iter_2;
        std::vector<int> edge_1, edge_2;
        edge_1.push_back(*triangle);
        edge_1.push_back(*(triangle + 1));
        edge_2.push_back(*(triangle + 1));
        edge_2.push_back(*triangle);

        bool contr;
        if( (iter_1 = mapTriangles.find(edge_1)) != mapTriangles.end()
                || (iter_2 = mapTriangles.find(edge_2)) != mapTriangles.end() ) {
            std::vector<int> n = (iter_1 != mapTriangles.end()) ? (*iter_1).second : (*iter_2).second;
            contr = false;
            size_t j;
            for(j = 0; j < n.size(); ++j) {
                if(n[j] != i) {
                    contr = true;
                    tempNeighbors.push_back(n[j] + 1);
                }
            }

            if(!contr)
                tempNeighbors.push_back(NAN);
        } else
            tempNeighbors.push_back(NAN);

        edge_1.clear();
        edge_2.clear();
        edge_1.push_back(*(triangle + 1));
        edge_1.push_back(*(triangle + 2));
        edge_2.push_back(*(triangle + 2));
        edge_2.push_back(*(triangle + 1));

        if( (iter_1 = mapTriangles.find(edge_1)) != mapTriangles.end()
                || (iter_2 = mapTriangles.find(edge_2)) != mapTriangles.end() ) {
            std::vector<int> n = (iter_1 != mapTriangles.end()) ? (*iter_1).second : (*iter_2).second;
            contr = false;
            size_t j;
            for(j = 0; j < n.size(); ++j) {
                if(n[j] != i) {
                    contr = true;
                    tempNeighbors.push_back(n[j] + 1);
                }
            }

            if(!contr)
                tempNeighbors.push_back(NAN);
        } else
            tempNeighbors.push_back(NAN);

        edge_1.clear();
        edge_2.clear();
        edge_1.push_back(*triangle);
        edge_1.push_back(*(triangle + 2));
        edge_2.push_back(*(triangle + 2));
        edge_2.push_back(*triangle);

        if( (iter_1 = mapTriangles.find(edge_1)) != mapTriangles.end()
                || (iter_2 = mapTriangles.find(edge_2)) != mapTriangles.end() ) {
            std::vector<int> n = (iter_1 != mapTriangles.end()) ? (*iter_1).second : (*iter_2).second;
            contr = false;
            size_t j;
            for(j = 0; j < n.size(); ++j) {
                if(n[j] != i) {
                    contr = true;
                    tempNeighbors.push_back(n[j] + 1);
                }
            }

            if(!contr)
                tempNeighbors.push_back(NAN);
        } else
            tempNeighbors.push_back(NAN);

        int j;
        for(j = 0; j < 3; j++) {
            neighbors[offset + j] = tempNeighbors[j];
        }
        offset += 3;
    }
}
