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

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <opencv2/opencv.hpp>

template<typename T>
cv::Mat vector2mat(const cv::Mat &_mat, const int &_rows, const int &_cols)
{
    cv::Mat v = cv::Mat(cv::Size(_cols, _rows), _mat.type());

    int x = 0, y = 0;
    for(int i = 0; i < _mat.rows; ++i)
    {
        v.at<T>(y, x) = _mat.at<T>(i, 0);
        y++;
        if(y == _rows)
        {
            y = 0;
            x++;
        /*    if(x == _cols) {
                break;
            }*/
        }
    }
    return v;
}

template<typename T>
cv::Mat mat2vector(const cv::Mat &_mat)
{
    int rows = _mat.rows;
    int cols = _mat.cols;
    int i = 0;

    cv::Mat v = cv::Mat(cv::Size(1, rows*cols), _mat.type());

    for(int x = 0; x < cols; ++x)
    {
        for(int y = 0; y < rows; ++y)
        {
            v.at<T>(i++) = _mat.at<T>(y, x);
        }
    }

    return v;
}


template<typename T>
T *allocateVector(const int &size, const bool &init, const T &initValue) {
    T* vec = new T[size];
    if(vec == NULL) {
        std::cerr << "Error to allocate Vector" << std::endl;
        exit (EXIT_FAILURE);
    }

    if(init) {
       int i;
       for(i = 0; i < size; ++i) {
           *(vec + i) = initValue;
       }
    }

    return vec;
}

template<typename T>
T *allocateMatrix(const int &rows, const int &cols,
                   const bool &init, const double &initValue) {
    int i;
    T *matrix;
    int size = rows*cols;
    matrix = new T[size];

    if(matrix == NULL) {
        std::cerr << "Error to allocate matrix" << std::endl;
        exit (EXIT_FAILURE);
    }

    if(init) {
        for(i = 0; i < size; ++i) {
            *(matrix + i ) = initValue;
        }
    }

    return matrix;
}

template<typename T>
T *copyMatrix(const T* _matrix, const int &rows,
               const int &cols) {

    T* matrix_ = allocateMatrix<T>(rows, cols, true, T(0));
    int i;
    int size = rows * cols;
    for(i = 0; i < size; ++i ) {
        *(matrix_ + i) = *(_matrix + i);
    }

    return matrix_;

}

template<typename T>
T *copyVector(const T* _vector, const int &size) {

    T* vector_ = allocateVector<T>(size, true, T());
    //T* vector_ = allocateMatrix<T>(rows, cols, true, T(0));
    int i;
    for(i = 0; i < size; ++i ) {
        *(vector_ + i) = *(_vector + i);
    }

    return vector_;

}

template<typename T>
bool deallocate(T *&elem) {
    if(elem != NULL) {
        delete[] elem;
        elem = NULL;
        if(elem != NULL)
            return false;
    }
    return true;
}

template<typename T>
inline T square(T p) {
    return (p * p);
}

template<typename T>
inline T sumSquare(T p1, T p2) {
    return ((p1 * p1) + (p2 * p2));
}

template<typename T>
inline T distPoint2Point(T p1_x, T p1_y, T p2_x, T p2_y) {
    return sqrt( square(p1_x - p2_x) +  square(p1_y - p2_y) );
}

template <class T>
class Point {
public:
    Point() {};

    Point(const T &x, const T &y) {
        this->x = x;
        this->y = y;
    }

    Point operator+(const Point& p) {
        Point result;
        result.x = p.x + this->x;
        result.y = p.y + this->y;
        return result;
    }

    Point operator-(const Point& p) {
        Point result;
        result.x = p.x - this->x;
        result.y = p.y - this->y;
        return result;
    }

    Point operator/(const Point& p) {
        Point result;
        result.x = p.x / this->x;
        result.y = p.y / this->y;
        return result;
    }

    Point operator*(const Point& p) {
        Point result;
        result.x = p.x * this->x;
        result.y = p.y * this->y;
        return result;
    }

    void operator=(const Point& p) {
        this->x = p.x;
        this->y = p.y;
    }

    T x,y;

};

#endif
