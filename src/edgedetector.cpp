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

#include "edgedetector.h"


EdgeDetector::EdgeDetector(const cv::Mat &src, cv::Mat &edgeImg,
                           const double &edge_thresh, const double &edge_sigma)
{
    if( !src.data ) {
        std::cerr << "Empty Image" << std::endl;
        exit(-1);
    }

    rows = src.rows;
    cols = src.cols;
    numPixels = rows * cols;

    computGradient(edge_sigma);

    //MAGIC NUMBERS
    _eps = 1.0e-10;
    PercentOfPixelsNotEdges = .7;
    ThresholdRatio = .4;

    std::vector<double> mat(numPixels, 0.);

    if( !im2single(src, mat) ) {
        std::cerr << "No B&W Image" << std::endl;
        exit(-1);
    }


    std::vector<double> dx(numPixels, 0.), dy(numPixels, 0.);

    if( !smoothGradient( mat, dx, dy ) ) {
        std::cout << "No mat!" << std::endl;
        exit(-1);
    }

    std::vector<double> mag(numPixels, 0.);
    hypot( dx, dy, mag );

    double _low=0.0, _high=0.0;
    selectThresholds( _low, _high, edge_thresh );

    std::vector<uchar> dst;
    thinAndThreshold( mag, dx, dy, _low, _high, dst );

    dtMat2Img( dst, edgeImg );

}

EdgeDetector::~EdgeDetector() {}

void EdgeDetector::computGradient(const double &edge_sigma) {

    dimKernel = 8*std::ceil(edge_sigma);

    double n = ((double)dimKernel - 1) / (double)2;

    gaussKernel.resize(dimKernel);
    derivGaussKernel.resize(dimKernel);

    double startKernel = -1 * n;

    double c = 1 / (std::sqrt(2 * CV_PI) * edge_sigma);
    double sum = 0.0;
    double den = 1. / double(2 * edge_sigma * edge_sigma);

    for(int i = 0; i < dimKernel; i++)
    {
        gaussKernel[i] = c * std::exp(-(startKernel*startKernel) * den);
        sum += gaussKernel[i];
        startKernel += 1.0;
    }

    den = 1. / sum;

    for(int i = 0; i < dimKernel; i++)
    {
        gaussKernel[i] *= den;
    }

    gradient();

    std::vector<int> posVals, negVals;
    double sumPosVals = .0, sumNegVals = .0;

    for(int i = 0; i < dimKernel; i++)
    {
        if(derivGaussKernel[i] > 0)
        {
            posVals.push_back(i);
            sumPosVals += derivGaussKernel[i];
        } else if(derivGaussKernel[i] < 0)
        {
            negVals.push_back(i);
            sumNegVals += derivGaussKernel[i];
        }
    }

    size_t i;
    den = 1. / sumPosVals;

    for(i = 0; i < posVals.size(); i++)
    {
        derivGaussKernel[posVals[i]] *= den;
    }

    den = 1. / std::abs(sumNegVals);
    for(i = 0; i < negVals.size(); i++)
    {
        derivGaussKernel[negVals[i]] *= den;
    }
}

void EdgeDetector::gradient() {
    derivGaussKernel[0] = gaussKernel[1] - gaussKernel[0];

    int end = dimKernel-1;

    double den = 1. / 2.;

    for(int i = 1; i < end; i++)
    {
        derivGaussKernel[i] = (gaussKernel[i+1] - gaussKernel[i-1]) * den;
    }

    derivGaussKernel[dimKernel - 1] = gaussKernel[dimKernel - 1] - gaussKernel[dimKernel - 2];
}

bool EdgeDetector::im2single( const cv::Mat &src, std::vector<double> &mat ) {

    if( src.channels() != 1 ) {
        return false;
    }

    int i, j;

    double max_range = 1.0;
    double half_range = 0.0;
    if( src.depth() == CV_8U ) {
        max_range = 255.0;
    } else if ( src.depth() == CV_8S ) {
        max_range = 255.0;
        half_range = 128.0;
    } else if( src.depth() == CV_16U ) {
        max_range = 65535.0;
    } else if (src.depth() == CV_16S ) {
        max_range = 65535.0;
        half_range = 32768.0;
    }


    double denominator = 1. / double(max_range + _eps);
    int offset = 0;

    for( i=0; i< rows; ++i)
    {
        offset = i * cols;
        for( j=0; j < cols; ++j)
        {
            mat[offset + j] = (double)((uchar) src.at<uchar>(i,j));
            mat[offset + j] = ( mat[offset + j] + half_range ) * double(denominator);
        }
    }

    return true;
}

bool EdgeDetector::imfilter(const std::vector<double> &src, const std::vector<double> &filter,
                             std::vector<double> &_mat, const bool &isTranspose) {

    std::vector<double> vec;
    int i, j, k;
    double sum;
    if( isTranspose ) {
        vec.resize(cols+16, 0.0);
        for(i=0; i < rows; ++i ) {
            for( j=0; j<cols; ++j )
                vec[j + 8] = src[i*cols + j];
            for( j=0; j<8; ++j )
                vec[j] = src[i*cols];
            int end = cols-1;
            for( j=0; j<8; ++j )
                vec[j + cols + 8] = src[i*cols + end];
            for( int j=0; j < cols; ++j )
            {
                sum = 0.0;
                for( k=0; k<16; ++k ) {
                    sum += vec[j + k + 1] * filter[15 - k];
                }
                _mat[i*cols + j] = sum;
            }
        }
    } else {
        vec.resize(rows+16, 0.0);
        for( i = 0; i < cols; ++i ) {
            for( j = 0; j < rows; ++j )
                vec[j + 8] = src[j*cols + i];
            for( j=0; j<8; ++j )
                vec[j] = src[i];
            int init = (rows - 1)*cols;
            for( j=0; j<8; ++j )
                vec[j + rows + 8] = src[init + i];
            for( j=0; j<rows; ++j ) {
                sum = 0.0;
                for( k=0; k<16; ++k ) {
                    sum += vec[j + k + 1] * filter[15 - k];
                }
                _mat[j*cols + i] = sum;
            }
        }
    }


    return true;
}

bool EdgeDetector::smoothGradient(const std::vector<double> &src, std::vector<double> &dx,
                                   std::vector<double> &dy) {

    imfilter( src, gaussKernel, dx, false );
    imfilter( dx, derivGaussKernel, dx, true );
    imfilter( src, gaussKernel, dy, true );
    imfilter( dy, derivGaussKernel, dy, false );

    return true;
}

void EdgeDetector::hypot(const std::vector<double> &dx, const std::vector<double> &dy,
                         std::vector<double> &magGrad) {

    int i, j;
    double maxValue = 0.0;
    double p1, p2, t;
    int offset = 0;

    for( i = 0; i < rows; ++i )
    {
        offset = i * cols;
        for( j = 0; j < cols; ++j ) {
            p1 = dx[offset + j];
            p2 = dy[offset + j];
            t = std::sqrt( sumSquare<double>(p1, p2) );
            if( maxValue < t ) maxValue = t;
            magGrad[offset + j] = t;
        }
        //offset += cols;
    }

    if( maxValue > _eps ) {
        offset = 0;
        double den = 1. / maxValue;
        for( i = 0; i < rows; i++ ) {
            offset = i * cols;
            for( j = 0; j < cols; j++ ) {
                magGrad[offset + j] *= den;
            }
            //offset += cols;
        }
    }
}

void EdgeDetector::selectThresholds(double& lowThresh, double& highThresh, const double &thresh ) {

    highThresh = thresh;

    lowThresh = ThresholdRatio * highThresh;
}

bool EdgeDetector::cannyFindLocalMaxima(const std::vector<double> &src, const std::vector<double> &dx,
                                         const std::vector<double> &dy, const cv::Point2i &pos) {

    int ix = pos.x;
    int iy = pos.y;

    if( ix <= 0 || ix >= rows-1 || iy <= 0 || iy >= cols-1 ) {
        return false;
    }
    int ix_cols = ix*cols;
    double dx_val = dx[ix_cols + iy];
    double dy_val = dy[ix_cols + iy];
    double gradmag1, gradmag2, gradmag = src[ix_cols + iy];
    bool case1 = false, case2 = false, case3 = false, case4 = false;

    if( (dy_val <= 0 && dx_val + dy_val > 0) || (dy_val >= 0 && dx_val + dy_val < 0) ) {
        double dam = dx_val >= 0 ? ( dx_val + _eps ) : ( dx_val - _eps );
        double d = fabs( dy_val / dam );
        gradmag1 = src[ix_cols + iy+1] * (1-d) + src[ix_cols - cols + iy + 1] * d;
        gradmag2 = src[ix_cols + iy-1] * (1-d) + src[ix_cols + cols + iy - 1] * d;
        case1 = (gradmag >= gradmag1 && gradmag >= gradmag2);
    }

    if( (dx_val>0&&dx_val+dy_val<=0) || (dx_val<0&&dx_val+dy_val>=0) ) {

        double dam = dy_val >= 0 ? ( dy_val + _eps ) : ( dy_val - _eps );
        double d = fabs( dx_val / dam );
        gradmag1 = src[ix_cols - cols + iy] * (1-d) + src[ix_cols - cols + iy + 1] * d;
        gradmag2 = src[ix_cols + cols + iy] * (1-d) + src[ix_cols + cols + iy - 1] * d;
        case2 = (gradmag >= gradmag1 && gradmag >= gradmag2);

    }

    if( (dx_val<=0&&dx_val>dy_val) || (dx_val>=0&&dx_val<dy_val) ) {

        double dam = dy_val >= 0 ? ( dy_val + _eps ) : ( dy_val - _eps );
        double d = fabs( dx_val / dam );
        gradmag1 = src[ix_cols - cols + iy] * (1-d) + src[ix_cols - cols + iy - 1] * d;
        gradmag2 = src[ix_cols + cols + iy] * (1-d) + src[ix_cols + cols + iy + 1] * d;
        case3 = (gradmag >= gradmag1 && gradmag >= gradmag2);

    }

    if( (dy_val<0&&dx_val<=dy_val) || (dy_val>0&&dx_val>=dy_val ) ) {

        double dam = dx_val >= 0 ? ( dx_val + _eps ) : ( dx_val - _eps );
        double d = fabs( dy_val / dam );
        gradmag1 = src[ix_cols + iy - 1] * (1-d) + src[ix_cols - cols + iy - 1] * d;
        gradmag2 = src[ix_cols + iy + 1] * (1-d) + src[ix_cols + cols + iy + 1] * d;
        case4 = (gradmag >= gradmag1 && gradmag >= gradmag2);

    }

    return ( case1 || case2 || case3 || case4 );

}

void EdgeDetector::thinAndThreshold( const std::vector<double> &src, const std::vector<double> &dx, const std::vector<double> &dy,
                                     const double& lowThresh, const double& highThresh, std::vector<uchar> &dst) {

    std::vector<uchar> idxWeak(numPixels, 0);
    std::vector<uchar> idxWeakNeg(numPixels, 1);

    //Position pp;
    //Position *idxStrongPts = allocateVector<Position>(std::ceil(rows * cols * PercentOfPixelsNotEdges), false, pp);
    std::vector<cv::Point2i> idxStrongPts(std::ceil(numPixels*PercentOfPixelsNotEdges));

    int i, j, k = 0;
    int offset = 0;

    for( i = 0; i < rows; ++i ) {
        for( j = 0; j < cols; ++j ) {
            cv::Point2i p;
            p.x = i;
            p.y = j;
            if( cannyFindLocalMaxima( src, dx, dy,  p) ) {
                if ( src[offset + j] > lowThresh ) {
                    idxWeak[offset + j] = 1;
                    idxWeakNeg[offset + j] = 0;
                }

                if( src[offset + j] > highThresh ) {
                    idxStrongPts[k] = p;
                    ++k;
                }
            }
        }
        offset += cols;
    }

    bwselect( idxWeak, idxWeakNeg, idxStrongPts, k );
    dst = idxWeak;
}

void EdgeDetector::bwselect(std::vector<uchar> &idxWeak, std::vector<uchar> &idxWeakNeg,
                            const std::vector<cv::Point2i> &idxStrongPts, const int &size) {

    const int ptSize = size;
    const int d_x[8] = {0, 0,1,1, 1,-1,-1, -1};
    const int d_y[8] = {1,-1,0,1,-1, 1, 0, -1};

    for( int i=0; i<ptSize; ++i ) {
        int ix = idxStrongPts[i].x;
        int iy = idxStrongPts[i].y;
        int ix_cols = ix*cols;
        if( idxWeakNeg[ix_cols + iy] == 0 ) {
            std::queue<cv::Point2i> Q;
            Q.push( idxStrongPts[i] );
            while( !Q.empty() ) {
                cv::Point2i pos  = Q.front();
                Q.pop();
                ix = pos.x;
                iy = pos.y;
                idxWeakNeg[ix_cols + iy] = 1;
                for( int j=0; j<8; ++j ) {
                    if( 0 == idxWeakNeg[(ix + *(d_x + j))*cols + iy + *(d_y + j) ] ) {
                        cv::Point2i pos;
                        pos.x = ix + *(d_x + j);
                        pos.y = iy + *(d_y + j);
                        Q.push( pos );
                        idxWeakNeg[(ix + *(d_x + j))*cols + iy + *(d_y + j) ] = 1;
                    }
                }
            }
        }
    }

    int offset = 0;

    for( int i=0; i< numPixels; ++i )
    {
        idxWeak[i] &= idxWeakNeg[i];
    }
}

void EdgeDetector::dtMat2Img( const std::vector<uchar> &mat, cv::Mat &img ) {

    img = cv::Mat( cv::Size(cols, rows), CV_8UC1 );

    int offset = 0;
    for( int i=0; i< rows; i++ ) {
        for( int j=0; j< cols; j++ ) {
            img.at<uchar>(i,j) = (mat[offset + j] == 0 ? 0 : 255 );
            if(mat[offset + j] != 0)
                points.push_back(cv::Point2i(j,i));
        }
        offset += cols;
    }
}
