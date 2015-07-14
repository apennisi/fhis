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

#include "rgbSegmentation.h"
#include <fstream>
#include <sstream>
#include <time.h>


RgbSegmentation::RgbSegmentation(const cv::Mat &img, std::vector<cv::Point2i> XYPoints,
                                 const double &sigma_blurring, const double &sigma_c,
                                 const double &merge_threshold, cv::Mat &RgbLabels) {
	
    rows = img.rows;
    cols = img.cols;
    pixelNumber = rows * cols;

    int i, newSize = XYPoints.size() + 4, l, old_size = XYPoints.size();


    std::sort(XYPoints.begin(), XYPoints.end(), RgbSegmentation::compare_xy());
    std::vector<cv::Point2i> newXYPoints(newSize);
    newXYPoints[0] = cv::Point2i(1, 1);
    newXYPoints[1] = cv::Point2i(cols, 1);
    newXYPoints[2] = cv::Point2i(1, rows);
    newXYPoints[3] = cv::Point2i(cols, rows);

    std::copy(XYPoints.begin(), XYPoints.end(), newXYPoints.begin() + 4);

    //TRIANGULATION
    DelaunayTri d(newXYPoints, cols, rows, newSize);

    sizeTriangles = d.getSize();
    std::vector<double> radii = d.getRadii();
    std::vector<int> indices = d.getIndices();
    neighbors = d.getNeighbors();
    std::vector<cv::Point2d> centers = d.getCenters();

    //CIRCLE OVERLAPPING
    this->computeCircumcircleOverlap( radii, centers, neighbors );


    int width = 2*std::ceil(sigma_blurring) + 1;

    int total_width = width*2 + 1;

    gaussKernel.resize(total_width);

    double sum = .0;
    int x;
    double den = 1. / sigma_blurring;
    for(x = -width, i = 0; x < width && i < total_width; ++x, ++i)
    {
        gaussKernel[i]  = std::exp(-(square<double>(x * den)));
        sum += gaussKernel[i];
    }

    den = 1. / sum;
    
    for(i = 0; i < total_width; ++i)
    {
        gaussKernel[i] *= den;
    }


    std::vector<cv::Mat> bgr_planes;
    cv::split( img, bgr_planes );

    std::vector<double> bgr_dtMats(pixelNumber * 3, 0.);
    std::vector<double> bgr_gaussDtMats(pixelNumber * 3, 0.);

    size_t j;
    
    for(j = 0; j < 3; j++) {

        if(!im2single(bgr_planes[j], bgr_dtMats, j)) {
            std::cerr << "No B&W Image[" << j << "]" << std::endl;
            exit(-1);
        }

        imfilter( bgr_dtMats, gaussKernel, bgr_gaussDtMats, true, j );
        imfilter( bgr_gaussDtMats, gaussKernel, bgr_gaussDtMats, false, j );
    }

    std::vector<double> RGBb = dtMat2mat(bgr_gaussDtMats);


    this->mxCircleStats(radii, centers, RGBb);
    this->rgbtolab();
    this->computeSquaredEdgeDistances(sigma_c);

    MergeTriangles *mergeTriangles = new MergeTriangles(neighbors, color_weights, edge_lengths, sizeTriangles);
    mergeTriangles->getNewLabels(merge_costs, old_labels, new_labels, merge_lengths);
    delete mergeTriangles;

    std::vector<int> tri_labels = this->integrate_merges(old_labels, new_labels,
                                                merge_threshold);

    //LABELING
    RasterizeTriangles rt(rows, cols, newXYPoints, indices, sizeTriangles);
    std::vector<double> labels = rt.getRasterizedImage();
    labels2 = this->compactLabels(tri_labels, labels);

    RgbLabels = this->colorImageSegments(img, labels2, cv::Scalar(0, 0, 255));
}

RgbSegmentation::~RgbSegmentation() {

}

cv::Mat RgbSegmentation::colorImageSegments(const cv::Mat &bgr, const std::vector<double> &labels,
                                            const cv::Scalar &color) {

    int i, j, k;
    std::vector<cv::Mat> bgr_planes;
    cv::split(bgr.clone(), bgr_planes);

    std::vector<double> b(pixelNumber), g(pixelNumber), r(pixelNumber);
    std::vector<int> labelVec(pixelNumber);


    k = 0;
    double den = 1. / 255.;
    int offset = 0;
    for(j = 0; j < cols; ++j) {
        offset = 0;
        for(i = 0; i < rows; ++i) {
            b[k] = (double)bgr_planes[0].at<uchar>(i,j) * den;
            g[k] = (double)bgr_planes[1].at<uchar>(i,j) * den;
            r[k] = (double)bgr_planes[2].at<uchar>(i,j) * den;
            labelVec[k] = labels[offset + j];
            ++k;
            offset += cols;
        }
    }

    std::vector<double> red_sums = accumarray(labelVec, r);
    std::vector<double> green_sums = accumarray(labelVec, g);
    std::vector<double> blue_sums = accumarray(labelVec, b);

    int size = labelSize + 1;

    std::vector<double> lut(size * 3);


    lut[0] = (double)color[0] / 255.;
    lut[1] = (double)color[1] / 255.;
    lut[2] = (double)color[2] / 255.;

    offset = 3;

    for(i = 0; i < labelSize; ++i) {
        lut[offset] = blue_sums[i];
        lut[offset + 1] = green_sums[i];
        lut[offset + 2] = red_sums[i];
        offset += 3;
    }


    cv::Mat erodeLabels, dilateLabels;

    cv::Mat kernel = cv::Mat::ones(cv::Size(3, 3), CV_64F);

    cv::Mat labelsCopy = cv::Mat(cv::Size(cols, rows), CV_64FC1);

    offset = 0;
    
    for(i = 0; i < rows; ++i) {
        offset = i*cols;
        for(j = 0; j < cols; ++j) {
            labelsCopy.at<double>(i, j) = labels[offset + j];
        }
    }

    cv::erode(labelsCopy, erodeLabels, kernel);
    cv::dilate(labelsCopy, dilateLabels, kernel);

    std::vector<double> newLabels(pixelNumber);


    offset = 0;
    
    for(i = 0; i < rows; ++i) {
        offset = i*cols;
        for(j = 0; j < cols; ++j) {
            if((int)erodeLabels.at<double>(i,j) != (int)labels[offset + j] ||
                    (int)dilateLabels.at<double>(i,j) != (int)labels[offset + j]) {
                newLabels[offset + j] = 1;
            } else {
                newLabels[offset + j] = labels[i*cols + j] + 1;
            }
            newLabels[offset + j] = std::max(1., std::min(newLabels[offset + j],
                                                                (double) size));
        }
    }



    cv::Mat bgrNew = cv::Mat::zeros(bgr.size(), CV_8UC3);

    offset = 0;

    for(i = 0; i < rows; ++i)
    {
        for(j = 0; j < cols; ++j)
        {
            k = newLabels[offset + j] - 1;
            bgrNew.at<cv::Vec3b>(i,j)[0] = lut[k * 3] * 255;
            bgrNew.at<cv::Vec3b>(i,j)[1] = lut[k * 3 + 1] * 255;
            bgrNew.at<cv::Vec3b>(i,j)[2] = lut[k * 3 + 2] * 255;
        }
        offset += cols;
    }

    return bgrNew;

}

int RgbSegmentation::maxVal(const int *vec) {
    int maxV = -1000;
    int i, end = pixelNumber;
    for(i = 0; i < end; ++i) {
        if(*(vec + i) > maxV)
            maxV = *(vec + i);
    }
    return maxV;
}

int RgbSegmentation::findVal(const int *vec, const int &elem, const int &start) {
    int idx = -1000;
    int i;
    for( i = start; i < pixelNumber; ++i) {
        if(elem == *(vec + i))
            return i;
    }

    return idx;
}

std::vector<double> RgbSegmentation::accumarray(const std::vector<int> &subs, const std::vector<double> &vals) {
    int i;

    labelSize = *std::max_element(subs.begin(), subs.end());
            //maxVal(subs);

    std::vector<AccumElem> elements(labelSize);

    int index = 0;

    for(i = 0; i < pixelNumber; ++i)
    {
        index = (subs[i] == 0) ? subs[i] : subs[i] - 1;
        elements[index].index++;
        elements[index].sum += vals[i];
    }

    std::vector<double> out(labelSize);
    for(i = 0; i < labelSize; ++i) {
        out[i] = elements[i].sum / (double)elements[i].index;
    }
    return out;
}

std::vector<double> RgbSegmentation::compactLabels(const std::vector<int> &tri_labels,
                                       const std::vector<double> &labels) {

    std::vector<double> newLabels = labels;

    int i, j, k;
    int sizeLabels = -1;
    
    for(i = 0; i < sizeTriangles; ++i)
    {
        if(tri_labels[i] > sizeLabels)
        {
            sizeLabels = tri_labels[i];
        }
    }

    bool checked[sizeLabels + 1];
    std::fill_n(checked, sizeLabels + 1, 0);
    int labelArrayC[sizeLabels + 1];
    std::fill_n(labelArrayC, sizeLabels + 1, -1);

    int labelArrayA[pixelNumber];

    k = 0;

    int offset;

    for(j = 0; j < cols; ++j) {
        offset = 0;
        for(i = 0; i < rows; ++i) {
            int idx = labels[offset + j] - 1;
            int elem = tri_labels[idx];
            if(!checked[elem]) {
                checked[elem] = true;
            }
            labelArrayA[k] = elem;
            ++k;
            offset += cols;
        }
    }

    int counter = 0;
    for(i = 0; i < sizeLabels; ++i) {
        if(checked[i]) {
            *(labelArrayC + i) = counter;
            counter++;
        }
    }

    int labelArrayIc[pixelNumber];

    for(i = 0; i < pixelNumber; i++) {
        *(labelArrayIc + i) = *(labelArrayC + *(labelArrayA + i) );
    }

    i = 0, j = 0;
    offset = 0;
    for(k = 0; k < pixelNumber && j < cols; k++) {
        newLabels[offset + j] = *(labelArrayIc + k) + 1;
        ++i;
        offset += cols;
        if(i == rows) {
            i = 0;
            offset = 0;
            ++j;
        }
    }

    return newLabels;
}

std::vector<int> RgbSegmentation::integrate_merges(const std::vector<int> &old_labels,
                                       const std::vector<int> &new_labels,
                                       const double &merge_threshold) {

    std::vector<int> lut(sizeTriangles), tail(sizeTriangles), next(sizeTriangles, 0);

    
    for(int i = 0; i < sizeTriangles; ++i) {
        lut[i] = tail[i] = i + 1;
    }


    int size = sizeTriangles - 1;
    int old_label = -1;
    int new_label = -1;

    for(int i = 0; i < size; ++i)
    {
        if(merge_costs[i] < merge_threshold)
        {
            old_label = old_labels[i];
            new_label = new_labels[i];
            next[tail[new_label - 1] - 1] = old_label;
            tail[new_label - 1] = tail[old_label - 1];
            tail[old_label - 1] = 0;
        }
    }


    for(int i = 0; i < sizeTriangles; ++i) {
        if(tail[i]) {
            int idx = i + 1;
            while(idx) {
                lut[idx - 1] = i + 1;
                idx = next[idx - 1];
            }
        }
    }
    return lut;
}

std::vector<double> RgbSegmentation::dtMat2mat(const std::vector<double> &src) {
    std::vector<double> bgr(pixelNumber * 3, 0.);
    int i, j;
    int doublePixelNumber = pixelNumber*2;

    int offset = 0;
    
    for(i = 0; i < rows; ++i) {
        offset = i*cols;
        for(j = 0; j < cols; ++j) {
            bgr[offset + j] = src[offset + j]*255.;
            bgr[offset + pixelNumber + j] = src[offset + pixelNumber + j]*255.;
            bgr[offset + doublePixelNumber + j] = src[offset + doublePixelNumber + j]*255.;
        }
    }

    return bgr;
}

bool RgbSegmentation::imfilter(const std::vector<double> &src, const std::vector<double> &filter,
                             std::vector<double> &_mat, const bool &isTranspose, const int &plane) {

    int factor = pixelNumber*plane;
    int i, j, k;
    std::vector<double> vec;
    double sum;

    if( isTranspose ) {
        vec.resize(cols + 15, 0.0);
        for( i=0; i<rows; ++i ) {
            for(j=0; j<cols; ++j )
                vec[j + 7] = src[i*cols + factor + j];
            for( j=0; j<7; ++j )
                vec[j] = src[i*cols + factor];
            for( j=0; j<7; ++j )
                vec[j + cols + 7] = src[i*cols + factor + cols - 1];

            for( j=0; j<cols; ++j ) {
                sum = 0.0;
                for( k=0; k<15; ++k ) {
                    sum += (vec[j + k]) * (filter[14 - k]);
                }
                _mat[i*cols + factor + j] = sum;
            }
        }
    } else {
        vec.resize(rows+15, 0.0);
        for( i=0; i < cols; ++i ) {
            for( j=0; j<rows; ++j )
                vec[j + 7] = src[j*cols + factor + i];
            for( j=0; j<7; ++j )
                vec[j] = src[factor + i];
            for( j=0; j<7; ++j )
                vec[j + rows + 7] = src[(rows - 1)*cols + factor + i];
            for( j=0; j<rows; ++j ) {
                sum = 0.0;
                for( int k=0; k<15; ++k ) {
                    sum += vec[j + k] * filter[14 - k];
                }
                _mat[j*cols + factor + i] = sum;
            }
        }
    }

    return true;
}

bool RgbSegmentation::im2single( const cv::Mat &src, std::vector<double> &mat, const int &plane ) {

    if( src.channels() != 1 ) {
        return false;
    }

    unsigned char *img = src.data;

    double _eps = 1.0e-10;

    double max_range = 1.0;
    double half_range = 0.0;
    if( src.depth() == CV_8U ) {
        max_range = 255.;
    } else if ( src.depth() == CV_8S ) {
        max_range = 255.;
        half_range = 128.;
    } else if( src.depth() == CV_16U ) {
        max_range = 65535.;
    } else if (src.depth() == CV_16S ) {
        max_range = 65535.;
        half_range = 32768.;
    }

    int factor = pixelNumber*plane;
    double denominator = 1. / (max_range + _eps);

    int offset = 0;
    
    for( int i=0; i<rows; ++i) {
        offset = i*cols;
        for( int j=0; j<cols; ++j) {
            mat[offset + factor + j] = (( *(img + offset + j) + half_range )
                                            * denominator);
        }
        //offset += cols;
    }

    return true;
}

void RgbSegmentation::rgbtolab() {

    //int nsize = mean_color_rgb.size();
    int i, j;
    mean_color_lab.resize(sizeTriangles * 3, 0.);
    const int size = 3;
    int offset = 0;
    
    for(i = 0; i < sizeTriangles; ++i) {

        double rgb[3], XYZ[3];
        for(j = 0; j < size; ++j) {
            double tmp_value = mean_color_rgb[offset + j] / 255.;

            (tmp_value > 0.0404482362771076) ? tmp_value = std::pow(((tmp_value + 0.055) / 1.055), 2.4) :
                    tmp_value /= 12.92;

            *(rgb + j) = tmp_value;

        }

        double X, Y, Z;

        X = *rgb * 0.436052025 + *(rgb + 1) * 0.385081593 + *(rgb + 2) * 0.143087414;
        Y = *rgb * 0.222491598 + *(rgb + 1) * 0.716886060 + *(rgb + 2) * 0.060621486;
        Z = *rgb * 0.013929122 + *(rgb + 1) * 0.097097002 + *(rgb + 2) * 0.714185470;


        *XYZ =  X / 0.964296 ;
        *(XYZ + 1) = Y;
        *(XYZ + 2) = Z / 0.825106 ;


        for(j = 0; j < size; ++j) {

            double tmp_value = *(XYZ + j);

            tmp_value = (tmp_value > 0.008856452) ?  std::pow(tmp_value, 0.333333333333333) :
                     (tmp_value * 7.787037037) + (0.137931034);

            *(XYZ + j) = tmp_value;

        }


        mean_color_lab[offset] = (116. * *(XYZ + 1)) - 16.;
        mean_color_lab[offset + 1] = 500. * (*XYZ - *(XYZ + 1));
        mean_color_lab[offset + 2] = 200. * (*(XYZ + 1) - *(XYZ + 2));

        offset += size;

    }
}

void RgbSegmentation::computeSquaredEdgeDistances(const double &sigma_c) {

    int size = 3, k, i, j;


    std::vector<double> dc2(sizeTriangles * size, 0.);
    color_weights.resize(sizeTriangles * size);
    for(k = 0; k < size; k++) {
        for(i = 0; i < sizeTriangles; i++) {
            for(j = 0; j < size; j++) {
                if(!std::isnan(neighbors[i*size + j])) {
                    int currentIdx = neighbors[i*size + j] - 1;
                    double delta = mean_color_lab[currentIdx*size + k] -
                            mean_color_lab[i*size + k];
                    dc2[i*size + j] += square<double>(delta);
                } else {
                    dc2[i*size + j] = 0.;
                }
            }

        }
    }

    double elem = 1./(sigma_c*sigma_c);

    int offset = 0;
    
    for(i = 0; i < sizeTriangles; i++) {
        offset = i*size;
        for(j = 0; j < size; j++) {
            color_weights[offset + j] = 1. - (std::exp(-elem * dc2[offset + j]));
        }
    }
}

void RgbSegmentation::mxCircleStats(const std::vector<double> &radii, const std::vector<cv::Point2d> &centers,
                                    const std::vector<double> &RGBb) {

    const int ndims = 3; //CHANGE FOR BLACK AND WHITE

    mean_color_rgb.resize(sizeTriangles * 3, 0.);

    int min_col, max_col, min_row, max_row, count, row, col, j, k;
    double radius, dx, dy;
    cv::Point2d current_center;

    int offset = 0;
    for(int i = 0; i < sizeTriangles; i++) {

        current_center.x = centers[i].x - 1;
        current_center.y = centers[i].y - 1;

        radius = radii[i];

        min_col = (int)(current_center.x - radius);
        if (min_col < 0) min_col = 0;

        max_col = (int)(current_center.x + radius);
        if (max_col >= cols) max_col = cols - 1;

        min_row = (int)(current_center.y - radius);
        if (min_row < 0) min_row = 0;

        max_row = (int)(current_center.y + radius);
        if (max_row >= rows) max_row = rows - 1;

        radius = square<double>(radius) + 0.25;

        count = 0;

        for (row = min_row; row <= max_row; ++row) {
            for (col = min_col; col <= max_col; ++col) {
                dx = col - current_center.x;
                dy = row - current_center.y;

                if ( (dx*dx + dy*dy) < radius ) {

                    for (j=0, k = 2; j < ndims && k >= 0; ++j, k--) {
                        mean_color_rgb[i*ndims + j] += RGBb[row*cols + k*pixelNumber + col];
                    }

                    ++count;
                }
            }
        }

        for (j=0; j < ndims; ++j) mean_color_rgb[offset + j] /= count;
        offset += 3;
    }
}

void RgbSegmentation::quicksort(std::vector<cv::Point2i> &XYpoints, int p, int r) {
    if ( p < r ) {
        int j = partition(XYpoints, p, r);
        quicksort(XYpoints, p, j-1);
        quicksort(XYpoints, j+1, r);
    }
}

int RgbSegmentation::partition(std::vector<cv::Point2i> &XYpoints, int p, int r) {

    int pivot = XYpoints[r].x;

    while ( p < r ) {

        while ( XYpoints[p].x  < pivot ) {
            p++;
        }

        while (XYpoints[r].x > pivot ) {
            r--;
        }

        if (XYpoints[p].x == XYpoints[r].x )
            p++;
        else if ( p < r ) {
            cv::Point2i temp = XYpoints[p];
            XYpoints[p] = XYpoints[r];
            XYpoints[r] = temp;
        }
    }

    return r;
}


void RgbSegmentation::computeCircumcircleOverlap(const std::vector<double> &radii, const std::vector<cv::Point2d> &centers,
                                                 const std::vector<double> &neighbors) {


       edge_lengths.resize(sizeTriangles * 3, 0.);

       for(int i = 0; i < 3; i++) {

           double radii_;
           double r, R, c1, c2;
           double distance;
           double theta1, theta2, s1, s2;

           int offset = 0;
           for(int j = 0; j < sizeTriangles; ++j) {

               int idx = std::isnan(neighbors[offset + i]) ? 0 : neighbors[offset + i] - 1;

               radii_ = radii[idx];

               r = std::min(radii[j], radii_);
               R = std::max(radii[j], radii_);


               cv::Point2d p1 = centers[j];
               cv::Point2d p2 = centers[idx];

               distance = distPoint2Point<double>(p1.x, p1.y, p2.x, p2.y);


               if(distance == 0.) {
                   edge_lengths[offset + i] = 2 * (double)R;

               } else {
                   double squareR = square<double>(R);
                   double squareDistance = square<double>(distance);
                   double squarer = square<double>(r);
                   double denominator = 1. / (2 * R * distance);
                   c1 = ((squareR + squareDistance) -
                                 squarer ) * denominator ;
                   c2 = ( (squarer + squareDistance) -
                                 squareR) * squareDistance;

                   c1 = std::max((double)-1., std::min((double)c1, (double)1.));
                   c2 = std::max((double)-1., std::min((double)c2, (double)1.));

                   theta1 = std::acos(c1);
                   theta2 = std::acos(c2);

                   s1 = std::sqrt(1 - square<double>(c1));
                   s2 = std::sqrt(1 - square<double>(c2));

                   edge_lengths[offset + i] = 2*(s1*R);
               }

               offset += 3;
           }
       }




}
