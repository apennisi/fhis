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

//C
#include <stdio.h>
//C++
#include <iostream>
#include <iomanip>
//OpenCV
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>

//fhis
#include "rgbSegmentation.h"
#include "edgedetector.h"

static const double thresh = 0.04;
static const double sigma = 2.;
static const double sigma_blurring = 5.;
static const double sigma_c = 4.f;
static const double merge_threshold = 0.6;

cv::Mat frame;     //current frame
double fps;   //frame per second for the input video sequence
int keyboard;  //input from keyboard

/**
 * Function Headers
*/
void help();
void processVideo(char* videoFilename);
void processImages(char* firstFrameFilename);
void processRtsp(char* address);
void processSingleImg(char* image);


void help()
{
    std::cout
            << "--------------------------------------------------------------------------" << std::endl
            << "FHIS Image Segmentation Library "                                           << std::endl
            << "This file main.cpp contains an example of usage for"                        << std::endl
            << "FHIS algorithm described in"                                                << std::endl
            << "Camillo Taylor AND Anthony Cowley"                                          << std::endl
            << "\"Parsing Indoor Scenes Using RGB-D Imagery\""                              << std::endl
            << "In Proceedings of Robotics: Science and Systems"			    << std::endl
            << "July, 2012."                                                                << std::endl
            << std::endl
            << "written by Andrea Pennisi"                                                  << std::endl
            << "andrea.pennisi@gmail.com"                                                   << std::endl
            << "--------------------------------------------------------------------------" << std::endl
            << "You can process both videos (-vid), images (-img), "
            << "and rtsp streams (-rtsp)."                                                  << std::endl
            << std::endl
            << "Usage:"                                                                     << std::endl
            << "fhis {-vid <video filename>|-img <image filename [-fps <value>]|"
            << "-rtsp <RTSP address>}"                                                      << std::endl
            << "for example: fhis -vid video.avi"                                           << std::endl
            << "or: fhis -singleimg /data/images/1.png"                                     << std::endl
            << "or: fhis -img /data/images/1.png"                                           << std::endl
            << "or: fhis -img /data/images/1.png -fps 7"                                    << std::endl
            << "or: fhis -rtsp rtsp://example.com/video.mp4"                                << std::endl
            << "--------------------------------------------------------------------------" << std::endl
            << std::endl;
}

int main(int argc, char* argv[]) {

    //print help information
    help();

    //check for the input parameter correctness
    if(argc < 3) {
        std::cerr <<"Incorrect input list" << std::endl;
        std::cerr <<"exiting..." << std::endl;
        return EXIT_FAILURE;
    }

    cv::namedWindow("Frame");

    if(strcmp(argv[1], "-vid") == 0) {

        //input data coming from a video
        processVideo(argv[2]);
    }
    else if(strcmp(argv[1], "-singleimg") == 0) {
        processSingleImg(argv[2]);
    }
    else if(strcmp(argv[1], "-img") == 0) {
        //input data coming from a sequence of images
        if(argc > 3 && strcmp(argv[3], "-fps") == 0) {
            fps = atof(argv[4]);
        }
        else {
            fps = 25;
        }

        processImages(argv[2]);
    }
    if(strcmp(argv[1], "-rtsp") == 0) {
        //input data coming from an rtsp stream
        processRtsp(argv[2]);
    }
    else {
        //error in reading input parameters
        std::cerr <<"Please, check the input parameters." << std::endl;
        std::cerr <<"Exiting..." << std::endl;
        return EXIT_FAILURE;
    }
    //destroy GUI windows
    cv::destroyAllWindows();
    return EXIT_SUCCESS;
}

/**
* @function processImages
* WARNING: this function can read only image sequences in the form
* <n>.<ext>
* where <n> is a number without zeros as prefix, e.g., 7
* and <ext> is the extension provided by command line, e.g., png, jpg, etc.
* Example of sequence:
*             1.png, 2.png, ..., 15.png, 16.png, ..., 128.png, 129.png, ...
*/
void processImages(char* fistFrameFilename) {
    //read the first file of the sequence
    std::cout << "qui" << std::endl;
    frame = cv::imread(fistFrameFilename);
    if(!frame.data){
        //error in opening the first image
        std::cerr << "Unable to open first image frame: " << fistFrameFilename << std::endl;
        exit(EXIT_FAILURE);
    }

    //current image filename
    std::string fn(fistFrameFilename);
    int frameNumber = 0;

    //read input data. ESC or 'q' for quitting
    while( (char)keyboard != 'q' && (char)keyboard != 27 ){

        //Edge Detection
        cv::Mat grayImg;
        cv::cvtColor(frame.clone(), grayImg, CV_BGR2GRAY);

        cv::Mat edge;
        EdgeDetector edges(grayImg, edge, thresh, sigma);

        cv::imshow("Edge Detection", edge);

        std::vector<cv::Point2i> XYpoints = edges.getEdgesPoints();

        //RGB Segmentation
        cv::Mat rgbSegmentation;

        RgbSegmentation p(frame, XYpoints, sigma_blurring, sigma_c, merge_threshold, rgbSegmentation);

        cv::imshow("Rgb Segmented", rgbSegmentation);

        //get the frame number and write it on the current frame
        size_t index = fn.find_last_of("/");
        if(index == std::string::npos) {
            index = fn.find_last_of("\\");
        }
        size_t index2 = fn.find_last_of(".");
        std::string prefix = fn.substr(0,index+1);
        std::string suffix = fn.substr(index2);
        std::string frameNumberString = fn.substr(index+1, index2-index-1);

        cv::rectangle(frame, cv::Point(10, 2), cv::Point(100,20),
                  cv::Scalar(255,255,255), -1);
        cv::putText(frame, frameNumberString.c_str(), cv::Point(15, 15),
                cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));
        //show the current frame and the fg masks
        cv::imshow("Frame", frame);
        
        //get the input from the keyboard
        keyboard = cv::waitKey( 10 );

        //search for the next image in the sequence
        std::ostringstream oss;
        oss << (frameNumber++);
        std::string nextFrameNumberString = oss.str();

        std::stringstream nextFrameStream;

        nextFrameStream << std::setw(4) << std::setfill('0') << nextFrameNumberString;

        std::string nextFrameFilename = prefix + "frame_" + nextFrameStream.str() + suffix;


        //read the next frame
        frame = cv::imread(nextFrameFilename);
        if(!frame.data){
            //error in opening the next image in the sequence
            std::cerr << "Unable to open image frame: " << nextFrameFilename << std::endl;
            exit(EXIT_FAILURE);
        }
        //update the path of the current frame
        fn.assign(nextFrameFilename);
    }
}

void processRtsp(char* address) {
    //create the capture object
    cv::VideoCapture capture(address);
    if(!capture.isOpened()){
        //error in opening the video input
        std::cerr << "Unable to open rtsp stream: " << address << std::endl;
        exit(EXIT_FAILURE);
    }

    //read input data. ESC or 'q' for quitting
    while( (char)keyboard != 'q' && (char)keyboard != 27 ){
        //read the current frame
        if(!capture.read(frame)) {
            std::cerr << "Unable to read next frame." << std::endl;
            std::cerr << "Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
        //Edge Detection
        cv::Mat grayImg;
        cv::cvtColor(frame.clone(), grayImg, CV_BGR2GRAY);

        cv::Mat edge;
        EdgeDetector edges(grayImg, edge, thresh, sigma);

        cv::imshow("Edge Detection", edge);

        std::vector<cv::Point2i> XYpoints = edges.getEdgesPoints();

        //RGB Segmentation
        cv::Mat rgbSegmentation;

        RgbSegmentation p(frame, XYpoints, sigma_blurring, sigma_c, merge_threshold, rgbSegmentation);

        cv::imshow("Rgb Segmented", rgbSegmentation);

        std::stringstream ss;
        cv::rectangle(frame, cv::Point(10, 2), cv::Point(100,20),
                  cv::Scalar(255,255,255), -1);
        ss << capture.get(1); //CV_CAP_PROP_POS_FRAMES
        std::string frameNumberString = ss.str();
        cv::putText(frame, frameNumberString.c_str(), cv::Point(15, 15),
                cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));



        //get the input from the keyboard
        keyboard = cv::waitKey( 30 );
    }
    //delete capture object
    capture.release();
}

void processSingleImg(char *image)  {
    frame = cv::imread(image);
    
    cv::imshow("Frame", frame);

    //Edge Detection
    cv::Mat grayImg;
    cv::cvtColor(frame.clone(), grayImg, CV_BGR2GRAY);

    cv::Mat edge;
    EdgeDetector edges(grayImg, edge, thresh, sigma);
    cv::imshow("Edge Detection", edge);

    std::vector<cv::Point2i> XYpoints = edges.getEdgesPoints();

    //RGB Segmentation
    cv::Mat rgbSegmentation;

    RgbSegmentation p(frame, XYpoints, sigma_blurring, sigma_c, merge_threshold, rgbSegmentation);

    cv::imshow("Rgb Segmented", rgbSegmentation);

    cv::waitKey(0);
}

/**
* @function processVideo
*/
void processVideo(char* videoFilename) {
    //create the capture object
    cv::VideoCapture capture(videoFilename);
    if(!capture.isOpened()){
        //error in opening the video input
        std::cerr << "Unable to open video file: " << videoFilename << std::endl;
        exit(EXIT_FAILURE);
    }


    //read input data. ESC or 'q' for quitting
    while( (char)keyboard != 'q' && (char)keyboard != 27 ){
        //read the current frame
        if(!capture.read(frame)) {
            std::cerr << "Unable to read next frame." << std::endl;
            std::cerr << "Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }


        cv::Mat resizedFrame(cv::Size(640, 480), CV_8UC3);
        cv::resize(frame, resizedFrame, resizedFrame.size());

        cv::imshow("Frame", resizedFrame);

        //Edge Detection
        cv::Mat grayImg;
        cv::cvtColor(resizedFrame.clone(), grayImg, CV_BGR2GRAY);

        cv::Mat edge;
        EdgeDetector edges(grayImg, edge, thresh, sigma);

        cv::imshow("Edge Detection", edge);

        std::vector<cv::Point2i> XYpoints = edges.getEdgesPoints();

        //RGB Segmentation
        cv::Mat rgbSegmentation;

        RgbSegmentation p(resizedFrame, XYpoints, sigma_blurring, sigma_c, merge_threshold, rgbSegmentation);

        cv::imshow("Rgb Segmented", rgbSegmentation);


        //get the frame number and write it on the current frame
        std::stringstream ss;
        cv::rectangle(frame, cv::Point(10, 2), cv::Point(100,20),
                  cv::Scalar(255,255,255), -1);
        ss << capture.get(1); //CV_CAP_PROP_POS_FRAMES
        std::string frameNumberString = ss.str();
        cv::putText(frame, frameNumberString.c_str(), cv::Point(15, 15),
                cv::FONT_HERSHEY_SIMPLEX, 0.5 , cv::Scalar(0,0,0));

        //get the input from the keyboard
        keyboard = cv::waitKey( 30 );
    }
    //delete capture object
    capture.release();
}


