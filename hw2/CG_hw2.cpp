/*
 * OITest.c
 *
 *  Created on: Sep 15, 2018
 *      Author: joe
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include "Node.h"

using namespace std;

bool debug = false;

struct point {
    double x;
    double y;
    double z;

    point operator+(point a) {
        return {x+a.x,y+a.y,z+a.z};
    }

    point operator-(point a) {
        return {x-a.x,y-a.y,z-a.z};
    }

    point operator*(double a) {
        return {a*x,a*y,a*z};
    }

};

point pointMult(double factor, point srcPt) {
    point newPoint;
    newPoint.x = srcPt.x * factor;
    newPoint.y = srcPt.y * factor;
    newPoint.z = srcPt.z * factor;
    return newPoint;
}

point pointAdd(point point1, point point2) {
    point point3;
    point3.x = point1.x + point2.x;
    point3.y = point1.y + point2.y;
    point3.z = point1.z + point2.z;
    return point3;
}

map<int,float> factorialMap;


float fact(int k) {
    if (k == 0) {
        return 1;
    }

    // Memoize for performance
    auto iter = factorialMap.find(k);
    if (iter != factorialMap.end()) {
       return iter->second;
    }
    else {
        float retVal = k * fact(k-1);
        factorialMap.insert(std::pair<int,int>(k,retVal));
        return retVal;
    }

}



float kchoosei(int k, int i) {
    return fact(k) / (fact(i)*fact(k-i));
}


int
main(int argc, char ** argv)
{

    string fName;
    float du = 0.09;
    float radius = 0.1;

    for(int i=0; i < argc; i++) {
        if (std::string(argv[i]) == "-f") {

            if (i + 1 < argc) {
                fName = std::string(argv[i + 1]);
            }
            else {
                std::cerr << "Must provide file name after -f argument" << std::endl;
            }
        }
        else if (std::string(argv[i]) == "-u") {

            if (i + 1 < argc) {
                du = std::stof(std::string(argv[i + 1]));
            }
            else {
                std::cerr << "Must provide du value after -u argument" << std::endl;
            }
        }
        else if (std::string(argv[i]) == "-r") {

            if (i + 1 < argc) {
                radius = std::stof(std::string(argv[i + 1]));
            }
            else {
                std::cerr << "Must provide radius value after -r argument" << std::endl;
            }
        }
    }

    if (fName.empty()) {
        fName = "cpts_in.txt";
    }

    if (debug) {
        std::cout << "Reading from file: " << fName << std::endl;
    }
    std::ifstream input(fName.c_str());
    if (input.fail()) {
        std::cerr << "Failed to open" << std::endl;
    }
    point tangent0;
    point tangent1;
    vector<point> parsedPoints;

    string currentLine;
    int i = 0;
    while(std::getline(input, currentLine)) {

        if(debug) {
            std::cout << currentLine << std::endl;
        }
        point point1;
        std::istringstream ss(currentLine);
        ss >> point1.x >> point1.y >> point1.z;

        if (i == 0) {
            tangent0 = point1;

        }
        else if(i == 1) {
            tangent1 = point1;

        }
        else {
            parsedPoints.push_back(point1);
        }

        i = i + 1;

        if (debug) {
            std::cout << "Parsed: " << point1.x << " " << point1.y << " " << point1.z << std::endl;
        }
    }

    if (debug) {
        std::cout << "Done" << std::endl;
    }

    // Convert the tangents and points (which are in Hermite form) into Bezier form
    // We specify the tangents at the first and last input points
    vector<vector<point>> curves;

    int parsedLen = (int) parsedPoints.size();

    for(int i = 0; i < parsedLen; i++) {

        vector<point> points;

        if (i < parsedLen - 1) {

            point pk = parsedPoints.at(i);

            point pkplus1 = parsedPoints.at(i+1);

            point t0;
            point t1;

            if (i == 0) {
                t0 = tangent0;
            }
            else {
                point pkminus1 = parsedPoints.at(i-1);
                t0 = (pkplus1-pkminus1)*0.5;
            }

            if (i + 2 == parsedLen) {
                t1 = tangent1;
            }
            else {
                point pkplus2 = parsedPoints.at(i+2);
                t1 = (pkplus2-pk)*0.5;
            }

            points.push_back(pk);
            points.push_back(pk + (t0*(1/3.0)));
            points.push_back(pkplus1 - (t1*(1/3.0)));
            points.push_back(pkplus1);

        }

        if (points.size() > 0) {
            curves.push_back(points);
        }
    }


    Node* root = new Node("","","");

    vector<point> calcPoints;

    for (auto cIt = curves.begin(); cIt != curves.end(); cIt++) {

        vector<point> points = *cIt;

        float u = 0.0;
        int k = (int) points.size() - 1;

        if(debug) {
            cout << "Number of points is: " << k << endl;
        }

        bool interpolatedEnd = false;
        while(u <= 1.0 && !interpolatedEnd) {

        point currentPoint;
        currentPoint.x = 0.0;
        currentPoint.y = 0.0;
        currentPoint.z = 0.0;

        for(int i = 0; i <= k; i++) {
            point controlPoint = points.at(i);

            double factor = kchoosei(k, i) * pow(1-u,k-i) * pow(u,i);
            if(debug) {
                cout << "     i is: " << i << " factor is: " << factor << endl;
            }

            point calcPoint = pointMult(factor,controlPoint);
            currentPoint = pointAdd(currentPoint,calcPoint);

        }

        calcPoints.push_back(currentPoint);

        if (u == 1.0) {
            interpolatedEnd = true;
        }
        if (u + du > 1 && !interpolatedEnd) {
            u = 1.0;
        }
        else {
            u += du;
        }
    }

    }


    std::ostringstream pointVals;
    std::ostringstream coordIndexSetVals;

    pointVals << "point [" << std::endl;
    coordIndexSetVals << "coordIndex [" << std::endl;

    i = 0;
    for(auto it = calcPoints.begin(); it != calcPoints.end(); it++) {
        pointVals << it->x << " " << it->y << " " << it->z << "," << std::endl;
        coordIndexSetVals << i << ", ";

        i++;
    }

    pointVals << "]" << std::endl;
    coordIndexSetVals << -1 << ", " << std::endl;
    coordIndexSetVals << "]" << std::endl;

    Node* interpolatedSeperator = new Node("Separator {","","}");

    Node* interpolatedPoints = new Node("Coordinate3 {",pointVals.str(),"}"); // SoCoordinate3
    Node* indexedLineSet = new Node("IndexedLineSet {",coordIndexSetVals.str(),"}"); // SoIndexedLineSet

    interpolatedSeperator->addChild(interpolatedPoints);
    interpolatedSeperator->addChild(indexedLineSet);
    root->addChild(interpolatedSeperator);

    // Plot the control points
    for (auto it = parsedPoints.begin(); it != parsedPoints.end(); it++) {

        Node* cpSep = new Node("Separator {","","}");

        double x = it->x;
        double y = it->y;
        double z = it->z;

        if(debug) {
            std::cout << "Placing: " << x << " " << y << " " << z << std::endl;
        }

        std::ostringstream transformStr;
        transformStr << "translation " << x << " " << y << " " << z << std::endl;
        Node* translation = new Node("Transform {",transformStr.str(),"}");
        cpSep->addChild(translation);

        std::ostringstream radiusStr;
        radiusStr << "radius  " << radius << std::endl;
        Node* sphere = new Node("Sphere {",radiusStr.str(),"}");
        cpSep->addChild(sphere);

        root->addChild(cpSep);

    }

    string ivContent = root->getString();
    ostringstream outputContent;
    outputContent << "#Inventor V2.0 ascii" << endl << ivContent;

    std::cout << outputContent.str() << std::endl;

    return 0;
}
