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



float binomial(int k, int i) {
    return fact(k) / (fact(i)*fact(k-i));
}


int
main(int argc, char ** argv)
{

    string fName;
    int num_u = 11;
    int num_v = 11;
    float radius = 0.1;
    bool shadeWithNormals = false;


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
                num_u = std::stoi(std::string(argv[i + 1]));
            }
            else {
                std::cerr << "Must provide num_u value after -u argument" << std::endl;
            }
        }
        else if (std::string(argv[i]) == "-v") {

            if (i + 1 < argc) {
                num_v = std::stoi(std::string(argv[i + 1]));
            }
            else {
                std::cerr << "Must provide num_v value after -v argument" << std::endl;
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
        else if (std::string(argv[i]) == "-F") {

            shadeWithNormals = false;
        }
        else if (std::string(argv[i]) == "-S") {
            shadeWithNormals = true;
        }
    }

    double du = 1.0/(float)(num_u-1);
    double dv = 1.0/(float)(num_v-1);

    if (fName.empty()) {
        fName = "patchPoints.txt";
    }

    if (debug) {
        std::cout << "Reading from file: " << fName << std::endl;
    }
    std::ifstream input(fName.c_str());
    if (input.fail()) {
        std::cerr << "Failed to open" << std::endl;
    }

    point* k = new point[4][4];


    vector<point> currentJ;
    string currentLine;
    int i = 0;
    int j = 0;
    while(std::getline(input, currentLine)) {

        if(debug) {
            std::cout << currentLine << std::endl;
        }

        point point1;
        std::istringstream ss(currentLine);
        ss >> point1.x >> point1.y >> point1.z;

        k[i][j] = point1;

        i = i + 1;

        if (i == 3) {
            i = 0;
            j = j + 1;
        }

        if (debug) {
            std::cout << "Parsed: " << point1.x << " " << point1.y << " " << point1.z << std::endl;
        }
    }

    if (debug) {
        std::cout << "Done" << std::endl;
    }

    Node* root = new Node("","","");

    // Interpolate the curve


    if(debug) {
        cout << "Number of points is: " << k << endl;
    }
    vector<point> calcPoints;

    int n = 3;
    int m = 3;

    bool interpolatedUEnd = false;

    for (int u = 1; u < num_u; u++) {

        float uParam = (float)u/float(num_u);

        bool interpolatedVEnd = false;
        while(v <= 1.0 && !interpolatedVEnd) {

            point currentPoint;
            currentPoint.x = 0.0;
            currentPoint.y = 0.0;
            currentPoint.z = 0.0;

            for(int i = 0; i <= 3; i++) {

                double factorU = binomial(n, i) * pow(u,i) * pow(1-u,n-i) ;

                for(int j=0; j <= 3; j++) {

                    double factorV = binomial(n, j) * pow(v,j) * pow(1-v,n-j);
                    point controlPoint = k.at(j).at(i);
                    currentPoint = currentPoint + (controlPoint * factorU * factorV);
                }
            }

            calcPoints.push_back(currentPoint);

            if (v == 1.0) {
                interpolatedVEnd = true;
            }
            if (v + dv > 1 && !interpolatedVEnd) {
                v = 1.0;
            }
            else {
                v += dv;
            }

        }


        if (u == 1.0) {
            interpolatedUEnd = true;
        }
        if (u + du > 1 && !interpolatedUEnd) {
            u = 1.0;
        }
        else {
            u += du;
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

    for(int j = 0; j <= 3; j++) {

        for(int i = 0; i <=3; i++) {

            Node* cpSep = new Node("Separator {","","}");

            point point1 = k.at(j).at(i);
            double x = point1.x;
            double y = point1.y;
            double z = point1.z;

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

    }





    string ivContent = root->getString();
    ostringstream outputContent;
    outputContent << "#Inventor V2.0 ascii" << endl << ivContent;

    std::cout << outputContent.str() << std::endl;

    return 0;
}
