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

    if(debug) {
        cout << "num_u is: " << num_u;
        cout << "num_v is: " << num_u;

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

    vector<vector<point>> k;
    // Allocate space for 16 control points
    for(int i = 0; i < 4; i++) {

        vector<point> iVec;

        for(int j = 0; j < 4; j++) {
            point p1;
            iVec.push_back(p1);
        }

        k.push_back(iVec);
    }


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

        if(debug) {
            cout << "Set at index: " << i << " and " << j << endl;
        }
        k[i][j] = point1;

        if (i == 3) {
            i = 0;
            j = j + 1;
            continue;
        }
        else {
            i = i + 1;
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

    // Ensure we have num_u * num_v points in interpolatedPoints
    vector<vector<point>> interpolatedPoints;
    for (int u = 0; u < num_u; u++) {

        vector<point> uVec;

        for(int v = 0; v < num_v; v++) {
            point p;
            uVec.push_back(p);
        }

        interpolatedPoints.push_back(uVec);
    }

    if(debug) {
        cout << "Created placeholder for points" << endl;
    }


    int n = 3;
    int m = 3;

    for (int u = 0; u < num_u; u++) {

        float uParam = (float)u/float(num_u);

        for(int v = 0; v < num_v; v++) {

            float vParam = float(v)/float(num_v);

            point currentPoint;
            currentPoint.x = 0.0;
            currentPoint.y = 0.0;
            currentPoint.z = 0.0;

            for(int i = 0; i <= n; i++) {

                double factorU = binomial(n, i) * pow(uParam,i) * pow(1-uParam,n-i) ;

                for(int j=0; j <= m; j++) {

                    double factorV = binomial(m, j) * pow(vParam,j) * pow(1-vParam,m-j);
                    if(debug) {
                        cout << " Get control point at: " << j << "," << i << endl;
                    }
                    point controlPoint = k[j][i];
                    currentPoint = currentPoint + (controlPoint * factorU * factorV);
                }
            }

            if(debug) {
                cout << " Set interpolated point at: " << u << "," << v << endl;
            }

            interpolatedPoints[u][v] = currentPoint;
        }

    }

    if(debug) {
        cout << "Writing to OpenInventor format" << endl;
    }

    std::ostringstream pointVals;
    std::ostringstream coordIndexSetVals;

    pointVals << "point [" << std::endl;
    coordIndexSetVals << "coordIndex [" << std::endl;


    int vertexIndex = 0;
    for(int i = 0; i < num_u - 1; i++) {
        for(int j = 0; j < num_v - 1; j++) {
            point vertex1 =  interpolatedPoints[i][j]; // 1, or 0,0
            point vertex2 = interpolatedPoints[i+1][j]; // 2, 1,0
            point vertex3 = interpolatedPoints[i][j+1]; // 3, 0,1

            point vertex4 = interpolatedPoints[i+1][j]; // 4, 1,0
            point vertex5 = interpolatedPoints[i+1][j+1]; // 5, 1,1
            point vertex6 = interpolatedPoints[i][j+1]; // 6, 0,1

            pointVals << vertex1.x << " " << vertex1.y << " " << vertex1.z << "," << std::endl;
            coordIndexSetVals << vertexIndex << ", ";
            vertexIndex++;

            pointVals << vertex2.x << " " << vertex2.y << " " << vertex2.z << "," << std::endl;
            coordIndexSetVals << vertexIndex << ", ";
            vertexIndex++;

            pointVals << vertex3.x << " " << vertex3.y << " " << vertex3.z << "," << std::endl;
            coordIndexSetVals << vertexIndex << ", ";
            vertexIndex++;

            coordIndexSetVals << -1 << ", " << endl;

            pointVals << vertex4.x << " " << vertex4.y << " " << vertex4.z << "," << std::endl;
            coordIndexSetVals << vertexIndex << ", ";
            vertexIndex++;

            pointVals << vertex5.x << " " << vertex5.y << " " << vertex5.z << "," << std::endl;
            coordIndexSetVals << vertexIndex << ", ";
            vertexIndex++;

            pointVals << vertex6.x << " " << vertex6.y << " " << vertex6.z << "," << std::endl;
            coordIndexSetVals << vertexIndex << ", ";
            vertexIndex++;

            coordIndexSetVals << -1 << ", " << endl;


        }
    }

    if(debug) {
        cout << "Serialized to OI format" << endl;
    }

    i = 0;

    pointVals << "]" << std::endl;
    coordIndexSetVals << "]" << std::endl;

    Node* shapeHints = new Node("ShapeHints {","vertexOrdering        COUNTERCLOCKWISE","}");
    root->addChild(shapeHints);

    Node* interpolatedSeperator = new Node("Separator {","","}");

    Node* interpolatedPointsNode = new Node("Coordinate3 {",pointVals.str(),"}"); // SoCoordinate3
    Node* indexedFaceSet = new Node("IndexedFaceSet {",coordIndexSetVals.str(),"}"); // SoIndexedLineSet

    interpolatedSeperator->addChild(interpolatedPointsNode);
    interpolatedSeperator->addChild(indexedFaceSet);
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
