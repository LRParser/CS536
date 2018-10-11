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
bool memoize = false;
class point {
public:
    double x;
    double y;
    double z;

    point(double x_i, double y_i, double z_i) {
        x = x_i;
        y = y_i;
        z = z_i;
    }

    point() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }

    point operator+(point a) {
        return {x+a.x,y+a.y,z+a.z};
    }

    point operator-(point a) {
        return {x-a.x,y-a.y,z-a.z};
    }



    point operator*(double a) {
        return {a*x,a*y,a*z};
    }

    point operator/(double a) {
        return {x/a,y/a,z/a};
    }


};

bool operator<(point a, point b) {
return std::make_tuple(a.x,a.y,a.z) <  std::make_tuple(b.x, b.y, b.z);
}

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

map<int,double> factorialMap;



double fact(int k) {
    if (k == 0) {
        return 1;
    }
    else {

        if(memoize) {
            // Memoize for performance
            auto iter = factorialMap.find(k);
            if (iter != factorialMap.end()) {
                return iter->second;
            }
            else {
                double retVal = k * fact(k-1);
                factorialMap.insert(std::pair<int,double>(k,retVal));
                return retVal;
            }
        }
        else {
            double retVal = k * fact(k-1);
            return retVal;
        }


    }



}

point calculateNormal(point point1, point point2, point point3) {

    point vec1 = point1 - point2;
    point vec2 = point1 - point3;

    if(debug) {
        cerr << "Cross of: (" << vec1.x << ", " << vec1.y << ", " << vec1.z << ")" << endl;
        cerr << "  and of: (" << vec2.x << ", " << vec2.y << ", " << vec2.z << ")" <<  endl;
    }

    point cross;
    cross.z = (vec1.x*vec2.y)-(vec2.x*vec1.y);
    cross.y = (vec1.z*vec2.x)-(vec2.z*vec1.x);
    cross.x = (vec1.y*vec2.z)-(vec2.y*vec1.z);

    if(debug) {
        cerr << "Is      : " << cross.x << ", " << cross.y << ", " << cross.z << endl;
    }

    // Now we need to normalize so it becomes a unit vector
    point normalizePt;
    normalizePt.x = pow(cross.x,2);
    normalizePt.y = pow(cross.y,2);
    normalizePt.z = pow(cross.z,2);
    double pointSums = sqrt(normalizePt.x + normalizePt.y + normalizePt.z);

    cross.x = cross.x / pointSums;
    cross.y = cross.y / pointSums;
    cross.z = cross.z / pointSums;

    return cross;
}


double binomial(int k, int i) {
    return fact(k) / (fact(i)*fact(k-i));
}

double sgn(double s) {
    if(s < 0) {
        return -1;
    }
    else if(s == 0) {
        return 0;
    }
    else {
        return 1;
    }
}

double myAbs(double a) {
    if (a < 0) {
        return -1*a;
    }
    else {
        return a;
    }
}

double c(double w, double m) {
    return sgn(cos(w))*pow(myAbs(cos(w)),m);
}

double s(double w, double m) {
    return sgn(sin(w))*pow(myAbs(sin(w)),m);
}

point calcNormal(double A, double B, double C, double s1, double s2, double u, double v) {
    point normal;
    normal.x = ((float)1/A)*c(v,2-s1)*c(u,2-s2);
    normal.y = ((float)1/B)*c(v,2-s1)*s(u,2-s2);
    normal.z = ((float)1/C)*s(v,2-s1);
    return normal;

}

static int vTexIndex = 0;

int mapVertex(point p1, vector<point> & vertices, map<point,int>& vertexIndexMapping, map<int,point> & indexVertexMapping) {
    vertices.push_back(p1);
    int assignedIndex = vTexIndex;
    vertexIndexMapping[p1] = assignedIndex;
    indexVertexMapping[assignedIndex] = p1;
    vTexIndex++;
    return assignedIndex;
}

int
main(int argc, char ** argv) {

    double theta1 = -51.0;
    double theta2 = 39.0;
    double theta3 = 65.0;
    double link1 = 4.0;
    double link2 = 3.0;
    double link3 = 2.5;
    bool debug = false;
    map<point,int> vertexIndexMapping;
    map<int,point> indexVertexMapping;

    map<point,point> vertexNormalMapping;


    for (int i = 0; i < argc; i++) {
        if (std::string(argv[i]) == "-t") {

            if (i + 1 < argc) {
                theta1 = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide num_u value after -u argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-u") {

            if (i + 1 < argc) {
                theta2 = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide num_v value after -v argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-v") {

            if (i + 1 < argc) {
                theta3 = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide r value after -r argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-l") {

            if (i + 1 < argc) {
                link1 = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide r value after -r argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-m") {

            if (i + 1 < argc) {
                link2 = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide A value after -A argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-n") {

            if (i + 1 < argc) {
                link3 = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide radius value after -B argument" << std::endl;
            }
        }
        else if (std::string(argv[i]) == "-d") {
            debug = true;
        }
    }

    Node *root = new Node("", "", "");
    vector<point> baseVertices;

    // Nodes for the base

    point p0(2.000000,2.000000,1.000000);
    point p1(-2.000000, 2.000000, 1.000000);
    point p2(-2.000000, -2.000000, 1.000000);
    point p3(2.000000, -2.000000, 1.000000);
    point p4(2.000000, 2.000000, 0.000000);
    point p5(-2.000000, 2.000000, 0.000000);
    point p6(-2.000000, -2.000000, 0.000000);
    point p7(2.000000, -2.000000, 0.000000);

    int index0 = mapVertex(p0,baseVertices,vertexIndexMapping,indexVertexMapping);
    int index1 = mapVertex(p1,baseVertices,vertexIndexMapping,indexVertexMapping);
    int index2 = mapVertex(p2,baseVertices,vertexIndexMapping,indexVertexMapping);
    int index3 = mapVertex(p3,baseVertices,vertexIndexMapping,indexVertexMapping);
    int index4 = mapVertex(p4,baseVertices,vertexIndexMapping,indexVertexMapping);
    int index5 = mapVertex(p5,baseVertices,vertexIndexMapping,indexVertexMapping);
    int index6 = mapVertex(p6,baseVertices,vertexIndexMapping,indexVertexMapping);
    int index7 = mapVertex(p7,baseVertices,vertexIndexMapping,indexVertexMapping);
    vTexIndex = 0;

    std::ostringstream basePointVals;
    std::ostringstream baseIndexSetVals;
    Node* baseSeparator = new Node("Separator {","","}");
    basePointVals << "point [" << std::endl;

    for(auto pt = baseVertices.begin(); pt < baseVertices.end(); pt++) {
        point currentPoint = *pt;
        basePointVals << currentPoint.x << " " << currentPoint.y << " " << currentPoint.z << "," << std::endl;
    }

    baseIndexSetVals << "coordIndex [" << std::endl;
    baseIndexSetVals << index0 << ", " << index1 << ", " << index2 << ", " << index0 << ", " << -1 << ", " << endl;
    baseIndexSetVals << index0 << ", " << index2 << ", " << index3 << ", " << index0 << ", " << -1 << ", " << endl;
    baseIndexSetVals << index7 << ", " << index6 << ", " << index5 << ", " << index7 << ", " << -1 << ", " << endl;
    baseIndexSetVals << index7 << ", " << index5 << ", " << index4 << ", " << index7 << ", " << -1 << ", " << endl;
    baseIndexSetVals << index0 << ", " << index3 << ", " << index7 << ", " << index0 << ", " << -1 << ", " << endl;
    baseIndexSetVals << index0 << ", " << index7 << ", " << index4 << ", " << index0 << ", " << -1 << ", " << endl;
    baseIndexSetVals << index1 << ", " << index5 << ", " << index6 << ", " << index1 << ", " << -1 << ", " << endl;
    baseIndexSetVals << index1 << ", " << index6 << ", " << index2 << ", " << index1 << ", " << -1 << ", " << endl;
    baseIndexSetVals << index0 << ", " << index4 << ", " << index5 << ", " << index0 << ", " << -1 << ", " << endl;
    baseIndexSetVals << index0 << ", " << index5 << ", " << index1 << ", " << index0 << ", " << -1 << ", " << endl;
    baseIndexSetVals << index3 << ", " << index2 << ", " << index6 << ", " << index3 << ", " << -1 << ", " << endl;
    baseIndexSetVals << index3 << ", " << index6 << ", " << index7 << ", " << index3 << ", " << -1 << ", " << endl;

    basePointVals << "]" << std::endl;
    baseIndexSetVals << "]" << std::endl;

    Node* basePoints = new Node("Coordinate3 {",basePointVals.str(),"}"); // SoCoordinate3
    Node* baseFaces = new Node("IndexedFaceSet {",baseIndexSetVals.str(),"}"); // SoIndexedLineSet

    baseSeparator->addChild(basePoints);
    baseSeparator->addChild(baseFaces);
    root->addChild(baseSeparator);


    // Nodes for link 1
    std::ostringstream link1PointVals;
    std::ostringstream link1IndexSetVals;
    vector<point> link1Vertices;
    Node* link1Separator = new Node("Separator {","","}");

    point p8(0.703233, -0.073913, 5.000000);
    point p9(0.073913, 0.703233, 5.000000);
    point p10(-0.703233, 0.073913, 5.000000);
    point p11(-0.073913, -0.703233, 5.000000);
    point p12(0.703233, -0.073913, 1.000000);
    point p13(0.073913, 0.703233, 1.000000);
    point p14(-0.703233, 0.073913, 1.000000);
    point p15(-0.073913, -0.703233, 1.000000);

    int index8 = mapVertex(p8,link1Vertices,vertexIndexMapping,indexVertexMapping);
    int index9 = mapVertex(p9,link1Vertices,vertexIndexMapping,indexVertexMapping);
    int index10 = mapVertex(p10,link1Vertices,vertexIndexMapping,indexVertexMapping);
    int index11 = mapVertex(p11,link1Vertices,vertexIndexMapping,indexVertexMapping);
    int index12 = mapVertex(p12,link1Vertices,vertexIndexMapping,indexVertexMapping);
    int index13 = mapVertex(p13,link1Vertices,vertexIndexMapping,indexVertexMapping);
    int index14 = mapVertex(p14,link1Vertices,vertexIndexMapping,indexVertexMapping);
    int index15 = mapVertex(p15,link1Vertices,vertexIndexMapping,indexVertexMapping);

    link1PointVals << "point [" << std::endl;
    for(auto pt = link1Vertices.begin(); pt < link1Vertices.end(); pt++) {
        point currentPoint = *pt;
        link1PointVals << currentPoint.x << " " << currentPoint.y << " " << currentPoint.z << "," << std::endl;
    }

    link1PointVals << "]" << std::endl;


    link1IndexSetVals << "coordIndex [" << std::endl;
    link1IndexSetVals << index8 << ", " << index9 << ", " << index10 << ", " << index8 << ", " << -1 << ", " << endl;
    link1IndexSetVals << index8 << ", " << index10 << ", " << index11 << ", " << index8 << ", " << -1 << ", " << endl;
    link1IndexSetVals << index15 << ", " << index14 << ", " << index13 << ", " << index15 << ", " << -1 << ", " << endl;
    link1IndexSetVals << index15 << ", " << index13 << ", " << index12 << ", " << index15 << ", " << -1 << ", " << endl;
    link1IndexSetVals << index8 << ", " << index11 << ", " << index15 << ", " << index8 << ", " << -1 << ", " << endl;
    link1IndexSetVals << index8 << ", " << index15 << ", " << index12 << ", " << index8 << ", " << -1 << ", " << endl;
    link1IndexSetVals << index9 << ", " << index13 << ", " << index14 << ", " << index9 << ", " << -1 << ", " << endl;
    link1IndexSetVals << index9 << ", " << index14 << ", " << index10 << ", " << index9 << ", " << -1 << ", " << endl;
    link1IndexSetVals << index8 << ", " << index12 << ", " << index13 << ", " << index8 << ", " << -1 << ", " << endl;
    link1IndexSetVals << index8 << ", " << index13 << ", " << index9 << ", " << index8 << ", " << -1 << ", " << endl;
    link1IndexSetVals << index11 << ", " << index10 << ", " << index14 << ", " << index11 << ", " << -1 << ", " << endl;
    link1IndexSetVals << index11 << ", " << index14 << ", " << index15 << ", " << index11 << ", " << -1 << ", " << endl;
    link1IndexSetVals << "]" << std::endl;

    Node* link1Points = new Node("Coordinate3 {",link1PointVals.str(),"}"); // SoCoordinate3
    Node* link1Faces = new Node("IndexedFaceSet {",link1IndexSetVals.str(),"}"); // SoIndexedLineSet

    link1Separator->addChild(link1Points);
    link1Separator->addChild(link1Faces);
    root->addChild(link1Separator);

    // Nodes for link 2

    // Nodes for link 3


    if(debug) {
        cerr << "Serializing to OI format" << endl;
    }



    Node* shapeHints = new Node("ShapeHints {","vertexOrdering        COUNTERCLOCKWISE\n","}");
    root->addChild(shapeHints);




    string ivContent = root->getString();
    ostringstream outputContent;
    outputContent << "#Inventor V2.0 ascii" << endl << ivContent;

    std::cout << outputContent.str() << std::endl;



    return 0;
}
