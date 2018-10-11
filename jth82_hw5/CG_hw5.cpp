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

    vector<point> normals;
    vector<point> vertices;

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


    int index0 = mapVertex(p0,vertices,vertexIndexMapping,indexVertexMapping);
    int index1 = mapVertex(p1,vertices,vertexIndexMapping,indexVertexMapping);
    int index2 = mapVertex(p2,vertices,vertexIndexMapping,indexVertexMapping);
    int index3 = mapVertex(p3,vertices,vertexIndexMapping,indexVertexMapping);
    int index4 = mapVertex(p4,vertices,vertexIndexMapping,indexVertexMapping);
    int index5 = mapVertex(p5,vertices,vertexIndexMapping,indexVertexMapping);
    int index6 = mapVertex(p6,vertices,vertexIndexMapping,indexVertexMapping);
    int index7 = mapVertex(p7,vertices,vertexIndexMapping,indexVertexMapping);

    std::ostringstream pointVals;
    std::ostringstream coordIndexSetVals;

    pointVals << "point [" << std::endl;

    for(auto pt = vertices.begin(); pt < vertices.end(); pt++) {
        point currentPoint = *pt;
        pointVals << currentPoint.x << " " << currentPoint.y << " " << currentPoint.z << "," << std::endl;

    }

    coordIndexSetVals << "coordIndex [" << std::endl;
    coordIndexSetVals << index0 << ", " << index1 << ", " << index2 << ", " << index0 << ", " << -1 << ", " << endl;;
    coordIndexSetVals << index0 << ", " << index2 << ", " << index3 << ", " << index0 << ", " << -1 << ", " << endl;;
    coordIndexSetVals << index7 << ", " << index6 << ", " << index5 << ", " << index7 << ", " << -1 << ", " << endl;;
    coordIndexSetVals << index7 << ", " << index5 << ", " << index4 << ", " << index7 << ", " << -1 << ", " << endl;;
    coordIndexSetVals << index0 << ", " << index3 << ", " << index7 << ", " << index0 << ", " << -1 << ", " << endl;;
    coordIndexSetVals << index0 << ", " << index7 << ", " << index4 << ", " << index0 << ", " << -1 << ", " << endl;;
    coordIndexSetVals << index1 << ", " << index5 << ", " << index6 << ", " << index1 << ", " << -1 << ", " << endl;;
    coordIndexSetVals << index1 << ", " << index6 << ", " << index2 << ", " << index1 << ", " << -1 << ", " << endl;;
    coordIndexSetVals << index0 << ", " << index4 << ", " << index5 << ", " << index0 << ", " << -1 << ", " << endl;;
    coordIndexSetVals << index0 << ", " << index5 << ", " << index1 << ", " << index0 << ", " << -1 << ", " << endl;;
    coordIndexSetVals << index3 << ", " << index2 << ", " << index6 << ", " << index3 << ", " << -1 << ", " << endl;;
    coordIndexSetVals << index3 << ", " << index6 << ", " << index7 << ", " << index3 << ", " << -1 << ", " << endl;;


    if(debug) {
        cerr << "Serializing to OI format" << endl;
    }


    pointVals << "]" << std::endl;
    coordIndexSetVals << "]" << std::endl;

    Node* shapeHints = new Node("ShapeHints {","vertexOrdering        COUNTERCLOCKWISE\n","}");
    root->addChild(shapeHints);

    Node* interpolatedSeperator = new Node("Separator {","","}");

    Node* interpolatedPointsNode = new Node("Coordinate3 {",pointVals.str(),"}"); // SoCoordinate3

    Node* indexedFaceSet = new Node("IndexedFaceSet {",coordIndexSetVals.str(),"}"); // SoIndexedLineSet

    interpolatedSeperator->addChild(interpolatedPointsNode);


    interpolatedSeperator->addChild(indexedFaceSet);
    root->addChild(interpolatedSeperator);

    string ivContent = root->getString();
    ostringstream outputContent;
    outputContent << "#Inventor V2.0 ascii" << endl << ivContent;

    std::cout << outputContent.str() << std::endl;



    return 0;
}
