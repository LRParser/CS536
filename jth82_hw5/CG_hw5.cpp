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

class quadresult {
public:
    Node* sep;
    vector<point> vertices;
    quadresult(Node* sep_i, vector<point> vertices_i) {
        sep=sep_i;
        vertices = vertices_i;
    }
};

Node* constructQuadNode(vector<point> vertices) {


    std::ostringstream pointVals;
    std::ostringstream indexSetVals;
    Node* baseSeparator = new Node("Separator {","","}");
    pointVals << "point [" << std::endl;

    for(auto pt = vertices.begin(); pt < vertices.end(); pt++) {
        point currentPoint = *pt;
        pointVals << currentPoint.x << " " << currentPoint.y << " " << currentPoint.z << "," << std::endl;
    }

    // For the 6 faces in a quad we have 12 triangles, as shown below
    indexSetVals << "coordIndex [" << std::endl;
    indexSetVals << 0 << ", " << 1 << ", " << 2 << ", " << 0 << ", " << -1 << ", " << endl;
    indexSetVals << 0 << ", " << 2 << ", " << 3 << ", " << 0 << ", " << -1 << ", " << endl;
    indexSetVals << 7 << ", " << 6 << ", " << 5 << ", " << 7 << ", " << -1 << ", " << endl;
    indexSetVals << 7 << ", " << 5 << ", " << 4 << ", " << 7 << ", " << -1 << ", " << endl;
    indexSetVals << 0 << ", " << 3 << ", " << 7 << ", " << 0 << ", " << -1 << ", " << endl;
    indexSetVals << 0 << ", " << 7 << ", " << 4 << ", " << 0 << ", " << -1 << ", " << endl;
    indexSetVals << 1 << ", " << 5 << ", " << 6 << ", " << 1 << ", " << -1 << ", " << endl;
    indexSetVals << 1 << ", " << 6 << ", " << 2 << ", " << 1 << ", " << -1 << ", " << endl;
    indexSetVals << 0 << ", " << 4 << ", " << 5 << ", " << 0 << ", " << -1 << ", " << endl;
    indexSetVals << 0 << ", " << 5 << ", " << 1 << ", " << 0 << ", " << -1 << ", " << endl;
    indexSetVals << 3 << ", " << 2 << ", " << 6 << ", " << 3 << ", " << -1 << ", " << endl;
    indexSetVals << 3 << ", " << 6 << ", " << 7 << ", " << 3 << ", " << -1 << ", " << endl;

    pointVals << "]" << std::endl;
    indexSetVals << "]" << std::endl;

    Node* basePoints = new Node("Coordinate3 {",pointVals.str(),"}"); // SoCoordinate3
    Node* baseFaces = new Node("IndexedFaceSet {",indexSetVals.str(),"}"); // SoIndexedLineSet

    baseSeparator->addChild(basePoints);
    baseSeparator->addChild(baseFaces);
    return baseSeparator;
}

point rotateVector(vector<vector<double>> matrix, point rhs) {
    point retval;
    retval.x = (matrix[0][0] * rhs.x) + (matrix[0][1] * rhs.y) + (matrix[0][2] * rhs.z);
    retval.y = (matrix[1][0] * rhs.x) + (matrix[1][1] * rhs.y) + (matrix[1][2] * rhs.z);
    retval.z = (matrix[2][0] * rhs.x) + (matrix[2][1] * rhs.y) + (matrix[2][2] * rhs.z);
    return retval;
}

double to_radians(double degree) {
    return degree * (M_PI/180.0);
}

vector<point> translateZ(vector<point> inputPoints, double zOffset) {
    vector<point> returnPoints;
    for(auto it = inputPoints.begin(); it < inputPoints.end(); it++) {
        point p = *it;
        p.z += zOffset;
        returnPoints.push_back(p);
    }
    return returnPoints;
}

vector<point> drawQuad(point ll, point ur) {
    vector<point> returnVec;

    point p0(ur.x, ur.y, ur.z); // point p0(2.000000,2.000000,1.000000); // upper right
    point p1(ll.x, ur.y, ur.z);
    point p2(ll.x, ll.y, ur.z);
    point p3(ur.x, ll.y, ur.z);
    point p4(ur.x, ur.y, ll.z);
    point p5(ll.x, ur.y, ll.z);
    point p6(ll.x, ll.y, ll.z); // point p6(-2.000000, -2.000000, 0.000000); // lower left
    point p7(ur.x, ll.y, ll.z);


    returnVec.push_back(p0);
    returnVec.push_back(p1);
    returnVec.push_back(p2);
    returnVec.push_back(p3);
    returnVec.push_back(p4);
    returnVec.push_back(p5);
    returnVec.push_back(p6);
    returnVec.push_back(p7);

    return returnVec;

}

int
main(int argc, char ** argv) {

    double theta1 = -51.0;
    double theta2 = 39.0;
    double theta3 = 65.0;
    double link1Length = 4.0;
    double link2Length = 3.0;
    double link3Length = 2.5;
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
                link1Length = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide r value after -r argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-m") {

            if (i + 1 < argc) {
                link2Length = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide A value after -A argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-n") {

            if (i + 1 < argc) {
                link3Length = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide radius value after -B argument" << std::endl;
            }
        }
        else if (std::string(argv[i]) == "-d") {
            debug = true;
        }
    }

    Node *root = new Node("", "", "");

    float zTranslation = 0.0;

    // Draw base
    point baseLL(-2,-2,0);
    point baseUR(2,2,1);
    vector<point> baseVertices = drawQuad(baseLL, baseUR);
    Node* baseSeparator = constructQuadNode(baseVertices);
    root->addChild(baseSeparator);
    zTranslation += 1;

    // Draw link 1
    point l1LL(-.5,-.5,0);
    point l1UR(.5,.5,link1Length);
    vector<point> link1Vertices = drawQuad(l1LL, l1UR);
    link1Vertices = translateZ(link1Vertices, zTranslation);
    Node* link1Separator = constructQuadNode(link1Vertices);
    root->addChild(link1Separator);
    zTranslation += link1Length;

    // Draw link 2
    float l2 = 3.0;
    theta2 = 37.0;
    double theta2Radians = to_radians(theta2);

    // Now to apply a rotation, create rotation matrix
    vector<double> row1;
    row1.push_back(cos(theta2Radians));
    row1.push_back(0);
    row1.push_back(sin(theta2Radians));

    vector<double> row2;
    row2.push_back(0);
    row2.push_back(1);
    row2.push_back(0);

    vector<double> row3;
    row3.push_back(-1 * sin(theta2Radians));
    row3.push_back(0);
    row3.push_back(cos(theta2Radians));

    vector<vector<double>> rotationMatrixLink2;
    rotationMatrixLink2.push_back(row1);
    rotationMatrixLink2.push_back(row2);
    rotationMatrixLink2.push_back(row3);

    point l2LL(-.5,-.5,0);
    point l2UR(.5,.5,l2);

    vector<point> link2Vertices = drawQuad(l2LL, l2UR);

    vector<point> link2VerticesRotated;
    for(auto it = link2Vertices.begin(); it < link2Vertices.end(); it++) {
        point p = *it;
        link2VerticesRotated.push_back(rotateVector(rotationMatrixLink2, p));
    }

    vector<point> link2VerticesRotatedTranslated = translateZ(link2VerticesRotated, zTranslation);
    
    Node* link2Separator = constructQuadNode(link2VerticesRotatedTranslated);
    root->addChild(link2Separator);
    zTranslation += l2;
    /*

    // Old method for base

    point p0(2.000000,2.000000,1.000000); // upper right
    point p1(-2.000000, 2.000000, 1.000000);
    point p2(-2.000000, -2.000000, 1.000000);
    point p3(2.000000, -2.000000, 1.000000);
    point p4(2.000000, 2.000000, 0.000000);
    point p5(-2.000000, 2.000000, 0.000000);
    point p6(-2.000000, -2.000000, 0.000000); // lower left
    point p7(2.000000, -2.000000, 0.000000);

    vector<point> baseVertices;

    baseVertices.push_back(p0);
    baseVertices.push_back(p1);
    baseVertices.push_back(p2);
    baseVertices.push_back(p3);
    baseVertices.push_back(p4);
    baseVertices.push_back(p5);
    baseVertices.push_back(p6);
    baseVertices.push_back(p7);

    // Nodes for link 1

    vector<point> link1Vertices;

    point p8(0.703233, -0.073913, 5.000000);
    point p9(0.073913, 0.703233, 5.000000);
    point p10(-0.703233, 0.073913, 5.000000);
    point p11(-0.073913, -0.703233, 5.000000);
    point p12(0.703233, -0.073913, 1.000000);
    point p13(0.073913, 0.703233, 1.000000);
    point p14(-0.703233, 0.073913, 1.000000);
    point p15(-0.073913, -0.703233, 1.000000);

    link1Vertices.push_back(p8);
    link1Vertices.push_back(p9);
    link1Vertices.push_back(p10);
    link1Vertices.push_back(p11);
    link1Vertices.push_back(p12);
    link1Vertices.push_back(p13);
    link1Vertices.push_back(p14);
    link1Vertices.push_back(p15);

    Node* link1Separator = constructQuadNode(link1Vertices);
    root->addChild(link1Separator);

    // Nodes for link 2

    vector<point> link2Vertices;

    point p16(1.821242, -1.454539, 7.016778);
    point p17(1.332169, -0.850583, 7.646098);
    point p18(0.555023, -1.479904, 7.646098);
    point p19(1.044096, -2.083860, 7.016778);
    point p20(0.633110, 0.012682, 4.685340);
    point p21(0.144036, 0.616638, 5.314660);
    point p22(-0.633110, -0.012682, 5.314660);
    point p23(-0.144036, -0.616638, 4.685340);

    link2Vertices.push_back(p16);
    link2Vertices.push_back(p17);
    link2Vertices.push_back(p18);
    link2Vertices.push_back(p19);
    link2Vertices.push_back(p20);
    link2Vertices.push_back(p21);
    link2Vertices.push_back(p22);
    link2Vertices.push_back(p23);

    Node* link2Separator = constructQuadNode(link2Vertices);
    root->addChild(link2Separator);

    // Nodes for link 3
    vector<point> link3Vertices;

    point p24(3.027149, -2.943710, 6.241485);
    point p25(3.179396, -3.131719, 7.211781);
    point p26(2.402250, -3.761039, 7.211781);
    point p27(2.250004, -3.573031, 6.241485);
    point p28(1.500582, -1.058557, 6.846290);
    point p29(1.652829, -1.246566, 7.816586);
    point p30(0.875683, -1.875886, 7.816586);
    point p31(0.723436, -1.687877, 6.846290);

    link3Vertices.push_back(p24);
    link3Vertices.push_back(p25);
    link3Vertices.push_back(p26);
    link3Vertices.push_back(p27);
    link3Vertices.push_back(p28);
    link3Vertices.push_back(p29);
    link3Vertices.push_back(p30);
    link3Vertices.push_back(p31);

    Node* link3Separator = constructQuadNode(link3Vertices);
    root->addChild(link3Separator);
    */

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
