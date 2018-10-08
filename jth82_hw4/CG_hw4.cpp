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
map<point,int> vertexIndexMapping;
map<point,point> vertexNormalMapping;


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

    /*
    point vec1 = point1 - point2;
    point vec2 = point1 - point3;


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
    */
    point p;
    p.x = 0;
    p.y = 0;
    p.z = 0;
    return p;
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

int
main(int argc, char ** argv)
{

    int num_u = 19;
    int num_v = 9;
    double s1 = 1.0;
    double s2 = 1.0;
    double A = 1.0;
    double B = 1.0;
    double C = 1.0;

    bool shadeWithNormals = false;


    for(int i=0; i < argc; i++) {
        if (std::string(argv[i]) == "-u") {

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
                s1 = std::stod(std::string(argv[i + 1]));
            }
            else {
                std::cerr << "Must provide r value after -r argument" << std::endl;
            }
        }
        else if (std::string(argv[i]) == "-t") {

            if (i + 1 < argc) {
                s2 = std::stod(std::string(argv[i + 1]));
            }
            else {
                std::cerr << "Must provide r value after -r argument" << std::endl;
            }
        }
        else if (std::string(argv[i]) == "-A") {

            if (i + 1 < argc) {
                A = std::stod(std::string(argv[i + 1]));
            }
            else {
                std::cerr << "Must provide A value after -A argument" << std::endl;
            }
        }
        else if (std::string(argv[i]) == "-B") {

            if (i + 1 < argc) {
                B = std::stod(std::string(argv[i + 1]));
            }
            else {
                std::cerr << "Must provide radius value after -B argument" << std::endl;
            }
        }
        else if (std::string(argv[i]) == "-C") {

            if (i + 1 < argc) {
                C = std::stod(std::string(argv[i + 1]));
            }
            else {
                std::cerr << "Must provide radius value after -C argument" << std::endl;
            }
        }
        else if (std::string(argv[i]) == "-F") {

            shadeWithNormals = false;
        }
        else if (std::string(argv[i]) == "-S") {
            shadeWithNormals = true;
        }
        else if (std::string(argv[i]) == "-d") {
            debug = true;
        }
    }

    if(debug) {
        cout << "num_u is: " << num_u;
        cout << "num_v is: " << num_v;

    }


    Node* root = new Node("","","");

    // Interpolate the shape

    double uRange = M_PI - -M_PI;
    double vRange = M_PI/2 - -(M_PI/2);
    double du = uRange / double(num_u);
    double dv = vRange / double(num_v);

    if(debug) {
        cout << "uRange is: " << uRange << "vRange is: " << vRange << "du is: " << du << " dv is: " << dv << endl;
    }


    int numPts = 0;
    int numUPts = 0;
    int numVPts = 0;
    vector<vector<point>> interpolatedPoints;

    while(numUPts < num_u) {
        vector<point> uVec;

        numVPts = 0;

        while(numVPts < num_v) {
            point p;
            uVec.push_back(p);
            numVPts++;
            numPts++;
        }
        numUPts++;
        interpolatedPoints.push_back(uVec);

    }

    if(debug) {
        cout << "Created placeholder for " << numPts << " points" << endl;
    }

    std::ostringstream pointVals;
    std::ostringstream coordIndexSetVals;

    pointVals << "point [" << std::endl;


    int uIdx = 0;
    int vIdx = 0;
    double u =  -M_PI;
    double v = -M_PI/2;

    while(u <= M_PI) {

        while(v <= (M_PI/2)) {


            if(debug) {
                cout << "u is: " << u << " ";
                cout << "v is: " << v << " ";
            }


            point currentPoint;
            currentPoint.x = A*c(v,s1)*c(u,s2);
            currentPoint.y = B*c(v,s1)*s(u,s2);
            currentPoint.z = C*s(v,s1);


            if(debug) {
                cout << " Set interpolated point at: " << uIdx << "," << vIdx << " to x: " << currentPoint.x << ", y: " << currentPoint.y << ", z: " << currentPoint.z << endl;
            }


            interpolatedPoints[uIdx][vIdx] = currentPoint;

            vIdx++;
            v = v+dv;
        }

        v = -M_PI/2;
        vIdx = 0;
        uIdx++;
        u=u+du;

    }

    if(debug) {
        cout << "numPts is" << numPts << endl;
        cout << "Writing to OpenInventor format" << endl;
    }


    coordIndexSetVals << "coordIndex [" << std::endl;


    vector<point> normals;
    vector<point> vertices;

    numPts = 0;
    for(int i = 0; i < num_u; i++) {
        for(int j = 0; j < num_v; j++) {

            point currentPoint =  interpolatedPoints[i][j]; // 0, or 0,0; 1
            pointVals << currentPoint.x << " " << currentPoint.y << " " << currentPoint.z << "," << std::endl;
            vertexIndexMapping.insert(std::pair<point,int>(currentPoint,numPts));
            vertices.push_back(currentPoint);

            numPts++;
        }
    }


    for(int i = 0; i < num_u; i++) {
        for(int j = 0; j < num_v; j++) {

            // Four distinct points become a patch (two tesselated triangles)
            point vertex0 =  interpolatedPoints[i][j]; // 0
            int index0 = vertexIndexMapping.at(vertex0);

            int indexiplus1 = i+1;
            if (indexiplus1 > num_u - 1) {
                indexiplus1 = 0;
            }

            int indexjplus1 = j+1;
            if (indexjplus1 > num_v - 1) {
                indexjplus1 = 0;
            }

            point vertex1 = interpolatedPoints[indexiplus1][j]; // 1
            int index1 = vertexIndexMapping.at(vertex1);

            point vertex2 = interpolatedPoints[i][indexjplus1]; // 2
            int index2 = vertexIndexMapping.at(vertex2);

            point vertex3 = interpolatedPoints[indexiplus1][indexjplus1]; // 3
            int index3 = vertexIndexMapping.at(vertex3);

            if(debug) {
                cout << "Indices at: " << index0 << ", " << index1 << ", " << index2 << ", " << index3 << endl;
            }

            vertexNormalMapping[vertex0] = calculateNormal(vertex0, vertex1, vertex2);
            vertexNormalMapping[vertex1] = calculateNormal(vertex1, vertex3, vertex2);
            vertexNormalMapping[vertex2] = calculateNormal(vertex2, vertex0, vertex1);
            vertexNormalMapping[vertex3] = calculateNormal(vertex3, vertex2, vertex1);

            coordIndexSetVals << index0 << ", ";
            coordIndexSetVals << index1 << ", ";
            coordIndexSetVals << index2 << ", ";
            coordIndexSetVals << -1 << ", " << endl;

            coordIndexSetVals << index1 << ", ";
            coordIndexSetVals << index3 << ", ";
            coordIndexSetVals << index2 << ", ";
            coordIndexSetVals << -1 << ", " << endl;

            numPts++;
        }
    }



    if(debug) {
        cout << "Serializing to OI format" << endl;
    }


    pointVals << "]" << std::endl;
    coordIndexSetVals << "]" << std::endl;

    Node* shapeHints = new Node("ShapeHints {","vertexOrdering        COUNTERCLOCKWISE\n","}");
    root->addChild(shapeHints);

    Node* interpolatedSeperator = new Node("Separator {","","}");

    Node* interpolatedPointsNode = new Node("Coordinate3 {",pointVals.str(),"}"); // SoCoordinate3

    Node* indexedFaceSet = new Node("IndexedFaceSet {",coordIndexSetVals.str(),"}"); // SoIndexedLineSet

    interpolatedSeperator->addChild(interpolatedPointsNode);

    if(shadeWithNormals) {
        // Calculate vertex normals as we are doing smooth shading.
        Node* normalBindingNode = new Node("NormalBinding {","value        PER_VERTEX_INDEXED","}");
        interpolatedSeperator->addChild(normalBindingNode);

        std::ostringstream normalVectorVals;

        normalVectorVals << "vector [" << std::endl;

        for(auto cit= vertices.begin(); cit != vertices.end(); cit++) {
            point p = *cit;
            point n = vertexNormalMapping[p];
            normalVectorVals << n.x << " " << n.y << " " << n.z << "," << std::endl;
        }

        normalVectorVals << "]" << std::endl;

        Node* normalVectorsNode = new Node("Normal {",normalVectorVals.str(),"}");
        interpolatedSeperator->addChild(normalVectorsNode);


    }

    interpolatedSeperator->addChild(indexedFaceSet);
    root->addChild(interpolatedSeperator);


    if(!debug) {
        string ivContent = root->getString();
        ostringstream outputContent;
        outputContent << "#Inventor V2.0 ascii" << endl << ivContent;

        std::cout << outputContent.str() << std::endl;
    }


    return 0;
}
