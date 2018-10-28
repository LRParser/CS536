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

double magnitude(point p) {
    return sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
}

double dot(point p, point v) {
    return p.x*v.x+p.y*v.y+p.z*v.z;
}

double cosineoftheta(point u, point v) {
    return dot(u,v)/(magnitude(u)*magnitude(v));
}

double angleoftheta(point u, point v) {
    return acos(cosineoftheta(u,v));
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
        else if (std::string(argv[i]) == "-d") {
            debug = true;
        }
    }

    if(debug) {
        cout << "num_u is: " << num_u;
        cout << "num_v is: " << num_v;

    }

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

    // We always interpolate at 0 and 1
    double du = 1 / (double(num_u) - 1);
    double dv = 1 / (double(num_v) - 1);

    if(debug) {
        cout << "du is: " << du << " dv is: " << dv << endl;
    }

    double u = 0.0;
    double v = 0.0;
    int numPts = 0;
    int numUPts = 0;
    int numVPts = 0;
    vector<vector<point>> interpolatedPoints;

    while(numUPts <= num_u) {
        vector<point> uVec;

        v = 0.0;
        numVPts = 0;

        while(numVPts <= num_v) {
            point p;
            uVec.push_back(p);
            numVPts++;
        }
        numUPts++;
        interpolatedPoints.push_back(uVec);

    }

    if(debug) {
        cout << "Created placeholder for " << numPts << " points" << endl;
        cout << "Size is: " << numUPts << " by " << numVPts << endl;
    }

    std::ostringstream pointVals;
    std::ostringstream coordIndexSetVals;

    pointVals << "point [" << std::endl;


    int n = 3;
    int m = 3;




    int uIdx = 0;
    int vIdx = 0;
    u = 0.0;
    v = 0.0;

    while(u <= 1 + 1e-9) {

        while(v <= 1) {

            if(debug) {
            }
            if(debug) {
                cout << "u is: " << u << " ";
                cout << "v is: " << v << " ";
            }

            point currentPoint;
            currentPoint.x = 0.0;
            currentPoint.y = 0.0;
            currentPoint.z = 0.0;

            for(int i = 0; i <= n; i++) {

                double factorU = binomial(n, i) * pow(u,i) * pow(1-u,n-i) ;

                for(int j=0; j <= m; j++) {

                    double factorV = binomial(m, j) * pow(v,j) * pow(1-v,m-j);
                    if(debug) {
                        // cout << " Get control point at: " << j << "," << i << endl;
                    }
                    point controlPoint = k[i][j]; // TBD
                    currentPoint = currentPoint + (controlPoint * factorU * factorV);
                }
            }

            if(debug) {
                cout << " Set interpolated point at: " << uIdx << "," << vIdx << " to x: " << currentPoint.x << ", y: " << currentPoint.y << ", z: " << currentPoint.z << endl;
            }


            interpolatedPoints[uIdx][vIdx] = currentPoint;

            vIdx++;
            v = v+dv;
        }

        v = 0.0;
        vIdx = 0;
        uIdx++;
        u=u+du;

    }

    if(debug) {
        cout << "numPts is" << numPts << endl;
        cout << "Writing to OpenInventor format" << endl;
    }


    coordIndexSetVals << "coordIndex [" << std::endl;


    int vertexIndex = 0;
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

    for(int i = 0; i < num_u - 1; i++) {
        for(int j = 0; j < num_v -1; j++) {

            if(debug) {
                cout << "Tesselate for u = " << i << " and v = " << j << endl;
            }

            // Four distinct points become a patch (two tesselated triangles)
            point vertex0 =  interpolatedPoints[i][j]; // 0
            int index0 = vertexIndexMapping.at(vertex0);

            point vertex1 = interpolatedPoints[i+1][j]; // 1
            int index1 = vertexIndexMapping.at(vertex1);

            point vertex2 = interpolatedPoints[i][j+1]; // 2
            int index2 = vertexIndexMapping.at(vertex2);

            point vertex3 = interpolatedPoints[i+1][j+1]; // 3
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



    // Plot the control points

    for(int i = 0; i <= 3; i++) {

        for(int j = 0; j <=3; j++) {

            Node* cpSep = new Node("Separator {","","}");

            point point1 = k.at(j).at(i);
            double x = point1.x;
            double y = point1.y;
            double z = point1.z;

            if(debug) {
                // std::cout << "Placing: " << x << " " << y << " " << z << std::endl;
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





    if(!debug) {
        string ivContent = root->getString();
        ostringstream outputContent;
        outputContent << "#Inventor V2.0 ascii" << endl << ivContent;

        std::cout << outputContent.str() << std::endl;
    }


    return 0;
}
