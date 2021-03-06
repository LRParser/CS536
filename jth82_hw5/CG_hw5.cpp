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
    double w;

    point(double x_i, double y_i, double z_i) {
        x = x_i;
        y = y_i;
        z = z_i;
        w = 1.0;
    }

    point() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        w = 1.0;
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
return std::make_tuple(a.x,a.y,a.z,a.w) <  std::make_tuple(b.x, b.y, b.z,b.w);
}



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



double to_radians(double degree) {
    return degree * (M_PI/180.0);
}


vector<point> drawQuad(point ll, point ur) {
    vector<point> returnVec;

    point p0(ur.x, ur.y, ur.z); // upper right
    point p1(ll.x, ur.y, ur.z);
    point p2(ll.x, ll.y, ur.z);
    point p3(ur.x, ll.y, ur.z);
    point p4(ur.x, ur.y, ll.z);
    point p5(ll.x, ur.y, ll.z);
    point p6(ll.x, ll.y, ll.z); // lower left
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

vector<vector<double>> getIdentityMatrix() {
    vector<vector<double>> result;

    for(int i = 0; i < 4; i++) {
        vector<double> row;
        for(int j = 0; j < 4; j++) {
            row.push_back(0.0);
        }
        result.push_back(row);
    }

    // Now make it an identity matrix
    result[0][0] = 1.0;
    result[1][1] = 1.0;
    result[2][2] = 1.0;
    result[3][3] = 1.0;

    return result;
}

vector<vector<double>> getYRotationMatrix(double thetaRadians) {

    vector<vector<double>> rotationMatrix = getIdentityMatrix();

    rotationMatrix[0][0] = cos(thetaRadians);
    rotationMatrix[0][2] = sin(thetaRadians);

    rotationMatrix[2][0] = -1*sin(thetaRadians);
    rotationMatrix[2][2] = cos(thetaRadians);

    return rotationMatrix;
}

vector<vector<double>> getZRotationMatrix(double thetaRadians) {

    vector<vector<double>> rotationMatrix = getIdentityMatrix();

    rotationMatrix[0][0] = cos(thetaRadians);
    rotationMatrix[0][1] = -1 * sin(thetaRadians);

    rotationMatrix[1][0] = sin(thetaRadians);
    rotationMatrix[1][1] = cos(thetaRadians);

    return rotationMatrix;
}

double dot(vector<double> row, vector<double> column) {
    double retVal = 0.0;
    for(int i = 0; i < 4; i++) {
        retVal += row[i]*column[i];
    }
    return retVal;
}

point matrixVectorMult(vector<vector<double>> matrix, point rhs) {
    point retval;

    retval.x = (matrix[0][0] * rhs.x) + (matrix[0][1] * rhs.y) + (matrix[0][2] * rhs.z) + (matrix[0][3] * rhs.w);
    retval.y = (matrix[1][0] * rhs.x) + (matrix[1][1] * rhs.y) + (matrix[1][2] * rhs.z) + (matrix[1][3] * rhs.w);
    retval.z = (matrix[2][0] * rhs.x) + (matrix[2][1] * rhs.y) + (matrix[2][2] * rhs.z) + (matrix[2][3] * rhs.w);
    retval.w = (matrix[3][0] * rhs.x) + (matrix[3][1] * rhs.y) + (matrix[3][2] * rhs.z) + (matrix[3][3] * rhs.w);

    return retval;
}

vector<vector<double>> matMult(vector<vector<double>> lhs, vector<vector<double>> rhs) {
    vector<vector<double>> result;

    for(int i = 0; i < 4; i++) {
        vector<double> row;
        for(int j = 0; j < 4; j++) {
            row.push_back(0.0);
        }
        result.push_back(row);
    }

    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            result[i][j] = dot(lhs[i],rhs[j]);
            // Find dot product of i row of lhs and j column of rhs
        }
    }
    return result;
}



vector<vector<double>> getTranslationMatrix(point p) {
    vector<vector<double>> translationMatrix = getIdentityMatrix();


    translationMatrix[0][3] = p.x;
    translationMatrix[1][3] = p.y;
    translationMatrix[2][3] = p.z;

    return translationMatrix;
}

vector<point> applyTransformationMatrixToPoints(vector<point> input, vector<vector<double>> matrix) {
    vector<point> retVal;
    for(auto it = input.begin(); it < input.end(); it++) {
        point p = *it;
        retVal.push_back(matrixVectorMult(matrix, p));
    }
    return retVal;
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
                std::cerr << "Must provide theta1 value after -t argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-u") {

            if (i + 1 < argc) {
                theta2 = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide theta2 value after -u argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-v") {

            if (i + 1 < argc) {
                theta3 = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide theta3 value after -v argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-l") {

            if (i + 1 < argc) {
                link1Length = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide link1Length value after -l argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-m") {

            if (i + 1 < argc) {
                link2Length = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide link2Length value after -m argument" << std::endl;
            }
        } else if (std::string(argv[i]) == "-n") {

            if (i + 1 < argc) {
                link3Length = std::stod(std::string(argv[i + 1]));
            } else {
                std::cerr << "Must provide link3Length value after -n argument" << std::endl;
            }
        }
        else if (std::string(argv[i]) == "-d") {
            debug = true;
        }
    }

    Node *root = new Node("", "", "");

    double baseHeight = 1.0;

    // Draw base
    point baseLL(-2,-2,0);
    point baseUR(2,2,1);
    vector<point> baseVertices = drawQuad(baseLL, baseUR);
    Node* baseSeparator = constructQuadNode(baseVertices);
    root->addChild(baseSeparator);

    // Draw link 1
    float theta1Radians = to_radians(theta1);
    point l1LL(-.5,-.5,0);
    point l1UR(.5,.5,link1Length);
    vector<point> link1Vertices = drawQuad(l1LL, l1UR);
    point link1Translation;
    link1Translation.z = baseHeight;
    vector<vector<double>> translationMatrixLink1 = getTranslationMatrix(link1Translation);
    vector<vector<double>> rotationMatrixLink1 = getZRotationMatrix(theta1Radians); //  getIdentityMatrix(); // Update later

    vector<vector<double>> transformationMatrixLink1 = matMult(translationMatrixLink1,rotationMatrixLink1); // matMult(rotationMatrixLink1,translationMatrixLink1); //
    link1Vertices = applyTransformationMatrixToPoints(link1Vertices, transformationMatrixLink1);
    Node* link1Separator = constructQuadNode(link1Vertices);
    root->addChild(link1Separator);

    // Draw link 2

    double theta2Radians = to_radians(theta2);

    // Now to apply a rotation, create rotation matrix
    point translationValsLink2;
    translationValsLink2.z = link1Length;
    vector<vector<double>> translationMatrixLink2 = getTranslationMatrix(translationValsLink2);
    vector<vector<double>> rotationMatrixLink2 = getYRotationMatrix(theta2Radians);

    vector<vector<double>> link2TransformationMatrix = matMult(transformationMatrixLink1,matMult(rotationMatrixLink2,translationMatrixLink2));

    point l2LL(-.5,-.5,0);
    point l2UR(.5,.5,link2Length);

    vector<point> link2Vertices = drawQuad(l2LL, l2UR);

    vector<point> link2VerticesRotatedTranslated = applyTransformationMatrixToPoints(link2Vertices, link2TransformationMatrix);

    Node* link2Separator = constructQuadNode(link2VerticesRotatedTranslated);
    root->addChild(link2Separator);

    // Draw link 3

    double theta3Radians = to_radians(theta3);
    point translationValsLink3;
    translationValsLink3.z = link2Length;
    vector<vector<double>> translationMatrixLink3 = getTranslationMatrix(translationValsLink3);
    vector<vector<double>> rotationMatrixLink3 = getYRotationMatrix(theta3Radians);

    vector<vector<double>> link3TransformationMatrix = matMult(link2TransformationMatrix,matMult(rotationMatrixLink3,translationMatrixLink3));

    point l3LL(-.5,-.5,0);
    point l3UR(.5,.5,link3Length);

    vector<point> link3Vertices = drawQuad(l3LL, l3UR);

    vector<point> link3VerticesRotatedTranslated = applyTransformationMatrixToPoints(link3Vertices, link3TransformationMatrix);


    Node* link3Separator = constructQuadNode(link3VerticesRotatedTranslated);
    root->addChild(link3Separator);


    // Now we need to place a sphere at the origin and at the end of the fourth link
    // Both have radius of 0.2
    std::ostringstream radiusStr;
    radiusStr << "radius  " << 0.2 << std::endl;
    Node* sphereNode = new Node("Sphere {",radiusStr.str(),"}");

    // Sphere at origin
    Node* originPointSep = new Node("Separator {","","}");
    point originPoint;
    std::ostringstream transformStrOrigin;
    transformStrOrigin << "translation " << originPoint.x << " " << originPoint.y << " " << originPoint.z << std::endl;
    Node* originTranslation = new Node("Transform {",transformStrOrigin.str(),"}");
    originPointSep->addChild(originTranslation);
    originPointSep->addChild(sphereNode);
    root->addChild(originPointSep);

    // Sphere at end of final link
    Node* finalPointSep = new Node("Separator {","","}");

    // Need a translation to the top of the final link

    vector<vector<double>> finalPointTransform = matMult(link3TransformationMatrix,getTranslationMatrix(point(0,0,0)));
    point finalPoint = matrixVectorMult(finalPointTransform,point(0,0,link3Length));
    // Translate it by z

    std::ostringstream transformStrFinal;
    transformStrFinal << "translation " << finalPoint.x << " " << finalPoint.y << " " << finalPoint.z << std::endl;
    Node* translationFinal = new Node("Transform {",transformStrFinal.str(),"}");
    finalPointSep->addChild(translationFinal);
    finalPointSep->addChild(sphereNode);
    root->addChild(finalPointSep);



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
