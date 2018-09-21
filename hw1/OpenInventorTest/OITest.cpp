/*
 * OITest.c
 *
 *  Created on: Sep 15, 2018
 *      Author: joe
 */
#include <Inventor/SoOutput.h>
#include <Inventor/SoDB.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/nodes/SoCone.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoIndexedLineSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/fields/SoMFVec3f.h>
#include <Inventor/nodes/SoCoordinate3.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

using namespace std;

// From OpenInventor library

static char * buffer;
static size_t buffer_size = 0;

static void * buffer_realloc(void * bufptr, size_t size)
{
  buffer = (char *)realloc(bufptr, size);
  buffer_size = size;
  return buffer;
}



struct point {
    float x;
    float y;
    float z;
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


float fact(int i) {
    if (i == 0) {
        return 1;
    }
    else {
        return i * fact(i-1);
    }
    /*
    auto iter = factorialMap.find(i);
    if (iter != factorialMap.end()) {
       return iter->second;
    }
    else {
        float retVal = i * fact(i-1);
        factorialMap.insert(std::pair<int,int>(i,retVal));
        return retVal;
    }
     */
}



float nchoosei(float n, float i) {
    return fact(n) / (fact(i)*fact(n-i));
}


int
main(int argc, char ** argv)
{
    string fName;
    vector<point> points;

    for(int i=0; i < argc; i++) {
        if (std::string(argv[i]) == "-f") {

            if (i + 1 < argc) {
                fName = std::string(argv[i]);
            }
            else {
                std::cerr << "Must provide file name after -f argument" << std::endl;
            }

        }
    }

    if (fName.empty()) {
        fName = "cpts_in.txt";
    }

    std::cout << "Reading from file: " << fName << std::endl;
    std::ifstream input(fName.c_str());
    if (input.fail()) {
        std::cerr << "Failed to open" << std::endl;
    }
    string currentLine;
    while(std::getline(input, currentLine)) {
        std::cout << currentLine << std::endl;
        point point1;
        std::istringstream ss(currentLine);
        if (!(ss >> point1.x >> point1.y >> point1.z)) {
            // std::cerr << "Cannot parse" << std::endl;
        }
        points.push_back(point1);

        std::cout << "Parsed: " << point1.x << " " << point1.y << " " << point1.z << std::endl;
    }

    std::cout << "Done" << std::endl;


    SoDB::init();

    SoSeparator* root = new SoSeparator;
    root->ref();

    std::cout << "Start" << std::endl;

    for (auto it = points.begin(); it != points.end(); it++) {

        SoSeparator* cpSep = new SoSeparator();
        root->addChild(cpSep);

        float x = it->x;
        float y = it->y;
        float z = it->z;
        std::cout << "Placing: " << x << " " << y << " " << z << std::endl;


        SoTranslation* translation = new SoTranslation;
        translation->translation.setValue(x,y,z);
        cpSep->addChild(translation);
        SoSphere* sphere = new SoSphere();
        sphere->radius = .1;
        cpSep->addChild(sphere);
    }

    SoSeparator* interpolatedSeperator = new SoSeparator;

    SoCoordinate3* interpolatedPoints = new SoCoordinate3();
    SoIndexedLineSet* indexedLineSet = new SoIndexedLineSet;


    float du = 0.05;
    float u = 0.0;
    int n = (int) points.size();

    vector<point> calcPoints;
    while(u <= 1.0) {
        point currentPoint;
        currentPoint.x = 0.0;
        currentPoint.y = 0.0;
        currentPoint.z = 0.0;

        for(int i = 0; i < n; i++) {
            point controlPoint = points.at(i);

            double factor = nchoosei(n,i) * pow(1-u,n-i) * pow(u,i);
            point calcPoint = pointMult(factor,controlPoint);
            currentPoint = pointAdd(currentPoint,calcPoint);

        }
        calcPoints.push_back(currentPoint);

        u += du;
    }


    int i = 0;
    for(auto it = calcPoints.begin(); it != calcPoints.end(); it++) {
        interpolatedPoints->point.set1Value(i,it->x,it->y,it->z);
        indexedLineSet->coordIndex.set1Value(i,i);

        i++;
    }
    indexedLineSet->coordIndex.set1Value(i,-1);
    

    interpolatedSeperator->addChild(interpolatedPoints);


    interpolatedSeperator->addChild(indexedLineSet);

    root->addChild(interpolatedSeperator);



    SoOutput out;
    buffer = (char *)malloc(1024);
    buffer_size = 1024;
    out.setBuffer(buffer, buffer_size, buffer_realloc);

    SoWriteAction wa(&out);
    wa.getOutput()->openFile( "output.iv" );
    //wa.getOutput()->setBinary( FALSE );  // Optional: write binary format

    wa.apply(root);
    wa.getOutput()->closeFile();

    free(buffer);

    root->unref();
   
    return 0;
}
