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
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

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
        sphere->radius = .05;
        cpSep->addChild(sphere);
    }

    SoSeparator* interpolatedSeperator = new SoSeparator;

    SoCoordinate3* interpolatedPoints = new SoCoordinate3();
    for(int i = 0; i < 5; i++) {
        interpolatedPoints->point.set1Value(i,i,i,i);
    }

    interpolatedSeperator->addChild(interpolatedPoints);

    SoIndexedLineSet* indexedLineSet = new SoIndexedLineSet;
    indexedLineSet->coordIndex.set1Value(0,0);
    indexedLineSet->coordIndex.set1Value(1,1);
    indexedLineSet->coordIndex.set1Value(2,2);
    indexedLineSet->coordIndex.set1Value(3,3);
    indexedLineSet->coordIndex.set1Value(4,4);
    indexedLineSet->coordIndex.set1Value(5,-1);

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
