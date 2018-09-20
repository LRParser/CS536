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
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoSeparator.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

static char * buffer;
static size_t buffer_size = 0;

static void *
buffer_realloc(void * bufptr, size_t size)
{
  buffer = (char *)realloc(bufptr, size);
  buffer_size = size;
  return buffer;
}

static SbString
buffer_writeaction(SoNode * root)
{
  SoOutput out;
  buffer = (char *)malloc(1024);
  buffer_size = 1024;
  out.setBuffer(buffer, buffer_size, buffer_realloc);

  SoWriteAction wa(&out);
  wa.getOutput()->openFile( "output.iv" );
  //wa.getOutput()->setBinary( FALSE );  // Optional: write binary format

  wa.apply(root);
  SbString s(buffer);
  wa.getOutput()->closeFile();

  free(buffer);
  return s;
}

int
main(int argc, char ** argv)
{
    string fName;
    vector<int> xs;
    vector<int> ys;
    vector<int> zs;
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
        int x, y, z;
        std::istringstream ss(currentLine);
        if (!(ss >> x >> y >> z)) {
            // std::cerr << "Cannot parse" << std::endl;
        }
        xs.push_back(x);
        ys.push_back(y);
        zs.push_back(z);
        std::cout << "Parsed: " << x << " " << y << " " << z << std::endl;
    }

    std::cout << "Done" << std::endl;

    /*
  SoDB::init();

  SoSeparator * root = new SoSeparator;
  root->ref();

  root->addChild(new SoSphere);

  SbString s = buffer_writeaction(root);
  (void)fprintf(stdout, "%s\n", s.getString());

  root->unref();
     */
  return 0;
}
