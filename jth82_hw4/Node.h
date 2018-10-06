//
// Created by joe on 9/26/18.
//

#ifndef CG_HW1_NODE_H
#define CG_HW1_NODE_H
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

class Node {

public:
    std::string prefix;
    std::string content;
    std::string postfix;
    std::vector<Node*> children;

    Node(std::string pre, std::string con, std::string post);
    void addChild(Node* n);
    std::string getString();

};


#endif //CG_HW1_NODE_H
