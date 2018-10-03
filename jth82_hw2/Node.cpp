//
// Created by joe on 9/26/18.
//

#include "Node.h"

Node::Node(std::string pre, std::string con, std::string post) {
    prefix = pre;
    content = con;
    postfix = post;
}

void Node::addChild(Node* n) {
    children.push_back(n);
}

std::string Node::getString() {
    std::ostringstream output;
    output << prefix << std::endl;
    for(auto it = children.begin(); it != children.end(); it++) {
        Node* node = *it;
        output << node->getString() << std::endl;
    }
    output << content << postfix << std::endl;
    return output.str();
}