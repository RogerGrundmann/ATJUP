#include <iostream>
#include <stdexcept>

#include "cJupiterModel.h"
#include "Config.h"
#include "tinyxml2.h"

void cJupiterModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    cout << "reading config-xml-file in cJupiterModel::LoadConfig ======= error message = err " << err << "\n\n";
    if(err){
        doc.PrintError();
        throw std::invalid_argument("couldn't load config file");
    }
    XMLElement* atjup = doc.FirstChildElement("atjup");
    if(!atjup){
        return;
    }
    XMLElement* elem_common = doc.FirstChildElement("atjup")->FirstChildElement("elem_common");
    if(!elem_common){
        return;
    }
    XMLElement* elem_atmosphere = doc.FirstChildElement("elem_common")->FirstChildElement("elem_jupiter");
    if(!elem_jupiter){
        return;
    }
    #include "JupiterLoadConfig.cpp.inc"
}
