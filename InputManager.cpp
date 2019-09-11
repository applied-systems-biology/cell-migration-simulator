/* 
 * File:   main.cpp
BSD 3-Clause License

Copyright (c) 2018-2019, Marco Blickensdorf
Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge  
https://www.leibniz-hki.de/en/applied-systems-biology.html  
HKI-Center for Systems Biology of Infection  
Leibniz Institute for Natural Product Research and Infection Biology -  
Hans Knöll Insitute (HKI)  
Adolf-Reichwein-Straße 23, 07745 Jena, Germany  

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Created on July 4, 2018, 3:45 PM
 */
#include "libs/pugixml/pugixml.hpp"
#include "libs/pugixml/pugixml.cpp"
#include "InputManager.h"
#include <map>
#include <string>
#include <iostream>
#include "Framework.h"

pugi::xml_document config;

/*
 Constructor is reading the input file and saving it to config.
 */
InputManager::InputManager() {


    
}

InputManager::InputManager(const InputManager& orig) {
}

InputManager::~InputManager() {
}


void InputManager::init(std::vector<float> args, std::string runname){
    pugi::xml_parse_result configXML = config.load_file("config.xml");
    std::cout << "Input file was read: " << configXML.description() << std::endl;
//    plog::init(plog::warning, "log.csv", 1000000, 5); 
//    for(int i = 0; i < args.size();i++){
//        std::cout << args[i] << runname << std::endl;
//    }
    params = args;
    runtitle = runname;
    
}

int InputManager::getIntParameter(std::string field, std::string identifier){
    for(pugi::xml_node_iterator it = config.children().begin(); it != config.children().end();it++){
        std::string test = (*it).name();    
        if(strcmp(test.c_str(),field.c_str())==0){
            //the right field was found
            std::string type = (*it).child(identifier.c_str()).attribute("type").as_string();
            if(type.compare("")==0){
                std::cout << "Parameter was not found: " << field << "|" << identifier << std::endl;
                LOG_ERROR << "Parameter was not found: " << field << "|" << identifier ;
            }else{
                if(type.compare("int")!=0){
                    std::cout << "Warning!You are reading a "<< type <<"-type parameter as int: " << identifier << std::endl;          //put some better log-management in here!
                    LOG_WARNING << "Warning!You are reading a "<< type <<"-type parameter as int: " << identifier;          //put some better log-management in here!
                }
                return (*it).child(identifier.c_str()).attribute("value").as_int();
            }
        }
    }    
    int a = 0;
}


//some parameters can be overwirtten by command line arguments. this is realized here
float InputManager::getFloatParameter(std::string field, std::string identifier){
    
    if(params.size()>0){
        if(identifier.compare("neighbour_weight")==0){
            LOG_WARNING << "I am overwriting the parameter neighbour_weight with value " << params[0];
            std::cout << "I am overwriting the parameter neighbour_weight with value " << params[0] << std::endl;

            return params[0];
        }
        if(identifier.compare("pos_weight")==0){
            LOG_WARNING << "I am overwriting the parameter pos_weight with value " << params[1];
            std::cout << "I am overwriting the parameter pos_weight with value " << params[1]<< std::endl;
            return params[1];
        }
        if(identifier.compare("distance_weight")==0){
            LOG_WARNING << "I am overwriting the parameter distance_weight with value " << params[2];
            std::cout << "I am overwriting the parameter distance_weight with value " << params[2]<< std::endl;
            return params[2];
        }
        if(identifier.compare("resolution")==0){
            LOG_WARNING << "I am overwriting the parameter resolution with value " << params[3];
            std::cout << "I am overwriting the parameter resolution with value " << params[3]<< std::endl;
            return params[3];
        }
        if(identifier.compare("frontback")==0){
            LOG_WARNING << "I am overwriting the parameter frontback with value " << params[4];
            std::cout << "I am overwriting the parameter frontback with value " << params[4]<< std::endl;
            return params[4];
        }
    }
    
    for(pugi::xml_node_iterator it = config.children().begin(); it != config.children().end();it++){
        std::string test = (*it).name();        
        if(strcmp(test.c_str(),field.c_str())==0){
            //the right field was found
            std::string type = (*it).child(identifier.c_str()).attribute("type").as_string();
            if(type.compare("")==0){
                std::cout << "Parameter was not found: " << field << "|" << identifier << std::endl;
                LOG_ERROR << "Parameter was not found: " << field << "|" << identifier ;
            }else{
                if(type.compare("float")!=0){
                    std::cout << "Warning!You are reading a "<< type <<"-type parameter as float!" << std::endl;          //put some better log-management in here!
                    LOG_WARNING << "Warning!You are reading a "<< type <<"-type parameter as float!" << std::endl;          //put some better log-management in here!
                }
                return (*it).child(identifier.c_str()).attribute("value").as_float();
            }
        }
    }    
    int a = 0;
}

bool InputManager::getBoolParameter(std::string field, std::string identifier){
    for(pugi::xml_node_iterator it = config.children().begin(); it != config.children().end();it++){
        std::string test = (*it).name();        
        if(strcmp(test.c_str(),field.c_str())==0){
            //the right field was found
            std::string type = (*it).child(identifier.c_str()).attribute("type").as_string();
            if(type.compare("")==0){
                std::cout << "Parameter was not found: " << field << "|" << identifier << std::endl;
                LOG_ERROR << "Parameter was not found: " << field << "|" << identifier ;
            }else{
                if(type.compare("bool")!=0){
                    std::cout << "Warning!You are reading a "<< type <<"-type parameter as bool!" << std::endl;          //put some better log-management in here!
                    LOG_WARNING << "Warning!You are reading a "<< type <<"-type parameter as bool!" << std::endl;          //put some better log-management in here!
                }
                return (*it).child(identifier.c_str()).attribute("value").as_bool();
            }
        }
    }    
    int a = 0;
}

std::string InputManager::getStringParameter(std::string field, std::string identifier){
    if(runtitle.compare("")!=0){
        if(identifier.compare("runname")==0){
            LOG_WARNING << "I am overwriting the parameter runname with value " << runtitle;
            return runtitle;
        }
    }
    
    
    for(pugi::xml_node_iterator it = config.children().begin(); it != config.children().end();it++){
        std::string test = (*it).name();        
        if(strcmp(test.c_str(),field.c_str())==0){
            //the right field was found
            std::string type = (*it).child(identifier.c_str()).attribute("type").as_string();
            if(type.compare("")==0){
                std::cout << "Parameter was not found: " << field << "|" << identifier << std::endl;
                LOG_ERROR << "Parameter was not found: " << field << "|" << identifier ;
            }else{
                if(type.compare("string")!=0){
                    std::cout << "Warning!You are reading a "<< type <<"-type parameter as string!" << std::endl;          //put some better log-management in here!
                    LOG_WARNING << "Warning!You are reading a "<< type <<"-type parameter as string!" << std::endl;          //put some better log-management in here!
                }
                return (*it).child(identifier.c_str()).attribute("value").as_string();
            }

        }
    }    
}

