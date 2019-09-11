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
#include "Framework.h"
#include <iostream>
#include "boost/filesystem.hpp"
#include "Environment.h"
#include <stdlib.h>     /* srand, rand */

InputManager Framework::input;
Environment  Framework::world;
boost::random::mt19937 Framework::rng;

Framework::Framework() {
}

Framework::Framework(const Framework& orig) {
}

Framework::~Framework() {
}

void Framework::init(std::vector<float> args, std::string runname){
    //init logger
    plog::init(plog::verbose, "log.csv", 10000000, 5);
    LOG_INFO << "Logger was initialized";
    
    //init random seed
    rng.seed(time(NULL));

    iteration = 0;  //the class internal clock/time

    
    
    //init input manager
    LOG_INFO << "Input manager was initialized";
    input.init(args,runname);
    
    //init time
    time_it = 0;            //time_it corresponds to model iterations and thus is not directly related to a real time scale
    
    //init the environment
    world.init();
    
    //init the output manager
    initOutput();
    LOG_INFO << "Output was initialized"; 

}

void Framework::initOutput(){
    //bafore using the output I should prepare the output
    /*
     create folder according to the run
     * save the information on the run
     * save the config file
     * create a dir to save the cells to
     */
    //create results folder if not existent
    boost::filesystem::path dir("results");
    if (boost::filesystem::create_directory(dir)){
        std::cout << "results folder was created" << "\n";
        LOG_INFO << "Results folder was created.";
    }
    ////create a folder for one specific run which name consist of the date and the  run_name (e.g. debug, test, runXYZ))
    //runname
    std::string runname = input.getStringParameter("outputSettings","runname");
    //time to str
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%F_%T",timeinfo);
    std::string str(buffer);
    //create dir
    float res = input.getFloatParameter("modelSettings","resolution");
    float nb = input.getFloatParameter("modelSettings","pos_weight");
    float pw = input.getFloatParameter("modelSettings","neighbour_weight");
    float dw = input.getFloatParameter("modelSettings","distance_weight");
    float fb = input.getFloatParameter("modelSettings","frontback");
    std::ostringstream sss;
      /* initialize random seed: */
    srand (time(NULL));

    /* generate secret number between 1 and 10: */
    int randNum = rand() % 90000000 + 10000000; //to create a unique file name
    sss <<randNum << "_" << nb << "_" << pw << "_"<< dw << "_" << res << "_" << fb;
    std::cout << "nw:"<< nb << "_pw:" << pw << "_dw:"<< dw << "_res:" << res << "_fb:" << fb << std::endl;
    
    outpath = "results/"+str+"_"+runname+"_"+sss.str();
    boost::filesystem::path dir2(outpath);
    bool createFolder = boost::filesystem::create_directory(dir2);
    int counter = 1;
    while(!createFolder){
        std::stringstream ss;
        ss << counter;
        outpath = "results/"+str+"_"+runname+"_"+sss.str() + "_R" +  ss.str();
        boost::filesystem::path dir2(outpath);
        createFolder = boost::filesystem::create_directory(dir2);
        counter++;
    }
    std::cout << "results folder was created: " << outpath << "\n";
    LOG_INFO << "Results folder was created." << outpath;
    boost::filesystem::path dir3(outpath+"/cells");
    if (boost::filesystem::create_directory(dir3)){
        std::cout << "result cells folder was created" << "\n";
        LOG_INFO << "Result cells folder was created.";
    }
    
    ////copy config file to this directory
    boost::filesystem::copy_file("config.xml",outpath+"/config.xml");
}

void Framework::output(std::string path, int time_it){
    //i now to output the whole world, or other, depending on output options
    world.output(path,time_it);
}


bool Framework::run(){
    //here now the mode should be iterated
//    int limit = input.getIntParameter("modelSettings","iteration_limit");
    int limit2 = 100000;//world.getSize()*5;  
    int intervall = input.getIntParameter("outputSettings","intervall");
//    int intervall = world.getSize()/6;

    for(; iteration < limit2; iteration++){
        if(iteration%intervall==0){
            output(outpath,iteration);        
        }
//        std::cout << "world:" << world.getSize() << "at iteration " << iteration << std::endl;
//        if(world.getSize() != 4139){
//            std::cout << "NOOOOOOOO" << std::endl;
//            return false;
//            
//            std::cout << "NOOOOOOOO" << std::endl;
//
//        }
        int succIteration = world.iterate(iteration);
        if(succIteration != world.getN_cells()){
            LOG_FATAL << "Could not iterate all cells!" ;
            std::cout  << "Could not iterate all cells!" << std::endl;

        }
//        world.testIntegrity();
    }
    return true;
    LOG_INFO << "simulation was finished successfull";
}

void Framework::finalize(){
    LOG_INFO << "starting finalization.";
    //write down run statistics;
    world.doRunAnalysis(outpath);
    world.doShapeAnalysis(outpath);
        boost::filesystem::copy_file("log.csv",outpath+"/log.csv");

}




