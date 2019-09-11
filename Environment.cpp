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
#include "Environment.h"
#include "InputManager.h"
#include "Framework.h"
#include <iostream>
#include "boost/date_time/posix_time/posix_time.hpp" 
//#include <vector>
 
typedef boost::posix_time::ptime Time;
typedef boost::posix_time::time_duration TimeDuration;
 
namespace io = boost::iostreams;



Environment::Environment() {
}

Environment::Environment(const Environment& orig) {
}

Environment::~Environment() {
}

/*
 i once deleted this whole cpp file. so i try to add it now
 */

/*
 init the cells and tissue
 */
void Environment::init(){
    //load parameter
    InputManager im = Framework::getInputManager();
    int ncells = im.getIntParameter("modelSettings","n_c");
    float border = im.getFloatParameter("modelSettings","space");
    float resolution = im.getFloatParameter("modelSettings","resolution");  
    std::string tissueInclude = im.getStringParameter("modelSettings","tissue");  
    std::cout << "starting with tissue " << tissueInclude << std::endl;
    
    
    //a surrounding sphere with radius 15
//        int size_circleOut = 20/resolution;
//        int size_circleIn = (20-2*resolution)/resolution;
//
//        for(int x = -size_circleOut; x < size_circleOut;x++){
//            for(int y = -size_circleOut; y < size_circleOut;y++){
//                for(int z = -size_circleOut; z < size_circleOut;z++){
//                    float dist = sqrt(x*x+y*y+z*z);
//                    if(dist < size_circleOut && dist > size_circleIn){
//                        Voxel* a = new Voxel(x,y,z,-1);
//    //                    std::cout << x <<"|"<< z <<"|"<< y << std::endl;
//                        bool add = addVoxel(a);
//                        if(add) tissue.push_back(a);
//                    }
//                }
//            }
//        }


    
    

    

    
    
    //init cells
    for(int i = 0; i < ncells; i++){
        Cell* c1 = new Cell(i);
        //draw a random position and init the cell until it does return true
        bool cellPositioned = false;
        while(!cellPositioned){
            boost::random::uniform_real_distribution<> pf_space(-border,border);

            float x = 0;//pf_space(Framework::rng);
            float y = 0;//pf_space(Framework::rng);
            float z = 0;//pf_space(Framework::rng);
            cellPositioned = c1->init(x,y,z);
        }
        cells.push_back(c1);
        std::cout << "was initing cell " << i << " with " << c1->getContent().size() << "Voxel"   << std::endl;
    }
    //init the neighbour count numbers correct for each voxel of the system
    std::map<boost::tuple<int,int,int>,Voxel*>::iterator it;
    for(it= space.begin(); it != space.end();it++){
        int nb = getNeighbours((*it).second, (*it).second->getID());
        (*it).second->setNB(nb);    
//        std::cout << "setting initial nb count to " << (*it).second->getNB() << std::endl;   
    }
    
    
    //randomly positioned spheres in space as obstacling tissue/cells
    if(tissueInclude.compare("spheres")==0){
        int nObs = 500; // number of obstacle spheres
        float rObs = 3/resolution; // radius of the obstacle sphere
        for(int i = 0; i < nObs;i++){
            //create a random position to place the sphere
            boost::random::uniform_real_distribution<> pf_space(-border/resolution,border/resolution);
            float x = pf_space(Framework::rng);
            float y = pf_space(Framework::rng);
            float z = pf_space(Framework::rng);
            for(int xL = -rObs; xL < rObs;xL++){
                for(int yL = -rObs; yL < rObs;yL++){
                    for(int zL = -rObs; zL < rObs;zL++){
                        float dist = sqrt(xL*xL+yL*yL+zL*zL);
                        if(dist < rObs && dist > 0.7*rObs){
                            Voxel* a = new Voxel(xL+x,yL+y,zL+z,-1);
    //                        std::cout << xL+x <<"|"<< yL+y <<"|"<< zL+z << std::endl;
                            bool add = addVoxel(a);
                            if(add) tissue.push_back(a);
                        }
                    }
                }
            }
        }
    }

    


    

        //tubular obstacles for celles as from e.g. veines
    if(tissueInclude.compare("tubule")==0){
        float tradius = 12/resolution;
        float tubeRadius = tradius + 5*resolution;
        float innerTubeRadius = tradius;
        boost::random::uniform_real_distribution<> Uni_ZeroOne(0,1);
        int middleX = 0;
        int middleY = 0;
        //assisting part to create random walk at center of tubule
        std::vector<int> shiftsX;
        std::vector<int> shiftsY;  
        for(int zL = -border/resolution; zL < border/resolution; zL++){     //this loop is creating the random walk
            float a = Uni_ZeroOne(Framework::rng);
            float b = Uni_ZeroOne(Framework::rng);
            if(a>0.5){shiftsX.push_back(1);} else{ shiftsX.push_back(-1);}
            if(b>0.5){shiftsY.push_back(1);} else{ shiftsY.push_back(-1);}
        }
        int counter = 0;
        for(int zL = -border/resolution; zL <1; zL++){     //this loop is getting the position of the random walk at z=0 (wher cell is initialized)
            middleX += shiftsX.at(counter);
            middleY += shiftsY.at(counter);
            counter++;
        }    
    //    std::cout << middleX << "|" << middleX << std::endl;
        counter = 0;
        int shiftX = 0;
        int shiftY = 0;
        for(int zL = -border/resolution; zL < border/resolution; zL++){
            for(int xL = -tubeRadius/resolution; xL < tubeRadius/resolution;xL++){
                for(int yL = -tubeRadius/resolution; yL < tubeRadius/resolution;yL++){
                    float dist = sqrt(yL*yL + xL * xL);
    //                std::cout << dist << std::endl;
    //                std::cout << xL <<"|"<< yL <<"|"<< zL << std::endl;
                    if(dist <= tubeRadius & dist > innerTubeRadius){
                        //insert the voxel in here
                        Voxel* a = new Voxel(xL+shiftX-middleX ,yL+shiftY-middleY,zL,-1);
    //                    std::cout << shiftsX.at(counter)-middleX <<"|"<< shiftsY.at(counter)-middleY <<"|"<< zL << ";:;:;" << dist << std::endl;
                        bool add = addVoxel(a);
                        if(add) tissue.push_back(a);
                    }
                }

            }
            shiftX+=shiftsX.at(counter);    //update shift variables
            shiftY+=shiftsY.at(counter);
            counter++;
    //        std::cout << a << ", " << b << std::endl;
        }   
    }


    
    
    //init the cell front data
    
    //a check if the neighbour comupations did really work!
//    Cell* a = cells.at(0);
//    bool works = true;
//    std::vector<Voxel*> cont = a->getContent();
//    for(int i =0; i < cont.size();i++){
//        boost::tuple<int,int,int> index(cont.at(i)->getX(),cont.at(i)->getY(),cont.at(i)->getZ());
//        int nb = getNeighbours(cont.at(i));
////        std::cout << "computed: " << nb << " vs. " << cont.at(i)->getNB() << " and world:" << space.at(index)->getNB() << std::endl;
//        if(nb != space.at(index)->getNB()){
//            works = false;
//            
//        }
//    }
//    if(!works){
//        std::cout << "Neighbourcount init did not work."  << std::endl;
//        LOG_FATAL << "Neighbourcount init did not work.";
//    }else{
//        std::cout << "Neighbourhood data was initialized correct"  << std::endl;
//        LOG_INFO << "Neighbourhood data was initialized correct";        
//    }
    for(int i = 0; i < cells.size();i++){
        cells.at(i)->initFront();
    }
   
};

bool Environment::addVoxel(Voxel* vox){
    //add cells
    if(vox->getX() == 0 && vox->getY()== 0 && vox->getZ() == 0){
        LOG_ERROR << "Fatal Error. You are trying to insert a voxel at occupied position";
        
        LOG_ERROR << "Fatal Error. You are trying to insert a voxel at occupied position";
    
    }
    
    
    if(voxelExists(vox)){
//        std::cout << "Voxel already exists:" << vox->getX() << "|" << vox->getY()<< "|" << vox->getZ() << std::endl;
        LOG_ERROR << "Fatal Error. You are trying to insert a voxel at occupied position";
        
        return false;
    }else{
        boost::tuple<int,int,int> index(vox->getX(),vox->getY(),vox->getZ());
        space.insert(std::make_pair( index, vox ));
        return true;
    }
    return false;
}

bool Environment::removeVoxel(Voxel* vox){
    //add cells
    if(voxelExists(vox)){
        boost::tuple<int,int,int> index(vox->getX(),vox->getY(),vox->getZ());
//        LOG_DEBUG << "space size before removing: " << space.size();
        space.erase(index);
//        LOG_DEBUG << "space size after removing: " << space.size();        
        return true;
    }else{
        LOG_ERROR << "You are trying to remove a voxel which is not existent!";
        return false;
    }
}

bool Environment::removeVoxel(boost::tuple<int,int,int> tuple){
    //add cells
    if(tupelExists(tuple)){
        space.erase(tuple);   
    }else{
        LOG_ERROR << "You are trying to remove a voxel which is not existent!";
        return false;
    }
}

bool Environment::voxelExists(Voxel* vox){
        //add cells
    boost::tuple<int,int,int> index(vox->getX(),vox->getY(),vox->getZ());
    if(space.find(index)==space.end()){
        return false;
    }else{
        return true;
    }     
}

bool Environment::tupelExists(boost::tuple<int,int,int> tupl){
        //add cells
    std::map<boost::tuple<int,int,int>,Voxel*>::const_iterator it = space.find(tupl);
    return it!=space.end();
}

bool Environment::tupelExistsWithID(boost::tuple<int,int,int> tupl, int ID){
        //add cells
    std::map<boost::tuple<int,int,int>,Voxel*>::const_iterator it = space.find(tupl);
    Voxel* a = it->second;
    return it!=space.end() && a->getID()==ID;
}

int Environment::getNeighbours(Voxel* vox, int ID){
    return getNeighbours(vox->getX(),vox->getY(),vox->getZ(), ID);
}

int Environment::getNeighbours(int ix, int iy, int iz, int ID){
    int counter = 0;
    int x = ix;
    int y = iy;
    int z = iz;
    boost::tuple<int,int,int> index1(x+1,y,z);
    if(tupelExists(index1)) counter++;
    boost::tuple<int,int,int> index2(x-1,y,z);
    if(tupelExists(index2)) counter++;
    boost::tuple<int,int,int> index3(x,y+1,z);
    if(tupelExists(index3)) counter++;
    boost::tuple<int,int,int> index4(x,y-1,z);
    if(tupelExists(index4)) counter++;
    boost::tuple<int,int,int> index5(x,y,z+1);
    if(tupelExists(index5)) counter++;
    boost::tuple<int,int,int> index6(x,y,z-1);
    if(tupelExists(index6)) counter++;
    return counter;
}

bool Environment::output(std::string path, int time_it){
    //output tissue
    std::string outpath = path+"/cells/tissue.csv";
    io::stream_buffer<io::file_sink> buf(outpath);
    std::ostream                     out(&buf);
    // out writes to cell_$ID_t_$time.txt
    for(int i = 0; i < tissue.size();i++){
        Voxel* vox = tissue.at(i);
        out << vox->getX() <<", " << vox->getY() << ", " << vox->getZ() << ", " << vox->getNB() << std::endl;
    }
    
    
    //output cells
    for(int i = 0; i < cells.size();i++){
        cells.at(i)->output(path, time_it);
    }
}


int Environment::iterate(int iteration){
//    std::cout << "starting an iteration of a cell from " << cells.size() << " cells." << std::endl;
    int succIterations = 0;
    for(int i = 0; i < cells.size();i++){
//        std::cout << "starting an iteration of cell" << i << std::endl;
        int counter = 0;
        Time t1(boost::posix_time::microsec_clock::local_time());
        bool succ = false;
        while(!succ && counter < 100){
            succ = cells.at(i)->iterate();
            counter++;
        }
//        LOG_DEBUG  << "Done iterating a cell with " << counter << " tries.";
        Time t2(boost::posix_time::microsec_clock::local_time());
        TimeDuration dt = t2 - t1;
        long msec = dt.total_milliseconds();
        addRunStatistic("iterationTries",iteration,counter);
        addRunStatistic("iterationTime",iteration,msec/1000.0);
        
        if(counter >=100){
            LOG_FATAL << "Could not iterate cell " << i  << " with front: " << cells.at(i)->getFrontCaniddatesSize() << " and back: " << cells.at(i)->getBackCaniddatesSize() <<  " now changing direction";
            std::cout << "Could not iterate cell " << i  << "changing direction" << std::endl;
            cells.at(i)->setRandomDirection();
            cells.at(i)->initFront();
            return false;
        }else{
            succIterations++;
            std::cout  << "[E] Iteration" << iteration << " | Was iterating cell" << i << " with " << counter <<  " tries in "<< msec/1000.0  << "msec" << std::endl;
//            return true;
        }


    }
//    std::cout << "[E] returning" << succIterations << " successfull iterations" << std::endl;
    return succIterations;
}

bool Environment::addRunStatistic(std::string key, int time, float val){
    runStatType.push_back(key);
    runStatTime.push_back(time);
    runStatVal.push_back(val);
}

void Environment::doRunAnalysis(std::string path){
    std::string outpath = path+"/runstatistics.csv";
    io::stream_buffer<io::file_sink> buf(outpath);
    std::ostream                     out(&buf);
    for(int i = 0; i < runStatType.size();i++){
        out << runStatType.at(i) << ", " << runStatTime.at(i) << ", " << runStatVal.at(i) << std::endl;
    }
    
    //run r script to analyse runstatistics.csv
    std::string cmd = "Rscript runStat.R "+outpath+" >> "+path+"/runanalysis.csv";
    system(cmd.c_str());
    //edit the path+"/runanalysis.csv" because theere are some fragment from the r ouptut
//    std::string cmd2 = "sed -i -- 's/\"//g' "+path+"/runanalysis.csv";
//    system(cmd2.c_str());
//    std::string cmd3 = "sed -i -- 's/\[1\] //g' "+path+"/runanalysis.csv";
//    system(cmd3.c_str());
    //not statistics can be found in the folder for the run, you can easily edit the r script for generation of e.g. plots
    LOG_INFO << "run analysis done";
} 

void Environment::doShapeAnalysis(std::string path){

    //run r script to analyse runstatistics.csv
    std::string cmd = "Rscript shapeAnalysis.R "+path+"/cells >> "+path+"/shapeanalysis.csv";
    system(cmd.c_str());
    //and clean up R output fragments
//    std::string cmd2 = "sed -i -- 's/\"//g' "+path+"/shapeanalysis.csv";
//    system(cmd2.c_str());
//    std::string cmd3 = "sed -i -- 's/\[1\] //g' "+path+"/shapeanalysis.csv";
//    system(cmd3.c_str());
    LOG_INFO << "shape analysis done";

}

float Environment::getRandomFloat(float limit){
    return ((float)(rand()% (int)(limit*1000000))/1000000);
}
