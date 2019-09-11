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

#include "Cell.h"
#include "Framework.h"
#include <queue>          // std::queue
#include <iostream>
#include <fstream>
#include "boost/date_time/posix_time/posix_time.hpp" 


#define PI 3.14159265
typedef boost::posix_time::ptime Time;
typedef boost::posix_time::time_duration TimeDuration;


namespace io = boost::iostreams;

Cell::Cell(int ID) {
    this->ID = ID;
}

Cell::Cell(const Cell& orig) {
}

Cell::~Cell() {
}

/*
 to initialize a cell. should be called after creating an instance
 * arguments are 3d position in vocel space
 */
bool Cell::init(float x, float y, float z){
        // i want to init the cells, so i need the radius, and the resolution
//    Framework fw;
    InputManager *im = &Framework::getInputManager();
    Environment *world = &Framework::getEnvironment();
    float size_c = im->getFloatParameter("modelSettings","size_c");
    float resolution = im->getFloatParameter("modelSettings","resolution");

    int vx = round(x/resolution);  //transfor float positions to voxel positions
    int vy = round(y/resolution);
    int vz = round(z/resolution);    
    this->cx =vx;
    this->cy =vy;
    this->cz =vz;
    this->cxinit =x;
    this->cyinit =y;
    this->czinit =z;
    
    int size_vox = ceil(size_c/resolution);
    std::cout << "Initializing a cell with id "<< this->ID << " and size " << size_c << " micrometer and resolution " << resolution << " (micrometer per voxel) which is equal to " << size_vox << " voxels at position " << x << "|" << y << "|" << z << "(" << this->cx << "|" << this->cy << "|" << this->cz << ")" << " now content size is" << content.size() << " and world size is " << world->getSize() << std::endl;
    for(int xr = -size_vox; xr < size_vox; xr++){
        for(int yr = -size_vox; yr < size_vox; yr++){
            for(int zr = -size_vox; zr < size_vox; zr++){
//                std::cout << "block distance" << sqrt(xr*xr+yr*yr+zr*zr) << std::endl;
                if(sqrt(xr*xr+yr*yr+zr*zr)*resolution < size_c){
//                    std::cout << "adding voxel" << i++ << std::endl;
                    //create the voxel object 
                    Voxel* vox = new Voxel(xr+vx,yr+vy,zr+vz, this->ID);
                    if(world->addVoxel(vox)){
                        content.push_back(vox);
//                        std::cout << "size of world:" << content.size() << std::endl;
                    }else{
                        //it failed to insert a voxel. so i have to remove the inserted ones and throw a false as return
                        for(int i = 0; i < content.size();i++){
                            world->removeVoxel(content.at(i));
                        }
                        content.clear();
//                        std::cout << "Init of cell did fail at voxel " << vox->getX() << "|"<< vox->getY() << "|"<< vox->getZ() << "now content size is" << content.size() << " and world size is " << world->getSize() << std::endl;
                        LOG_INFO << "Init of cell did fail at voxel " << vox->getX() << "|"<< vox->getY() << "|"<< vox->getZ() << std::endl;
                        return false;
                    }
                }
            }
        }
    }
//    std::cout << "after init the cell content has a size of " << content.size() << std::endl;
           
    //init a random direction  
    this->setRandomDirection();

    track_length = 0;
    
    //load parameters fro movement
    POSWEIGHT = im->getFloatParameter("modelSettings","pos_weight");
    NWEIGHT = im->getFloatParameter("modelSettings","neighbour_weight");
    DISTWEIGHT = im->getFloatParameter("modelSettings","distance_weight");
    fb = im->getFloatParameter("modelSettings","frontback");
    
    return true;
}

/*
 the algorithm needs to compute possible candidate positions at the cell front, these candidates are not computed new in each iteration but are dynamically updates.
 * this function initializes the information on front candidates
 */
void Cell::initFront(){
    Environment *world = &Framework::getEnvironment();
    //go through all voxel of the cell, check if it is at the border (nb is initialized) and at front
    //then put it on the     std::map<boost::tuple<int,int,int>,int> fnb; map when done, create an update function
    for(int i = 0; i < content.size(); i++){
        int nb = content.at(i)->getNB();
        if(nb<6){
            //get the direciton score 
            float posvalue = getDirectionValue(content.at(i)->getX(),content.at(i)->getY(),content.at(i)->getZ());
            if(posvalue >= fb){
                //we now have all voxel which are at the front rear
                //i now need to get all the neighbour positions of this, check if it is free and then add it to std::map<boost::tuple<int,int,int>,int> fnb
                int x = content.at(i)->getX();
                int y = content.at(i)->getY();
                int z = content.at(i)->getZ();
                
                int free = world->getNeighbours(x+1,y,z, this->ID);
                boost::tuple<int,int,int> ind1(x+1,y,z);
//                std::cout << posvalue << " vs. " << getDirectionValue(x+1,y,z)  << std::endl;
                if( ! world->tupelExists(ind1)){ fnb.insert(std::make_pair(ind1,free));}
                
                free = world->getNeighbours(x-1,y,z, this->ID);
                boost::tuple<int,int,int> ind2(x-1,y,z);
                if( ! world->tupelExists(ind2)){ fnb.insert(std::make_pair(ind2,free));}
                
                free = world->getNeighbours(x,y+1,z, this->ID);
                boost::tuple<int,int,int> ind3(x,y+1,z);
                if( ! world->tupelExists(ind3)){ fnb.insert(std::make_pair(ind3,free));}
                
                free = world->getNeighbours(x,y-1,z, this->ID);
                boost::tuple<int,int,int> ind4(x,y-1,z);
                if( ! world->tupelExists(ind4)){ fnb.insert(std::make_pair(ind4,free));}
                
                free = world->getNeighbours(x,y,z+1, this->ID);
                boost::tuple<int,int,int> ind5(x,y,z+1);
                if( ! world->tupelExists(ind5)){ fnb.insert(std::make_pair(ind5,free));}
                
                free = world->getNeighbours(x,y,z-1, this->ID);
                boost::tuple<int,int,int> ind6(x,y,z-1);
                if( ! world->tupelExists(ind6)){ fnb.insert(std::make_pair(ind6,free));}
            }
        }
    }
    int counter =0;
    int counter2 = 0;
    std::map<boost::tuple<int,int,int>,int>::iterator it;
    for(it= fnb.begin(); it != fnb.end();it++){ 
        //calculate new position score of index, if below fb, throw it our, else increment iterator
            float vecx = (float)boost::tuples::get<0>((*it).first) - cx;      //relational vector of xovel
            float vecy = (float)boost::tuples::get<1>((*it).first) - cy;
            float vecz = (float)boost::tuples::get<2>((*it).first) - cz;
//            std::cout << vecx << ", " << vecy << ", " << vecz << std::endl;
            //normalize to 1
            float fvecx = ((float)vecx)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float fvecy = ((float)vecy)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float fvecz = ((float)vecz)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float posvalue = fvecx*anglex + fvecy*angley + fvecz*anglez;
//            std::cout << posvalue << std::endl;
            if(posvalue!=posvalue){	//true if posvalue = NA
                    posvalue = 0;// set NA to zero
            }
            if(posvalue < 0){
                counter++;
            }else{
                counter2++;
            }
    }
//    std::cout << "after front-nb initialization there are " << fnb.size() << " positions " << counter << ", " << counter2 << std::endl;
}

/*
 the algorithm needs to compute possible candidate positions at the cell front, these candidates are not computed new in each iteration but are dynamically updates.
 * this function does the update
 */
void Cell::updateFront(Voxel* vox){
    Environment *world = &Framework::getEnvironment();
    //this function is trying to update the front-neighbour counts
    //for that i have to remove the one which was occupied, given as argument, and update the information for all of its surrounding
    int x = vox->getX();
    int y = vox->getY();
    int z = vox->getZ();
    //remove
    boost::tuple<int,int,int> index(x,y,z);
    fnb.erase(index);
    //update neighbour information

    int free = world->getNeighbours(x+1,y,z, this->ID);
    boost::tuple<int,int,int> ind1(x+1,y,z);
    if( ! world->tupelExists(ind1)){ fnb.insert(std::make_pair(ind1,free));}

    free = world->getNeighbours(x-1,y,z, this->ID);
    boost::tuple<int,int,int> ind2(x-1,y,z);
    if( ! world->tupelExists(ind2)){ fnb.insert(std::make_pair(ind2,free));}

    free = world->getNeighbours(x,y+1,z, this->ID);
    boost::tuple<int,int,int> ind3(x,y+1,z);
    if( ! world->tupelExists(ind3)){ fnb.insert(std::make_pair(ind3,free));}

    free = world->getNeighbours(x,y-1,z, this->ID);
    boost::tuple<int,int,int> ind4(x,y-1,z);
    if( ! world->tupelExists(ind4)){ fnb.insert(std::make_pair(ind4,free));}

    free = world->getNeighbours(x,y,z+1, this->ID);
    boost::tuple<int,int,int> ind5(x,y,z+1);
    if( ! world->tupelExists(ind5)){ fnb.insert(std::make_pair(ind5,free));}

    free = world->getNeighbours(x,y,z-1, this->ID);
    boost::tuple<int,int,int> ind6(x,y,z-1);
    if( ! world->tupelExists(ind6)){ fnb.insert(std::make_pair(ind6,free));}
//    std::cout << "center: " << cx << ", " << cy << ", " << cz << std::endl;

    //then i have to throw out all whichs positionvalue has changed due to a new center of the cell,
    std::map<boost::tuple<int,int,int>,int>::iterator it;
    std::vector<boost::tuple<int,int,int> > toremove; 
    for(it= fnb.begin(); it != fnb.end();it++){
        if(!world->tupelExists((*it).first)){
                //calculate new position score of index, if below fb, throw it our, else increment iterator
            float vecx = (float)boost::tuples::get<0>((*it).first) - cx;      //relational vector of xovel
            float vecy = (float)boost::tuples::get<1>((*it).first) - cy;
            float vecz = (float)boost::tuples::get<2>((*it).first) - cz;
//            std::cout << vecx << ", " << vecy << ", " << vecz << std::endl;
            //normalize to 1
            float fvecx = ((float)vecx)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float fvecy = ((float)vecy)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float fvecz = ((float)vecz)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float posvalue = fvecx*anglex + fvecy*angley + fvecz*anglez;
//            std::cout << posvalue << std::endl;
            if(posvalue!=posvalue){	//true if posvalue = NA
                    posvalue = 0;// set NA to zero
            }
            if(posvalue < fb){
//                fnb.erase(it);
                toremove.push_back((*it).first);
            }
        }else{
            toremove.push_back((*it).first);        
        }

    }
    //and finally delete all of them
    for(int i = 0; i < toremove.size();i++){
        fnb.erase(toremove.at(i));
    }
//    std::cout << "after updating the front-candidates-nb-count i have" << fnb.size() << " possible candidates" << " and i removed:" << toremove.size() << std::endl;
    toremove.clear();
}

/*
 give cell new random direction
 */
bool Cell::setRandomDirection(){
    //give the cell a random direction
    float theta = (rand() % ((int)(PI * 2000000)))/1000000.0f; //random angle in [0,2PI]#
    float z = ((rand() % 20000000)/10000000.0f)-1.0f;
    this->anglex = sqrt(1-(z*z))*cos(theta);
    this->angley = sqrt(1-(z*z))*sin(theta);
    this->anglez = z;  
}

/*
 write the cell to a given path
 * output format is x,y,z (all in voxel space), number of nieghbors
 */
bool Cell::output(std::string path, int time_it){
        Environment *world = &Framework::getEnvironment();
//    std::cout << "I now do ouput the cell-content to " << path <<" at timepoint " << time_it << std::endl;
    std::string outpath = path+"/cells/cell"+SSTR(ID)+"_t"+SSTR(time_it)+".csv";
    io::stream_buffer<io::file_sink> buf(outpath);
    std::ostream                     out(&buf);
//    std::cout << content.size() << std::endl;
    // out writes to cell_$ID_t_$time.txt
    for(int i = 0; i < this->content.size();i++){
        Voxel* vox = content.at(i);
        out << vox->getX() <<", " << vox->getY() << ", " << vox->getZ() << ", " << vox->getNB() << std::endl;
//        std::cout << vox->getX() <<", " << vox->getY() << ", " << vox->getZ() << std::endl;
    }
    
    
    //and also output the last cell for testing purpose
    std::string outpath2 = "last.csv";
    io::stream_buffer<io::file_sink> buf1(outpath2);
    std::ostream                     out2(&buf1);
    for(int i = 0; i < this->content.size();i++){
        Voxel* vox = content.at(i);
        out2 << vox->getX() <<", " << vox->getY() << ", " << vox->getZ() << std::endl;
//        std::cout << vox->getX() <<", " << vox->getY() << ", " << vox->getZ() << std::endl;
    }
}

/*
 do one iteration of the cell -> remove a SU (voxel) from the rear and put it on front
 */
bool Cell::iterate(){
    Environment *world = &Framework::getEnvironment();  
    prepareIteration(); //getting the front and back-positions and clear the collected data from last iteration
//    LOG_DEBUG  << "Done prparing an iteration";

    
    ////remove one pixel
    // do monte carlo sampling

    
//    std::cout << " back position was found: " << pos  <<" :" << a->getX() << "|" << a->getY() << "|" << a->getZ() << std::endl;
    //remove
    Voxel* a;   //the one to be removed later
    int pos = sampleBack();     //sample a position
    std::cout << "position for back sampling: " << pos << std::endl;
    a = back_positions.at(pos);
    bool remove = world->removeVoxel(a);      //remove from world.
//    LOG_DEBUG << "trying to remove " << a->print() << " with success: "  << remove;
    content.erase(back_contentIterator.at(pos));    //remove from cell content
    //update neighbourhood around removed voxel
    int x = a->getX();
    int y = a->getY();
    int z = a->getZ();
//    std::cout << x << ", " << y << ", " << z << std::endl;
    for(int a = x-1; a <=x+1;a++){
        for(int b = y-1; b <= y+1;b++){
            for(int c = z-1; c <= z+1;c++){
                //if there is a voxel, update its neighbour-number
                boost::tuple<int,int,int> nbindex(a,b,c);
//                std::cout << a << ", " << b << ", " << c << std::endl;
                if(world->tupelExists(nbindex)){
                    world->getVoxel(nbindex)->setNB(world->getNeighbours(a,b,c, this->ID));
                }
            }
        }
    }
    
    if(remove){
//        LOG_DEBUG  << "starting front part" << std::endl;
        //check for connectivity
        boost::tuple<int,int,int> removed(a->getX(),a->getY(),a->getZ());
//        bool conn = checkConnectivityEfficient(removed);
        bool conn = checkConnectivity();
//        bool conn = true;

        if(conn){
            boost::tuple<int,int,int> front_index;  
            front_index = sampleFront();
            bool occupied = world->tupelExists(front_index);
            int counter = 0;
            while(occupied || boost::tuples::get<0>(front_index)==INT_MAX){         //if the front-sampling has no candidates, it return a index with all value MAX_INT. sampleFront then automatically sets a new direction and initialises the front-candidates stack by itself. a new repeat of sampleFront in the next iteration should bring a different result.
                front_index = sampleFront();
                occupied = world->tupelExists(front_index);
                counter++;
                if(counter > 1000){
                    //the sampled front idex is already occupied. so i need to re-insert the removed voxel and return false   
                    LOG_WARNING << "A free front position could not be found! Turning the direction";
                    std::cout  << "A free front position could not be found! Turning the direction." << std::endl;
                    //so i need to reinsert the removed voxel and return a false because i have to re-try it.
                    world->addVoxel(a);
                    content.push_back(a);
                    setRandomDirection();
                    return false;  
                }
            }
            //insert one pixel           
//            std::cout << "done sampling front !" << std::endl;
            Voxel* vox = new Voxel(boost::tuples::get<0>(front_index),boost::tuples::get<1>(front_index),boost::tuples::get<2>(front_index), this->ID);
            bool add = world->addVoxel(vox);
            content.push_back(vox);
            //update neighbourhoodinformation of the voxels
            x = vox->getX();
            y = vox->getY();
            z = vox->getZ();
            for(int a = x-1; a <=x+1;a++){
                for(int b = y-1; b <= y+1;b++){
                    for(int c = z-1; c <= z+1;c++){
                        //if there is a voxel, update its neighbour-number
                        boost::tuple<int,int,int> nbindex(a,b,c);
                        if(world->tupelExists(nbindex)){
                            world->getVoxel(nbindex)->setNB(world->getNeighbours(a,b,c,this->ID));
                        }
                    }
                }
            }



            //maybe do some other things
            //the cell has moved slightly so i should update its center 
            float sumx = 0;
            float sumy = 0;
            float sumz = 0;
            for(int i = 0; i <content.size();i++){
                sumx += content.at(i)->getX();
                sumy += content.at(i)->getY();
                sumz += content.at(i)->getZ();

            }
            sumx /= content.size();
            sumy /= content.size();
            sumz /= content.size();
            //calculate dist to new center
            track_length= sqrt((cxinit-sumx)*(cxinit-sumx)+(cyinit-sumy)*(cyinit-sumy)+(czinit-sumz)*(czinit-sumz));
            cx = sumx;
            cy = sumy;
            cz = sumz;
            updateFront(vox);
            
//            std::cout << "new center: " << sumx << ", "<< sumx << ", " << sumz << std::endl;

            //reply
            if(add){
//            LOG_DEBUG << "Cwell iteration succeeded";                
                return true;
            }             
        }else{
            LOG_WARNING << "The cell has lost connectivity!";
            std::cout  << "The cell has lost connectivity! by moving " << a->getX() << a->getY() << a->getZ() << std::endl;;
            //so i need to reinsert the removed voxel and return a false because i have to re-try it.
            world->addVoxel(a);
            content.push_back(a);
            
            return false;
        }
    }
    std::cout << "could not remove a voxel" << std::endl;
    
    return false;
}

/* pos: position of voxel relative to migration direction, 1=front-pixel of cell, -1 = last vxl of cell
 free: number of free positions besided the voxel
 */
float Cell::assignScore(float pos, int free){
    if(pos<0)pos*=-1;
//    std::cout << "pos:" << pos <<"free:" << free<< std::endl;
    return pow(pos,POSWEIGHT) * pow(free,NWEIGHT);    
}

//this function in addition to the previous one is also accountign for the distance of a voxel so that ones with higher distance get higher scores.
float Cell::assignScore(float pos, int free,float dist){
    if(pos<0)pos*=-1;
//    std::cout << "pos:" << pos <<"free:" << free<< std::endl;
    return pow(pos,POSWEIGHT) * pow(free,NWEIGHT)*pow(dist,DISTWEIGHT);    
}

/*
 cell connectivity test, recursive version, fails for big cells due to hardware limits in the processor stack
 */
bool Cell::checkConnectivity(){
    Time t1(boost::posix_time::microsec_clock::local_time());
    //start with a random voxel of the cell
    Voxel* a = *content.begin();
    connTest.clear();
    fillList(a->getX(),a->getY(),a->getZ(),0,0);
//    std::cout << "found voxels of cell"<< connTest.size() << " and we compare it to content:" << content.size() << std::endl;
    Time t2(boost::posix_time::microsec_clock::local_time());
    TimeDuration dt = t2 - t1;
    long msec = dt.total_milliseconds();
//    std::cout << msec/1000.0 << "s needed for connectivity test" << std::endl;
    if(connTest.size()==content.size()){
        return true;
    }else {
        return false;
    }
}


/*
 cell connectivity test, non-recursive version
 */
bool Cell::checkConnectivityLoop(){
    Time t1(boost::posix_time::microsec_clock::local_time());
    std::map<boost::tuple<int,int,int>,bool> visited;       //this stores the voxel i already have visited!
    std::queue<boost::tuple<int,int,int>> candidates;       //candidates list for traversion
    Environment *world = &Framework::getEnvironment();     //environment
    Voxel* first = content.at(1);
    boost::tuple<int,int,int> index(first->getX(),first->getY(),first->getZ());               
    candidates.push(index);
    int counter = 0;
    do{ 
        //get the first of candidates,
        boost::tuple<int,int,int> current = candidates.front();
        int x = current.get<0>();
        int y = current.get<1>();
        int z = current.get<2>();
        //if not already visited
        std::map<boost::tuple<int,int,int>,bool>::const_iterator it = visited.find(current);    //create iterator over visited
        bool currentExists =  it!=visited.end();    //check if find did return the end, if so, the current tuple is not visited yet
        if(!currentExists){
            //add to visited---------
            visited.insert(std::make_pair( current, true ));
            //check neighbours of first element and put them on the candidates list when they do exist in world

            //top
            boost::tuple<int,int,int> top(x+1,y,z);
            if(world->tupelExistsWithID(top,this->ID)) candidates.push(top);
            boost::tuple<int,int,int> bottom(x-1,y,z);
            if(world->tupelExistsWithID(bottom,this->ID)) candidates.push(bottom);
            boost::tuple<int,int,int> left(x,y+1,z);
            if(world->tupelExistsWithID(left,this->ID)) candidates.push(left);
            boost::tuple<int,int,int> right(x,y-1,z);
            if(world->tupelExistsWithID(right,this->ID)) candidates.push(right);
            boost::tuple<int,int,int> front(x,y,z+1);
            if(world->tupelExistsWithID(front,this->ID)) candidates.push(front);
            boost::tuple<int,int,int> back(x,y,z-1);
            if(world->tupelExistsWithID(back,this->ID)) candidates.push(back);
//            std::cout << "doing this. " << std::endl;        
            counter = counter + 6;
        }
        candidates.pop();
//        std::cout << "done checking index " << x << ", " << y << ", " << z << ". Now i do have" << candidates.size() << " candidates ||" << counter << std::endl;  
        //pop first from list
    } while (candidates.size() > 0);
//    std::cout << "cell was traversed over" << 0 << " voxels" << std::endl;
    Time t2(boost::posix_time::microsec_clock::local_time());
    TimeDuration dt = t2 - t1;
    long msec = dt.total_milliseconds();
    std::cout << msec/1000.0 << "s needed for connectivity test" << std::endl;
    if( visited.size()==content.size()){
        return true;
    }else{
        return false;
    }
}

/*
 cell connectivity test, efficient version
 * needs to be testes more!
 */
bool Cell::checkConnectivityEfficient(boost::tuple<int,int,int> removed){
    //start with a random voxel of the cell
//    std::cout << "starting efficient conCheck after removal of voxel" << removed.get<0>() << "|"<< removed.get<1>() << "|"<< removed.get<2>() << std::endl;
    boost::tuple<int,int,int> top   (removed.get<0>()+1, removed.get<0>(), removed.get<0>());
    boost::tuple<int,int,int> bottom(removed.get<0>()-1, removed.get<0>(), removed.get<0>());
    boost::tuple<int,int,int> left  (removed.get<0>(), removed.get<0>()+1, removed.get<0>());
    boost::tuple<int,int,int> right (removed.get<0>(), removed.get<0>()-1, removed.get<0>());
    boost::tuple<int,int,int> front (removed.get<0>(), removed.get<0>(), removed.get<0>()+1);
    boost::tuple<int,int,int> back  (removed.get<0>(), removed.get<0>(), removed.get<0>()-1);
    Environment *world = &Framework::getEnvironment();     //environment
    //top to ...
    if(world->tupelExistsWithID(top, this->ID) && world->tupelExistsWithID(bottom, this->ID)){
        bool foundConnection = checkConnection(top, bottom);
        if(!foundConnection) return false;
    }
    if(world->tupelExistsWithID(top, this->ID) && world->tupelExistsWithID(left, this->ID)){
        bool foundConnection = checkConnection(top, left);
        if(!foundConnection) return false;
    }
    if(world->tupelExistsWithID(top, this->ID) && world->tupelExistsWithID(right, this->ID)){
        bool foundConnection = checkConnection(top, right);
        if(!foundConnection) return false;
    }
    if(world->tupelExistsWithID(top, this->ID) && world->tupelExistsWithID(front, this->ID)){
        bool foundConnection = checkConnection(top, front);
        if(!foundConnection) return false;
    }
    if(world->tupelExistsWithID(top, this->ID) && world->tupelExistsWithID(bottom, this->ID)){
        bool foundConnection = checkConnection(top, back);
        if(!foundConnection) return false;
    }
    //bottom to ...
    if(world->tupelExistsWithID(bottom, this->ID) && world->tupelExistsWithID(left, this->ID)){
        bool foundConnection = checkConnection(bottom, left);
        if(!foundConnection) return false;
    }
    if(world->tupelExistsWithID(bottom, this->ID) && world->tupelExistsWithID(right, this->ID)){
        bool foundConnection = checkConnection(bottom, right);
        if(!foundConnection) return false;
    }
    if(world->tupelExistsWithID(bottom, this->ID) && world->tupelExistsWithID(front, this->ID)){
        bool foundConnection = checkConnection(bottom, front);
        if(!foundConnection) return false;
    }
    if(world->tupelExistsWithID(bottom, this->ID) && world->tupelExistsWithID(back, this->ID)){
        bool foundConnection = checkConnection(bottom, back);
        if(!foundConnection) return false;
    }
    //left to ...
    if(world->tupelExistsWithID(left, this->ID) && world->tupelExistsWithID(right, this->ID)){
        bool foundConnection = checkConnection(left, right);
        if(!foundConnection) return false;
    }   
    if(world->tupelExistsWithID(left, this->ID) && world->tupelExistsWithID(front, this->ID)){
        bool foundConnection = checkConnection(left, front);
        if(!foundConnection) return false;
    }   
    if(world->tupelExistsWithID(left, this->ID) && world->tupelExistsWithID(back, this->ID)){
        bool foundConnection = checkConnection(left, back);
        if(!foundConnection) return false;
    }
    //right to ...
    if(world->tupelExistsWithID(right, this->ID) && world->tupelExistsWithID(front, this->ID)){
        bool foundConnection = checkConnection(right, front);
        if(!foundConnection) return false;
    }   
    if(world->tupelExistsWithID(right, this->ID) && world->tupelExistsWithID(back, this->ID)){
        bool foundConnection = checkConnection(right, back);
        if(!foundConnection) return false;
    }
    //front to ...
    if(world->tupelExistsWithID(front, this->ID) && world->tupelExistsWithID(back, this->ID)){
        bool foundConnection = checkConnection(front, back);
        if(!foundConnection) return false;
    }
    return true;
}
/*
 assisting function for checkConnectivityEfficient
 */
bool Cell::checkConnection(boost::tuple<int,int,int> start, boost::tuple<int,int,int> aim){
    std::string outpath = "/home/blickensdorf/connection.csv";
    std::ofstream outfile;
    outfile.open(outpath, std::ios_base::trunc);

    std::map<boost::tuple<int,int,int>,bool> visited;       //this stores the voxel i already have visited!
    std::queue<boost::tuple<int,int,int>> candidates;       //candidates list for traversion
    Environment *world = &Framework::getEnvironment();     //environment
    
    
    candidates.push(start);
    int counter = 1;
    do{ 
        //get the first of candidates,
        boost::tuple<int,int,int> current = candidates.front();
        int x = current.get<0>();
        int y = current.get<1>();
        int z = current.get<2>();
        //compare to aim 
        int ax = aim.get<0>();
        int ay = aim.get<1>();
        int az = aim.get<2>();        
        if(x==ax && y==ay &&  z==az){
            outfile  << counter << std::endl;
            return true;
        }
        
        //if not already visited
        std::map<boost::tuple<int,int,int>,bool>::const_iterator it = visited.find(current);    //create iterator over visited
        bool currentExists =  it!=visited.end();    //check if find did return the end, if so, the current tuple is not visited yet
        if(!currentExists){
            //add to visited---------
            visited.insert(std::make_pair( current, true ));
            //check neighbours of first element and put them on the candidates list when they do exist in world

            //top
            boost::tuple<int,int,int> top(x+1,y,z);
            if(world->tupelExistsWithID(top, this->ID)) candidates.push(top);
            boost::tuple<int,int,int> bottom(x-1,y,z);
            if(world->tupelExistsWithID(bottom, this->ID)) candidates.push(bottom);
            boost::tuple<int,int,int> left(x,y+1,z);
            if(world->tupelExistsWithID(left, this->ID)) candidates.push(left);
            boost::tuple<int,int,int> right(x,y-1,z);
            if(world->tupelExistsWithID(right, this->ID)) candidates.push(right);
            boost::tuple<int,int,int> front(x,y,z+1);
            if(world->tupelExistsWithID(front, this->ID)) candidates.push(front);
            boost::tuple<int,int,int> back(x,y,z-1);
            if(world->tupelExistsWithID(back, this->ID)) candidates.push(back);
//            std::cout << "doing this. " << std::endl;        
            counter = counter + 1;
        }
        candidates.pop();
//        std::cout << "done checking index " << x << ", " << y << ", " << z << ". Now i do have" << candidates.size() << " candidates ||" << counter << std::endl;  
        //pop first from list
    } while (candidates.size() > 0);
    outfile  << counter << std::endl;
    return false;
}

/*
 assisting function for recursive connectivity check
 */
void Cell::fillList(int x, int y, int z, int c, int l){
//    std::cout << "checking voxel" << a->print() << ", " << connTest.size() << std::endl;
    Environment *world = &Framework::getEnvironment();
    boost::tuple<int,int,int> index(x,y,z);
    //insert a into
    //check if the call voxel a is in the environment
    bool exists = world->tupelExistsWithID(index, this->ID);
    if(exists){
        fillCounter++;
//        std::cout << fillCounter << ", " << c << ", " << l <<  std::endl;
        //check if the voxel is already visited on the cell traversal
        if(!connTest.count(index)){
            //it is unknown so I should add it
            connTest.insert(std::map<boost::tuple<int,int,int>,bool>::value_type(index,true));
            //it is unknown, so I should visit its neighbours
            fillList(x+1,y,z,c+1,1);
            fillList(x-1,y,z,c+1,2);
            fillList(x,y+1,z,c+1,3);
            fillList(x,y-1,z,c+1,4);
            fillList(x,y,z+1,c+1,5);
            fillList(x,y,z-+1,c+1,6);
        }
    }
}

/*
 outdated function
 */
void Cell::saveFrontPos(Voxel* vox,float score){
//    boost::tuple<int,int,int> index(vox->getX(),vox->getY(),vox->getZ());
//    front_sum+=score;
//    frontCandidates.insert(std::map<boost::tuple<int,int,int>,float>::value_type(index,score));       //using a map makes removal of dublicates efficient, but i have to put all to a vector later for sampling
}


/*
 monte carle acceptance rejection sampling for the cell rear candidates.
 */
int Cell::sampleBack(){
    std::cout << "starting backsamling with candidates:" << back_scores.size() << std::endl;
    float p = pf01(Framework::rng);
    int pos = -1;
    double sum = 0;
    for(int i = 0; i < back_scores.size();i++){
        sum += back_scores.at(i)/back_sum;
        if(sum>p){
            pos = i;
            break;
        }
    }
    return pos;
}

/*
 monte carle acceptance rejection sampling for the cell front candidates.
 */
boost::tuple<int,int,int> Cell::sampleFront(){
    //first remove occupied positions from this list.
    Environment *world = &Framework::getEnvironment();
    std::map<boost::tuple<int,int,int>,int>::iterator it;
    std::vector<boost::tuple<int,int,int> > toremove; 
    for(it= fnb.begin(); it != fnb.end();it++){
        if(world->tupelExists((*it).first)){
            toremove.push_back((*it).first);        
        }
    }
    for(int i = 0; i < toremove.size();i++){
        fnb.erase(toremove.at(i));
    }
    

    float sum_scores = 0;
    std::map<boost::tuple<int,int,int>,int>::iterator ito;
    for(ito= fnb.begin(); ito != fnb.end();ito++){
        float score = assignScore(getDirectionValue((*ito).first), (*ito).second);
        sum_scores +=score; 
    }
//    std::cout << "have a fnb size of : " << fnb.size() << std::endl;
    if(fnb.size()==0){
        //the cell is stuck and can not move into this direction! There are no free positions to move.
        setRandomDirection();
        initFront();
        boost::tuple<int,int,int> front_index(INT_MAX,INT_MAX,INT_MAX);
        return(front_index);
    }
    boost::tuple<int,int,int> front_index;
    double sum = 0;
    float p = pf01(Framework::rng);
//            std::map<boost::tuple<int,int,int>,int>::iterator ito;
    for(ito= fnb.begin(); ito != fnb.end();ito++){
        float score = assignScore(getDirectionValue((*ito).first), (*ito).second)/sum_scores; 
        sum += score;
        if(sum > p){
            front_index = (*ito).first;
            break;
        }            

    }
    int x = boost::tuples::get<0>(front_index);
    return front_index;
}


/*
 necessary computations before a cell 'move'
 */
void Cell::prepareIteration(){
    std::cout << "preparing the iterations " << content.size() << std::endl;
    Environment *world = &Framework::getEnvironment();  
    
    //find back and front pixel and compute measures for voxel picking
    back_sum = 0;
    front_sum = 0;
    back_positions.clear();
    back_scores.clear();
    back_contentIterator.clear();
    
//    std::vector<Voxel*>::iterator it = content.begin();
    float pretime = 0;
    float backtime = 0;
    float fronttime = 0;
    
    
    float tcreate = 0;
    float tnei = 0;
    float texist = 0;
    float tsave = 0;
    
    
    int counter = 0;
        int plus = 0;
            int minus = 0;
    bool wrongNBcount = false;
    float sumx = 0;
    float sumy = 0;
    float sumz = 0;

    InputManager *im = &Framework::getInputManager();
    float resolution = im->getFloatParameter("modelSettings","resolution");

    float maxBackDist = 0;

    for(std::vector<Voxel*>::iterator it = content.begin(); it != content.end(); ++it) {
        int vecx = (*it)->getX() - cx;      //relational vector of xovel
        int vecy = (*it)->getY() - cy;
        int vecz = (*it)->getZ() - cz;
        sumx+= (*it)->getX();
        sumy+= (*it)->getY();
        sumz+= (*it)->getZ();                
        //normalize to 1
        float fvecx = ((float)vecx)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
        float fvecy = ((float)vecy)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
        float fvecz = ((float)vecz)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
        float posvalue = fvecx*anglex + fvecy*angley + fvecz*anglez;
        if(posvalue!=posvalue){	//true if posvalue = NA
                posvalue = 0;// set NA to zero

        }
        int distx = vecx*resolution;      //relational vector of xovel to center
        int disty = vecy*resolution;
        int distz = vecz*resolution;
        float dist = sqrt(pow(distx,2)+pow(disty,2)+pow(distz,2));
        if(dist > maxBackDist){maxBackDist = dist;}
    }
    
    for(std::vector<Voxel*>::iterator it = content.begin(); it != content.end(); ++it) {
        counter++;
//    BOOST_FOREACH(Voxel* now_voxel, content ){	//iterate over all voxel of content
	 	int vecx = (*it)->getX() - cx;      //relational vector of xovel
	 	int vecy = (*it)->getY() - cy;
	 	int vecz = (*it)->getZ() - cz;
                sumx+= (*it)->getX();
                sumy+= (*it)->getY();
                sumz+= (*it)->getZ();                
                //normalize to 1
                float fvecx = ((float)vecx)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
	 	float fvecy = ((float)vecy)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
	 	float fvecz = ((float)vecz)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
	 	float posvalue = fvecx*anglex + fvecy*angley + fvecz*anglez;
		if(posvalue!=posvalue){	//true if posvalue = NA
			posvalue = 0;// set NA to zero
                        
		}
                
                if(posvalue < fb){
                    minus++;
                    //voxel belongs to back
                    //save the voxel
                    //get number of neighbout
                    boost::tuple<int,int,int> index((*it)->getX(),(*it)->getY(),(*it)->getZ());                    
                    int nbs = world->getVoxel(index)->getNB();
                    int free = 6 - nbs;
//                    std::cout << "free:" << free << std::endl;
                    if(free>0){

                        int distx = vecx*resolution;      //relational vector of xovel to center
                        int disty = vecy*resolution;
                        int distz = vecz*resolution;
                        float dist = sqrt(pow(distx,2)+pow(disty,2)+pow(distz,2));
                        
                        
                        float score = assignScore(posvalue,free,dist/(maxBackDist));
                        back_positions.push_back((*it));
                        back_scores.push_back(score);
                        back_contentIterator.push_back(it);
                        back_sum+=score;
                    }
                }else{
                    plus++;
                }
                
                if(counter%100==0){
//                    std::cout << "asfdsad" << plus << "|" << minus << "posvalue:" << posvalue << " and direction" << anglex << "|" << anglex << "|" << anglex << "at center" << cx << "|" << cy << "|" << cz << std::endl;
                }

    }
    std::cout << "after preparing the iterations we have " << back_positions.size() << "back-candidates and calculated center of" << sumx/content.size() <<"|"<< sumy/content.size() <<"|"<< sumz/content.size() << "vs real center of " << cx << "|" << cy << "|" << cz  << std::endl;

}

/*
 get direction vector of a cell SU/voxel
 */
float Cell::getDirectionValue(int x, int y, int z){
            ;
            int vecx = x - cx;      //relational vector of xovel
            int vecy = y - cy;
            int vecz = z - cz;
            //normalize to 1
            float fvecx = ((float)vecx)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float fvecy = ((float)vecy)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float fvecz = ((float)vecz)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float posvalue = fvecx*anglex + fvecy*angley + fvecz*anglez;

            if(posvalue!=posvalue){	//true if posvalue = NA
                    posvalue = 0;// set NA to zero
            }
            return posvalue;
}
/*
 get direction vector of a cell SU/voxel
 */
float Cell::getDirectionValue(boost::tuple<int,int,int> index){
            int x = boost::tuples::get<0>(index);
            int y = boost::tuples::get<1>(index);
            int z = boost::tuples::get<2>(index);
            int vecx = x - cx;      //relational vector of xovel
            int vecy = y - cy;
            int vecz = z - cz;
            //normalize to 1
            float fvecx = ((float)vecx)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float fvecy = ((float)vecy)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float fvecz = ((float)vecz)/sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
            float posvalue = fvecx*anglex + fvecy*angley + fvecz*anglez;

            if(posvalue!=posvalue){	//true if posvalue = NA
                    posvalue = 0;// set NA to zero
            }
            return posvalue;
}
