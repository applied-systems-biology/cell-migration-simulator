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

#ifndef CELL_H
#define	CELL_H

#include <cstring>
#include <iostream>
#include <vector>
#include "Voxel.h"
#include <ostream>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/random.hpp>
#include <map>
#include "boost/tuple/tuple_comparison.hpp"
#include <sstream>
#include <boost/foreach.hpp>
#include <math.h>       /* ceil */
#include <set>




#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

class Cell {
public:
    Cell(int);
    Cell(const Cell& orig);
    virtual ~Cell();
    bool init(float, float, float);
    void initFront();
    void updateFront(Voxel*);
    bool output(std::string path, int time_it);
    int getSize(){
        return content.size();
    }
    bool iterate();
    bool setRandomDirection();
    float assignScore(float,int);
    float assignScore(float,int,float);
    
    /**
     * Please note the following.It is necessary to check connectivity of a cell from ime to time. For that i implemented a function
     * checkConnectivity. This function is traversing a cell (via recursive function fillList) and is returning if all voxel from cell could be traversed. then is is connected.
     * However, this function does not work for huge cells (e.g. radius of 6 mym and a resolution of 0.2 ) due to hardware stack limits.
     *  
     * That is why I came up with a checkConnectivityLoop which does the same but not recursive. 
     * It works but is much slower due to a much bigger usage of memory.  
     * 
     * Then i came up with a more efficient algorithm which in best cases only traversing local environment of cell around the removed voxel.
     * In worst case it is still slow but this does barely happen. It seems to work fine and is at least factor 20 faster than the recursive one.  
     */
    bool checkConnectivity();       //rekursive function, has its size limits! but working and moderately fast
    void fillList(int, int, int, int, int);
    
    bool checkConnectivityLoop();   //loop function which should run more stable, but slower
    
    bool checkConnectivityEfficient(boost::tuple<int,int,int>);       //effective version which only looks for specific voxel connections
    bool checkConnection(boost::tuple<int,int,int>, boost::tuple<int,int,int>); //checks for connections between two specific voxels only.
    
    bool removeFromContent(Voxel*);
    void saveFrontPos(Voxel*, float);
    int getID(){
        return this->ID;
    };
    std::vector<Voxel*> getContent(){
        return content;
    }
    float getTrackLength(){
        return track_length;
    }
    int getFrontCaniddatesSize(){
        return fnb.size();
    }
    int getBackCaniddatesSize(){
        return back_positions.size();
    }

    std::string printDirection(){
        std::string res = SSTR(anglex)+","+SSTR(angley)+","+SSTR(anglez);
        return res;
    }
    
    std::string printCenter(){
        std::string res = SSTR(cx)+","+SSTR(cy)+","+SSTR(cz);
        return res;    
    }

//    std::set< boost::tuple<int,int,int> > alibaba; 
    
    
private:
    std::vector<Voxel*> content;
    int ID;
    float cx,cy,cz,cxinit,czinit,cyinit;;                       //center of cell
    float anglex, angley, anglez;       //direction of cell
    float POSWEIGHT, NWEIGHT, fb, DISTWEIGHT;
    std::vector<Voxel*> back_positions;
    std::vector<float> back_scores;
    std::vector<std::vector<Voxel*>::iterator> back_contentIterator;


    float back_sum, front_sum;          //sume of scores of back/front; is needed for sampling;
    boost::random::uniform_real_distribution<> pf01;
    float track_length;
    
    std::map<boost::tuple<int,int,int>,int> fnb;
    //for connectivity test
    std::map<boost::tuple<int,int,int>,bool> connTest;
    
    
    int fillCounter;

    int sampleBack();
    boost::tuple<int,int,int> sampleFront();
    float getDirectionValue(int,int,int);
    float getDirectionValue(boost::tuple<int,int,int>); //return dotproduct between a voxel and the cells direction in [0,1]
    void prepareIteration();        //collect all front and back-positions and clear the storage datastructures before 

    
};

#endif	/* CELL_H */

