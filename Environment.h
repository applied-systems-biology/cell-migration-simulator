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
#include "Voxel.h"
#include <map>
#include "Cell.h"
#include <vector>
#include "boost/tuple/tuple_comparison.hpp"
#include <boost/random.hpp>
#include "string"
#include <cstdlib>      // std::rand, std::srand

#ifndef ENVIRONMENT_H
#define	ENVIRONMENT_H

class Environment {
public:
    Environment();
    Environment(const Environment& orig);
    virtual ~Environment();
    void init();
    bool addVoxel(Voxel*);
    bool removeVoxel(Voxel*);
    bool removeVoxel(boost::tuple<int,int,int>);
    bool voxelExists(Voxel*);
    bool tupelExists(boost::tuple<int,int,int>);
    bool tupelExistsWithID(boost::tuple<int,int,int>, int);
    bool output(std::string path, int time_it);
    int iterate(int iteration);
    int getSize(){
        return space.size();
    }
    bool testIntegrity();
    bool addRunStatistic(std::string, int, float);
    void doRunAnalysis(std::string);
    void doShapeAnalysis(std::string);    
    int getNeighbours(Voxel*, int);
    int getNeighbours(int, int, int, int);
    Voxel* getVoxel(boost::tuple<int,int,int> index){
        return space.at(index);
    }
    float getRandomFloat(float);
    int getN_cells(){
        return cells.size();
    };

private:
    std::map<boost::tuple<int,int,int>,Voxel*> space;
    int myrandom (int i) { return std::rand()%i;}
    std::vector<Voxel*> tissue;
   

    
    std::vector<Cell*> cells;    
    //vars for run statistics
    std::vector<std::string> runStatType;
    std::vector<int> runStatTime;
    std::vector<float> runStatVal;
};

#endif	/* ENVIRONMENT_H */

