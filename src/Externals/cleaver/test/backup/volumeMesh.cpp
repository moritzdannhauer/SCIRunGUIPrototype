/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2020 Scientific Computing and Imaging Institute,
   University of Utah.

   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/


#include <iostream>
#include <vector>
#include <string>

#include <Cleaver/Cleaver.h>
#include <Cleaver/InverseField.h>
#include <Cleaver/FloatField.h>
//#include "nrrd2cleaver.h"


using namespace std;
//using namespace Cleaver;

const std::string scirun = "scirun";
const std::string tetgen = "tetgen";
const std::string matlab = "matlab";


int main(int argc, char *argv[])
{


	cerr<<"THIS IS VOLUMEMESH"<<endl;

	int dim_x=30;
	int dim_y=30;
	int dim_z=30;
	float v_x=1.0;
	float v_y=1.0;
	float v_z=1.0;

	bool verbose=true;
	string format="tetgen";
	string outputFileName="meshOutput";

	//vector<float> outsideVec(dim_x*dim_y*dim_x,0);
	vector<float> insideVec(dim_x*dim_y*dim_x,0.0);
	int idx=0;
	//filling volumetric mateiral definition - a sphere centered in the middle of the lattice
	for (int x = 0  ; x < dim_x ; ++x ){
	  for (int y = 0  ; y < dim_y ; ++y ){
	    for (int z = 0  ; z < dim_z ; ++z ){
	      if(x < 2 || y < 2 || z < 2 ||
		 x > (dim_x - 2) || y > (dim_y - 2) || z > (dim_z - 2))
		{
		  insideVec[idx]=1.0;
		  ++idx;
			      continue;
		}

	      if ((x-dim_x/2)*(x-dim_x/2)+(y-dim_y/2)*(y-dim_y/2)+(z-dim_z/2)<0.4*dim_x){

		insideVec[idx]=1.0;
	      }else{
		insideVec[idx]=0.0;
		//insideVec[idx]=1;
		//outsideVec[idx]=0;
	      }

	      std::cout << insideVec[idx] << std::endl;

	      ++idx;
	    }
	  }
	}

	Cleaver::FloatField insideField=Cleaver::FloatField(dim_x,dim_y,dim_z,&insideVec[0]);
	//Cleaver::FloatField outsideField=Cleaver::FloatField(dim_x,dim_y,dim_z,&outsideVec[0]);

	Cleaver::InverseField inverseField=Cleaver::InverseField(&insideField);

	std::vector<Cleaver::ScalarField*> fields;
	fields.push_back(&insideField);
	fields.push_back(&inverseField);

	Cleaver::Volume volume(fields);
	Cleaver::TetMesh *mesh = Cleaver::createMeshFromVolume(volume, verbose);

    //------------------
    //  Compute Angles
    //------------------
    mesh->computeAngles();
    if(verbose){
        std::cout.precision(12);
        std::cout << "Worst Angles:" << std::endl;
        std::cout << "min: " << mesh->min_angle << std::endl;
        std::cout << "max: " << mesh->max_angle << std::endl;
    }


    //----------------------
    //  Write Info File
    //----------------------
    mesh->writeInfo(outputFileName, verbose);


    //----------------------
    // Write Tet Mesh Files
    //----------------------
    if(format == tetgen)
        mesh->writeNodeEle(outputFileName, verbose);
    else if(format == scirun)
        mesh->writePtsEle(outputFileName, verbose);
    else if(format == matlab)
        mesh->writeMatlab(outputFileName, verbose);

    //----------------------
    // Write Surface Files
    //----------------------
    mesh->constructFaces();
    mesh->writePly(outputFileName, verbose);



}
