/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2012 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
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

#include <Core/Algorithms/Base/AlgorithmPreconditions.h>
#include <Core/Algorithms/BrainStimulator/GenerateROIStatisticsAlgorithm.h>
#include <Core/GeometryPrimitives/Vector.h>
#include <Core/Datatypes/Legacy/Field/Field.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
#include <Core/Datatypes/Legacy/Field/VMesh.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/String.h>
//////////////////////////////////////////////////////////////////////////
/// @todo MORITZ
//////////////////////////////////////////////////////////////////////////
#include <iostream>
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Core::Algorithms::BrainStimulator;
using namespace SCIRun::Core::Geometry;
using namespace SCIRun;
    
const AlgorithmInputName GenerateROIStatisticsAlgorithm::MESH_DATA_ON_ELEMENTS("MESH_DATA_ON_ELEMENTS");
const AlgorithmInputName GenerateROIStatisticsAlgorithm::PHYSICAL_UNIT("PHYSICAL_UNIT");
const AlgorithmInputName GenerateROIStatisticsAlgorithm::ATLAS_MESH("ATLAS_MESH");
const AlgorithmInputName GenerateROIStatisticsAlgorithm::ATLAS_MESH_LABELS("ATLAS_MESH_LABELS");
const AlgorithmInputName GenerateROIStatisticsAlgorithm::COORDINATE_SPACE("COORDINATE_SPACE");
const AlgorithmOutputName GenerateROIStatisticsAlgorithm::STATISTICAL_RESULTS("STATISTICAL_RESULTS");


DenseMatrixHandle GenerateROIStatisticsAlgorithm::run(FieldHandle mesh, FieldHandle atlas_mesh) const
{
 DenseMatrixHandle output;
 
 return output;
}

AlgorithmOutput GenerateROIStatisticsAlgorithm::run_generic(const AlgorithmInput& input) const
{
  auto mesh = input.get<Field>(MESH_DATA_ON_ELEMENTS);
  auto physical_unit = input.get<Datatypes::String>(PHYSICAL_UNIT);
  auto atlas_mesh = input.get<Field>(ATLAS_MESH);
  auto atlas_mesh_labels = input.get<Datatypes::String>(ATLAS_MESH_LABELS);
  auto coordinate = input.get<Field>(COORDINATE_SPACE);
  
  DenseMatrixHandle statistics = run(mesh, atlas_mesh);
  
 // getOptionalInput
 /* auto pos_orient = input.get<Field>(MESH_DATA_ON_ELEMENTS);
  auto tri = input.get<Field>(PHYSICAL_UNIT);
  auto tri2 = input.get<Field>(ELECTRODE_TRIANGULATION2);
  auto coil = input.get<Field>(COIL);
  auto COORDINATE_SPACE = input.get<Field>(COORDINATE_SPACE);*/
 
  //old-style run call, just put algorithm code here
  //auto outputs = run(boost::make_tuple(lhs, rhs), Option(get(Variables::AppendMatrixOption).getInt()));
  // CODE HERE
  FieldHandle out1;
  Datatypes::String out2,out3;

  //Algorithm starts here:
  //VField* vfield = elc_coil_pos_and_normal->vfield();
  // VMesh*  vmesh  = pos_orient->vmesh();
 
   //std::cout << "a: " << vmesh->num_nodes() << std::endl;
   //for (int i=0;i<vmesh->num_nodes();;i++)
   //{
   
   
   //}
  //


  AlgorithmOutput output;
  output[STATISTICAL_RESULTS] = statistics;
  
  return output;
}
