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
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>
#include <Core/Datatypes/Legacy/Field/VMesh.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/String.h>
#include <boost/range/algorithm/count.hpp>
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
const AlgorithmInputName GenerateROIStatisticsAlgorithm::COORDINATE_SPACE_LABEL("COORDINATE_SPACE_LABEL");
const AlgorithmOutputName GenerateROIStatisticsAlgorithm::STATISTICAL_RESULTS("STATISTICAL_RESULTS");

GenerateROIStatisticsAlgorithm::GenerateROIStatisticsAlgorithm()
{
 
}

DenseMatrixHandle GenerateROIStatisticsAlgorithm::run(FieldHandle mesh, FieldHandle atlas_mesh) const
{
 DenseMatrixHandle output;
 
 VField* vfield1 = mesh->vfield();
 VField* vfield2 = atlas_mesh->vfield();
 
 std::vector<int>  label_vector(vfield2->vmesh()->num_elems());
 std::vector<double> value_vector(vfield1->vmesh()->num_elems());
 
 for (VMesh::Elem::index_type i=0; i < vfield2->vmesh()->num_elems(); i++) // loop over all tetrahedral elements (mesh)
 {
   int label = 0;
   vfield2->get_value(label, i);
   label_vector[i]=label;
   double val = 0;
   vfield1->get_value(val, i);
   value_vector[i]=val;
 }
  
 std::vector<int>::iterator it;
 it = std::unique (label_vector.begin(), label_vector.end());
 label_vector.resize( std::distance(label_vector.begin(),it) );
 
 long number_of_atlas_materials = label_vector.size();
  
 std::vector<double> value_avr(number_of_atlas_materials);
 std::vector<int> value_count(number_of_atlas_materials);
 std::vector<double> value_min(number_of_atlas_materials);
 std::vector<double> value_max(number_of_atlas_materials);
 std::vector<double> value_std(number_of_atlas_materials);
 std::vector<double> Sxsqr(number_of_atlas_materials);
 std::vector<double> stddev(number_of_atlas_materials);
 
 for (VMesh::Elem::index_type i=0; i < vfield1->vmesh()->num_elems(); i++) // loop over all tetrahedral elements (atlas_mesh)
 {
   double value = 0;
   vfield1->get_value(value, i); 
   
   int label = 0;
   vfield2->get_value(label, i);
   
   for (VMesh::Elem::index_type j=0; j < number_of_atlas_materials; j++)
   {
     if (label==label_vector[j]) 
        {
	   value_avr[j] += value; 
	   value_count[j]++;   
	   if (value>value_max[j]) value_max[j]=value;
	   if (value<value_min[j]) value_min[j]=value;
	   Sxsqr[j]+=value*value;
	}
   }
 }
 
 output = DenseMatrixHandle(new DenseMatrix(number_of_atlas_materials, 4));
 
 //efficient way to compute std dev. in just one loop over all mesh elements: sqrt ( 1/(n-1) (Sx^2 - avr Sx + n avr^2 )
 for (VMesh::Elem::index_type j=0; j < number_of_atlas_materials; j++)
 {
   double Sx=value_avr[j];
   value_avr[j]/=value_count[j]; 
   stddev[j]=sqrt(1/(value_count[j]-1)*(Sxsqr[j]-2*value_avr[j]*Sx+value_count[j]*value_avr[j]*value_avr[j]));
   (*output)(j,0)=value_avr[j];
   (*output)(j,1)=stddev[j];
   (*output)(j,2)=value_min[j];
   (*output)(j,3)=value_max[j];
 }
  
 return output;
}

std::vector<std::string> GenerateROIStatisticsAlgorithm::ConvertInputAtlasStringIntoVector(const std::string atlas_labels) const
{
  std::string s = atlas_labels;

  int cnt = boost::count(atlas_labels, ';')+1;  /// Place a ";" after each region name, e.g. ROI1;ROI2; 

  std::vector<std::string> result;
  result.reserve(cnt);
  std::string delimiter = ";";

  size_t pos = 0;
  std::string token;   

  while ((pos = s.find(delimiter)) != std::string::npos) 
  {
    token = s.substr(0, pos);
    result.push_back(token);
    s.erase(0, pos + delimiter.length());
  }  
    
  return result;
}

AlgorithmOutput GenerateROIStatisticsAlgorithm::run_generic(const AlgorithmInput& input) const
{
  auto mesh = input.get<Field>(MESH_DATA_ON_ELEMENTS);
  auto physical_unit_ = input.get<Datatypes::String>(PHYSICAL_UNIT);
  auto atlas_mesh = input.get<Field>(ATLAS_MESH);
  auto atlas_mesh_labels = (input.get<Datatypes::String>(ATLAS_MESH_LABELS))->value();
  auto coordinate = input.get<Field>(COORDINATE_SPACE);
  auto coordinate_label = input.get<Datatypes::String>(COORDINATE_SPACE_LABEL);
  
  if (!mesh)  
     THROW_ALGORITHM_INPUT_ERROR("First input (mesh) is empty.");
  
  if (!atlas_mesh)  
     THROW_ALGORITHM_INPUT_ERROR("Third input (atlas mesh) is empty.");
  
  FieldInformation fi(mesh);
  
  if (!fi.is_constantdata())
    THROW_ALGORITHM_INPUT_ERROR("First input (mesh) requires the data to be on the elements.");
  
  // making sure the field contains data
  VField* vfield1 = mesh->vfield();
  if (vfield1->is_nodata())
    THROW_ALGORITHM_INPUT_ERROR("First input field (mesh) contained no data.");
  
  // making sure the field is not in vector format
  if (!vfield1->is_scalar())
    THROW_ALGORITHM_INPUT_ERROR("First input field needs to have scalar data.");      
  
  FieldInformation fi2(atlas_mesh);
  
  if (!fi2.is_constantdata())
    THROW_ALGORITHM_INPUT_ERROR("First input (mesh) requires the data to be on the elements.");
  
  // making sure the field contains data
  VField* vfield2 = atlas_mesh->vfield();
  if (vfield2->is_nodata())
    THROW_ALGORITHM_INPUT_ERROR("First input field (mesh) contained no data.");
  
  // making sure the field is not in vector format
  if (!vfield2->is_scalar())
    THROW_ALGORITHM_INPUT_ERROR("First input field needs to have scalar data."); 
    
  if(vfield1->vmesh()->num_elems()<1 && vfield2->vmesh()->num_elems()<1)
    THROW_ALGORITHM_INPUT_ERROR("First (mesh) or second (atlas_mesh) input field does not contain elements."); 
  
  if(vfield2->vmesh()->num_elems() !=  vfield1->vmesh()->num_elems())
    THROW_ALGORITHM_INPUT_ERROR(" Number of mesh elements of first input and third input does not match.");  
  
  std::vector<std::string> atlas_mesh_labels_vector;
  if (!atlas_mesh_labels.empty())
  {
   atlas_mesh_labels_vector = ConvertInputAtlasStringIntoVector(atlas_mesh_labels); 
  }
   
  DenseMatrixHandle statistics = run(mesh, atlas_mesh);

  AlgorithmOutput output;
  output[STATISTICAL_RESULTS] = statistics;
  
  return output;
}
