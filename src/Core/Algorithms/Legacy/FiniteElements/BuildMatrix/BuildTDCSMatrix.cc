/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2009 Scientific Computing and Imaging Institute,
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
   
   Author            : Moritz Dannhauer
   Last modification : July 31 2012
*/



#include <Core/Algorithms/Legacy/FiniteElements/BuildMatrix/BuildTDCSMatrix.h>

#include <Core/Datatypes/Legacy/Field/Mesh.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/Legacy/Matrix/MatrixOperations.h>
#include <Core/Datatypes/Legacy/Matrix/SparseRowMatrixFromMap.h>

#include <Core/GeometryPrimitives/Point.h>
#include <Core/GeometryPrimitives/Tensor.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace SCIRun;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Algorithms::FiniteElements;
using namespace SCIRun::Core::Geometry;
using namespace SCIRun::Core::Utilities;
using namespace SCIRun::Core::Algorithms;

class TDCSMatrixBuilder
{
public:
  TDCSMatrixBuilder(AlgorithmBase* algo) :
  ref_cnt(0),
  algo_(algo),
  electrodes_(1)
  {
      
  }
  SparseRowMatrixHandle getOutput();
  void initialize_mesh(FieldHandle mesh);
  bool initialize_inputs(MatrixHandle stiff, MatrixHandle ElectrodeElements, MatrixHandle ElectrodeElementType, MatrixHandle ElectrodeElementDefinition, MatrixHandle contactimpedance);
  bool build_matrix(MatrixHandle& output);
  bool singlethread();
  int ref_cnt;

private:
  VMesh *mesh_;
  AlgorithmBase* algo_;

  std::vector<unsigned int> electrodes_;
  
  SparseRowMatrixHandle stiffnessMatrix_, tdcs_;
  DenseMatrixHandle electrodeElements_, electrodeElementType_, electrodeElementDefinition_, contactImpedanceInformation_;
  unsigned int electrodeElementsRows_, electrodeElementTypeRows_, electrodeElementDefinitionRows_;
  unsigned int electrodeElementsCols_, electrodeElementTypeCols_, electrodeElementDefinitionCols_;
  unsigned int mesh_nrnodes_, number_electrodes_;

};
    
  
void TDCSMatrixBuilder::initialize_mesh(FieldHandle mesh)
{
  mesh_ = mesh->vmesh(); 
  mesh_nrnodes_=static_cast<unsigned int>(mesh_->num_nodes());
  number_electrodes_=0; 
}


bool  TDCSMatrixBuilder::singlethread()  //single threaded implementation to compute TDCS matrix for point, triangle and tetrahedral electrodes
{
  //the stiffness-matrix updates are written into SparseRowMatrixFromMap::Values additionalData
  //the additional matrices B,C are written into a efficient tree structure (SparseRowMatrixFromMap, thanks to Dan White) and later on added to the modified stiffness matrix to the resulting "tdcs_" - matrix
  //the final tdcs outputmatrix takes the original input stiffness and overwrites it with additionalData
  size_type m = static_cast<size_type>(mesh_nrnodes_);
  size_type n = static_cast<size_type>(mesh_nrnodes_);
  SparseRowMatrixFromMap::Values additionalData;
      
  unsigned int p1=0,p2=0,p3=0,p4=0;
  Point pos; VMesh::Elem::index_type idx; double detJ=0.0;
  double volume=0.0, volume1_4=0.0, volume1_10=0.0, volume1_20=0.0, elc=0.0, surface_impedance=0.0, v_14_imp=0.0, v_110_imp=0.0, v_120_imp=0.0;
  double x1=0.0,x2=0.0,x3=0.0,x4=0.0,y1=0.0,y2=0.0,y3=0.0,y4=0.0,z1=0.0,z2=0.0,z3=0.0,z4=0.0,tmp=0.0,tmp1=0.0,tmp2=0.0,triangle_area=0.0;

 for(index_type i = 0; i<electrodeElementDefinitionRows_; i++)
  {
    if(electrodeElementType_->get(i,0)==1) //point electrodes
    {
      surface_impedance=contactImpedanceInformation_->get(i,0);
      if (surface_impedance<=0){
       algo_->error("Contact surface impedance is negative or zero !");
       return false;
      }
      tmp2=1.0/surface_impedance;
      tmp1=tmp2/2.0;
      p1=static_cast<unsigned int>(electrodeElementDefinition_->get(i,0));
      tmp=stiffnessMatrix_->get(p1,p1);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); return false; 
      }

      if (tmp==0) 
      { 
//       std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p1 << "," << p1 << ")" << std::endl;  
//       algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }

      tmp+=tmp1; additionalData[p1][p1]+=tmp;      
      elc=electrodeElements_->get(i,0)+mesh_nrnodes_;
      SparseRowMatrixFromMap::Row& rowElc = additionalData[elc];
      rowElc[p1] -= tmp1;
      additionalData[p1][elc] -= tmp1;
      rowElc[elc] += tmp2*0.5;
    } else
    if(electrodeElementType_->get(i,0)==2) //electrode made of triangles
    {
      // compute triangle surface area
      mesh_->get_point(pos, VMesh::Node::index_type(electrodeElementDefinition_->get(i,0)));
      x1=pos.x(); y1=pos.y(); z1=pos.z();
      
      mesh_->get_point(pos, VMesh::Node::index_type(electrodeElementDefinition_->get(i,1)));
      x2=pos.x(); y2=pos.y(); z2=pos.z();
      
      mesh_->get_point(pos, VMesh::Node::index_type(electrodeElementDefinition_->get(i,2)));
      x3=pos.x(); y3=pos.y(); z3=pos.z();
      
      //compute triangle area in 3D
      tmp =y1*z2+z1*y3+y2*z3-z2*y3-z1*y2-y1*z3;
      tmp1=z1*x2+x1*z3+z2*x3-x2*z3-x1*z2-z1*x3;
      tmp2=x1*y2+y1*x3+x2*y3-y2*x3-y1*x2-x1*y3;
      
      triangle_area=0.5 * sqrt(tmp*tmp+tmp1*tmp1+tmp2*tmp2);

      if(triangle_area<=0) { 
       algo_->error(" Triangle area should be positive! "); 
       return false;
      }
      
      surface_impedance=contactImpedanceInformation_->get(i,0);
      if (surface_impedance<=0)
      {
       algo_->error("Contact surface impedance is negative or zeros !");
       return false;
      }
      tmp1 = (2.0*triangle_area) / surface_impedance;
      
      p1=static_cast<unsigned int>(electrodeElementDefinition_->get(i,0));
      p2=static_cast<unsigned int>(electrodeElementDefinition_->get(i,1));
      p3=static_cast<unsigned int>(electrodeElementDefinition_->get(i,2));
      tmp=tmp1/6.0; //prepare value for B and Bt
      
      elc=electrodeElements_->get(i,0)+mesh_nrnodes_;
      SparseRowMatrixFromMap::Row& rowElc = additionalData[elc];
      //Bt
      rowElc[p1] -= tmp;
      rowElc[p2] -= tmp;
      rowElc[p3] -= tmp;
      
      //B
      additionalData[p1][elc] = rowElc[p1];
      additionalData[p2][elc] = rowElc[p2];
      additionalData[p3][elc] = rowElc[p3];

      tmp=tmp1/12.0;
      tmp2=stiffnessMatrix_->get(p1,p1);
      if(IsNan(tmp2) || !IsFinite(tmp2) || IsInfinite(tmp2)) { 
       algo_->error(" No valid value "); 
       return false; 
      }
      if (tmp2==0) {
  //     std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p1 << "," << p1 << ")" << std::endl;  
  //     algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp2+=tmp; additionalData[p1][p1]=tmp2;

      tmp2=stiffnessMatrix_->get(p2,p2);
      if(IsNan(tmp2) || !IsFinite(tmp2) || IsInfinite(tmp2)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
    //   std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p2 << "," << p2 << ")" << std::endl;  
    //   algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }     
      tmp2+=tmp; additionalData[p2][p2]=tmp2;
      
      tmp2=stiffnessMatrix_->get(p3,p3);
      if(IsNan(tmp2) || !IsFinite(tmp2) || IsInfinite(tmp2)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
    //   std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p3 << "," << p3 << ")" << std::endl;  
    //   algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp2+=tmp; additionalData[p3][p3]=tmp2;
      
      tmp=tmp1/24.0;
      tmp2=stiffnessMatrix_->get(p1,p2);
      if(IsNan(tmp2) || !IsFinite(tmp2) || IsInfinite(tmp2)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
  //     std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p1 << "," << p2 << ")" << std::endl;  
  //     algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp2+=tmp; additionalData[p1][p2]=tmp2;
      
      tmp2=stiffnessMatrix_->get(p1,p3);
      if(IsNan(tmp2) || !IsFinite(tmp2) || IsInfinite(tmp2)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p1 << "," << p3 << ")" << std::endl;  
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp2+=tmp; additionalData[p1][p3]=tmp2;
      
      tmp2=stiffnessMatrix_->get(p2,p1);
      if(IsNan(tmp2) || !IsFinite(tmp2) || IsInfinite(tmp2)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
  //     std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p2 << "," << p1 << ")" << std::endl;  
  //     algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp2+=tmp; additionalData[p2][p1]=tmp2;
      
      tmp2=stiffnessMatrix_->get(p2,p3);
      if(IsNan(tmp2) || !IsFinite(tmp2) || IsInfinite(tmp2)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
  //     std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p2 << "," << p3 << ")" << std::endl;  
  //     algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp2+=tmp; additionalData[p2][p3]=tmp2;
      
      tmp2=stiffnessMatrix_->get(p3,p1);
      if(IsNan(tmp2) || !IsFinite(tmp2) || IsInfinite(tmp2)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p3 << "," << p1 << ")" << std::endl;  
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp2+=tmp; additionalData[p3][p1]=tmp2;
      
      tmp2=stiffnessMatrix_->get(p3,p2);
      if(IsNan(tmp2) || !IsFinite(tmp2) || IsInfinite(tmp2)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p3 << "," << p2 << ")" << std::endl;  
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp2+=tmp; additionalData[p3][p2]=tmp2;
     
      additionalData[elc][elc] += 0.5*tmp1;      
    } else
    if(electrodeElementType_->get(i,0)==3)  //electrode made of tetrahedral elements
    {
      mesh_->get_point(pos, VMesh::Node::index_type(electrodeElementDefinition_->get(i,0)));
      x1=pos.x(); y1=pos.y(); z1=pos.z();
        
      mesh_->get_point(pos, VMesh::Node::index_type(electrodeElementDefinition_->get(i,1)));
      x2=pos.x(); y2=pos.y(); z2=pos.z();
        
      mesh_->get_point(pos, VMesh::Node::index_type(electrodeElementDefinition_->get(i,2)));
      x3=pos.x(); y3=pos.y(); z3=pos.z();
       
      mesh_->get_point(pos, VMesh::Node::index_type(electrodeElementDefinition_->get(i,3)));
      x4=pos.x(); y4=pos.y(); z4=pos.z();
        
      //  compute determinant of jacobian which is needed for volume of tet	
      detJ=(x3*(y2*z1-y1*z2)+ x2*(y1*z3-y3*z1) - x1*(y2*z3-y3*z2))+(-x4*(y2*z1-y1*z2)-x2* ( y1*z4 - y4*z1)+x1* ( y2*z4 - y4*z2))-(-x4*(y3*z1-y1*z3)-x3*(y1*z4-y4*z1)+x1*( y3*z4-y4*z3))+(-x4* (y3*z2-y2*z3)-x3* ( y2*z4-y4*z2)+x2* (y3*z4-y4*z3));      
      
      if (detJ<=0)
      {
        tdcs_ = SparseRowMatrixFromMap::appendToSparseMatrix(m+number_electrodes_, n+number_electrodes_, *stiffnessMatrix_, additionalData);
        algo_->error("Mesh has elements with negative/zero jacobians, check the order of the nodes that define an element");
        return false;
      }
      volume=1.0/6.0*detJ; volume1_4 =volume/4.0; volume1_10=volume/10.0; volume1_20=volume/20.0;
      
      surface_impedance=contactImpedanceInformation_->get(i,0);
      if (surface_impedance<=0)
      {
        tdcs_ = SparseRowMatrixFromMap::appendToSparseMatrix(m+number_electrodes_, n+number_electrodes_, *stiffnessMatrix_, additionalData);
        algo_->error("Contact surface impedance is negative or zeros !");
        return false;
      }
      
      p1=static_cast<unsigned int>(electrodeElementDefinition_->get(i,0));
      p2=static_cast<unsigned int>(electrodeElementDefinition_->get(i,1));
      p3=static_cast<unsigned int>(electrodeElementDefinition_->get(i,2));
      p4=static_cast<unsigned int>(electrodeElementDefinition_->get(i,3));
      elc=electrodeElements_->get(i,0)+mesh_nrnodes_;
      v_14_imp=static_cast<double>(-volume1_4/surface_impedance);
      v_110_imp=static_cast<double>(volume1_10/surface_impedance);
      v_120_imp=static_cast<double>(volume1_20/surface_impedance);
      
      // B^t
      SparseRowMatrixFromMap::Row& rowElc = additionalData[elc];
      rowElc[p1] += v_14_imp;
      rowElc[p2] += v_14_imp;
      rowElc[p3] += v_14_imp;
      rowElc[p4] += v_14_imp;
      // B
      additionalData[p1][elc] = rowElc[p1];
      additionalData[p2][elc] = rowElc[p2];
      additionalData[p3][elc] = rowElc[p3];
      additionalData[p4][elc] = rowElc[p4];
      
      //A
      tmp=stiffnessMatrix_->get(p1,p1);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p1 << "," << p1 << ")" << std::endl;  
//       algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      } 
      tmp+=v_110_imp; additionalData[p1][p1]=tmp;
      
      tmp=stiffnessMatrix_->get(p2,p2);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p2 << "," << p2 << ")" << std::endl;  
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_110_imp; additionalData[p2][p2]=tmp;   
      
      tmp=stiffnessMatrix_->get(p3,p3);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      
      if (tmp2==0) {
  //     std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p3 << "," << p3 << ")" << std::endl;  
  //     algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_110_imp; additionalData[p3][p3]=tmp;
      
      tmp=stiffnessMatrix_->get(p4,p4);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); return false; 
      }
      
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p4 << "," << p4 << ")" << std::endl; 
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      } 
      tmp+=v_110_imp; additionalData[p4][p4]=tmp; 
      
      tmp=stiffnessMatrix_->get(p1,p2);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p1 << "," << p2 << ")" << std::endl;  
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_120_imp; additionalData[p1][p2]=tmp;
      
      tmp=stiffnessMatrix_->get(p2,p1);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p2 << "," << p1 << ")" << std::endl;  
//       algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_120_imp; additionalData[p2][p1]=tmp;
      
      tmp=stiffnessMatrix_->get(p1,p3);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      
      if (tmp2==0) {
//       std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p1 << "," << p3 << ")" << std::endl;  
//       algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_120_imp; additionalData[p1][p3]=tmp;
      tmp=stiffnessMatrix_->get(p3,p1);
      
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      
      if (tmp2==0) {
//       std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p3 << "," << p1 << ")" << std::endl;  
//       algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_120_imp; additionalData[p3][p1]=tmp;
      
      tmp=stiffnessMatrix_->get(p1,p4);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p1 << "," << p4 << ")" << std::endl;  
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      } 
      tmp+=v_120_imp; additionalData[p1][p4]=tmp;
      
      tmp=stiffnessMatrix_->get(p4,p1);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
//       std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p4 << "," << p1 << ")" << std::endl;  
//       algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_120_imp; additionalData[p4][p1]=tmp;
      
      tmp=stiffnessMatrix_->get(p2,p3);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p2 << "," << p3 << ")" << std::endl;  
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_120_imp; additionalData[p2][p3]=tmp;      
      
      tmp=stiffnessMatrix_->get(p3,p2);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p3 << "," << p2 << ")" << std::endl;  
//       algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_120_imp; additionalData[p3][p2]=tmp;
      
      tmp=stiffnessMatrix_->get(p2,p4);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }      
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p2 << "," << p4 << ")" << std::endl;  
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_120_imp; additionalData[p2][p4]=tmp;
      
      tmp=stiffnessMatrix_->get(p4,p2);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p4 << "," << p2 << ")" << std::endl;  
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_120_imp; additionalData[p4][p2]=tmp;
            
      tmp=stiffnessMatrix_->get(p3,p4);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); 
       return false; 
      }
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p3 << "," << p4 << ")" << std::endl;  
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_120_imp; additionalData[p3][p4]=tmp;
      
      tmp=stiffnessMatrix_->get(p4,p3);
      if(IsNan(tmp) || !IsFinite(tmp) || IsInfinite(tmp)) { 
       algo_->error(" No valid value in stiffnessmatrix !"); return false; 
      }
      if (tmp2==0) {
 //      std::cout<<"Unexpected zero entry in stiffness matrix at ("<< p4 << "," << p3 << ")" << std::endl;  
 //      algo_->remark("At least one entry in the stiffness matrix for an element (see console output) is zero. Please review the used units (conductivity, mesh location) and scale them up.");
      }
      tmp+=v_120_imp; additionalData[p4][p3]=tmp;
      
      //C
      rowElc[elc] += volume/surface_impedance;
    } else
   {
    tdcs_ = SparseRowMatrixFromMap::appendToSparseMatrix(m+number_electrodes_, n+number_electrodes_, *stiffnessMatrix_, additionalData);
    algo_->error(" This Electrode-type is not implemented, a electrode only consist of points(1), triangles(2) and tetrahedral elements(3)."); return false;
  }
  }
  tdcs_ = SparseRowMatrixFromMap::appendToSparseMatrix(m+number_electrodes_, n+number_electrodes_, *stiffnessMatrix_, additionalData);
    
  return true;
}  
    
SparseRowMatrixHandle TDCSMatrixBuilder::getOutput()
{
  return tdcs_;  
}

bool TDCSMatrixBuilder::initialize_inputs(MatrixHandle stiff, MatrixHandle ElectrodeElements, MatrixHandle ElectrodeElementType, MatrixHandle ElectrodeElementDefinition, MatrixHandle contactimpedance)
{
  electrodeElements_=(ElectrodeElements)->dense();
  electrodeElementType_=(ElectrodeElementType)->dense();
  electrodeElementDefinition_=(ElectrodeElementDefinition)->dense();
  
  //check matrices dimensions
  electrodeElementsRows_=static_cast<unsigned int>(electrodeElements_->nrows());
  electrodeElementTypeRows_=static_cast<unsigned int>(electrodeElementType_->nrows());
  electrodeElementDefinitionRows_=static_cast<unsigned int>(electrodeElementDefinition_->nrows());
  electrodeElementsCols_=static_cast<unsigned int>(electrodeElements_->ncols());
  electrodeElementTypeCols_=static_cast<unsigned int>(electrodeElementType_->ncols());
  electrodeElementDefinitionCols_=static_cast<unsigned int>(electrodeElementDefinition_->ncols());

  if( !((electrodeElementsRows_==electrodeElementTypeRows_) && (electrodeElementTypeRows_==electrodeElementDefinitionRows_)) ) {
    algo_->error("Number of Matrix-rows of Matrices for Electrode-Definition should be the same!");
    return (false);
  }
  
  if( !((electrodeElementsCols_==1) && (electrodeElementTypeCols_==1) && (electrodeElementDefinitionCols_>=1 && electrodeElementDefinitionCols_<=4) ) ) {
    algo_->error(" Number of Matrix-columns of Matrices in Electrode-Definition should be: ElectrodeElementType=1, ElectrodeElementDefinition=1,ElectrodeElementDefinition=4 !");
    return (false);
  }
  
  if (electrodeElementsRows_>mesh_nrnodes_) {
    algo_->error(" Number of Electrode-Definition Nodes can not excced number of input mesh (since it refers to that) !");
    return (false);
  } 
 
  if ( electrodeElementsRows_==0 ) {
    algo_->error(" Number of Electrode-Definition Nodes = 0 !");
    return (false);
  }

  //check matrices content
  electrodes_.push_back(electrodeElements_->get(0,0));
  number_electrodes_=1; bool numbering_ok=false;
  for (unsigned int i=0;i<electrodeElementsRows_;i++)
  {
    bool found = false;
    if (static_cast<unsigned int>(electrodeElements_->get(0,0))==0)  {
      numbering_ok=true;
    }
    
    unsigned int tmp=static_cast<unsigned int>(electrodeElements_->get(i,0));
    if( (tmp<0) || (tmp>=mesh_nrnodes_) ) {
      algo_->error(" Specified ElectrodeElement node out of range (0..#inputmeshnodes-1) !");
      return (false);
    }
    
    unsigned int tmp1=static_cast<unsigned int>(electrodeElementType_->get(i,0));
    if ( !(tmp1>=1 && tmp1<=3) ) {
      algo_->error(" Specified ElectrodeElementType out of range, allowed range: 1 (point), 2 (triangle), 3 (tetrahedra) !");
      return (false);
    }
    
    for(unsigned int j=0;j<electrodeElementDefinitionCols_;j++)
    {
      unsigned int tmp2=static_cast<unsigned int>(electrodeElementDefinition_->get(i,j));
      if ( (tmp2<0) || (tmp2>=mesh_nrnodes_) ) {
        algo_->error(" Specified ElectrodeElementDefinition is out of range (> mesh nodes-1) - allowed range is (allowed range: 0..#meshnodes-1) !");
        return (false);
      }  
    }
    
    for (unsigned int j=0;j<electrodes_.size();j++)
    {
      if (electrodes_[j]==tmp)
      {
        found = true;
        break;
      }
    }

    if (! found)
    {
      electrodes_.push_back(electrodeElements_->get(i,0));
      number_electrodes_++;
    }
  }
  
  if (!numbering_ok)
  {
   algo_->error(" The electrode numbering should should start at 0 (allowed range: 0..#meshnodes-1) !");
   return (false);
  }
  
   // get surface impedance
  if(contactimpedance.get_rep())
  {
    contactImpedanceInformation_ = contactimpedance->dense();
    if(static_cast<unsigned int>(contactImpedanceInformation_->nrows())!=electrodeElementsRows_) {
     algo_->error(" Contact surface impedance vector and electrode definition does not fit !"); 
     return (false);
    }
  }
  else 
  { 
    contactImpedanceInformation_=new DenseMatrix(electrodeElementsRows_,1);  
    for(unsigned int i=0;i<electrodeElementsRows_;i++) 
      contactImpedanceInformation_->put(i,0,1); 
  }

  if (  !( (contactImpedanceInformation_->nrows()==electrodeElementsRows_) && (contactImpedanceInformation_->ncols()==1))  ) {
    algo_->error(" ContactImpedanceMatrix should have matrix dimensions: #electrodeselementsx1 !");
    return (false);
  } 
  
  stiffnessMatrix_ = stiff->sparse();
  if ( !( (stiffnessMatrix_->nrows()>0) && (stiffnessMatrix_->nrows()==stiffnessMatrix_->ncols()))) {
    algo_->error(" StiffnessMatrix should be square and non-empty !");
    return (false);
  }
  if ( !( stiffnessMatrix_->nrows()==mesh_nrnodes_  ) ) {
    algo_->error(" StiffnessMatrix should have the same number of nodes as inputmesh !");
    return (false);
  }
  if (! stiffnessMatrix_->validate()) {
    algo_->error(" StiffnessMatrix is broken  !"); 
    return (false);
  }
  
  return true;
}
  

bool TDCSMatrixBuilder::build_matrix(MatrixHandle& output)
{
  return singlethread();
}

bool BuildTDCSMatrix::run(MatrixHandle stiff, FieldHandle mesh, MatrixHandle ElectrodeElements, MatrixHandle ElectrodeElementType, MatrixHandle ElectrodeElementDefinition, MatrixHandle contactimpedance, MatrixHandle& output)
{
  algo_start("TDCSMatrixBuilder");
    
  if (! mesh.get_rep())
  {
    error(" Without a mesh there is nothing to do !");
    algo_end();
    return (false);   
  }
 
  if (! ElectrodeElements.get_rep()) //get electrode definition
  {
    error("ElectrodeElements object not available....");
    algo_end();
    return false;
  }   
   
  if (! ElectrodeElementType.get_rep()) //get electrode definition
  {
    error("ElectrodeElements object not available....");
    algo_end();
    return false;
  }  
   
  if (! ElectrodeElementDefinition.get_rep()) //get electrode definition
  {
    error("ElectrodeElements object not available....");
    algo_end();
    return false;
  }  
 
  Handle<TDCSMatrixBuilder> builder = new TDCSMatrixBuilder(this);

  builder->initialize_mesh(mesh); //set mesh 

  if (! builder->initialize_inputs(stiff, ElectrodeElements, ElectrodeElementType, ElectrodeElementDefinition, contactimpedance) ) // set other inputs
  {
    algo_end();
    return false;
  }

  if (! builder->build_matrix(output))
  {
    algo_end();
    return false;
  }

  output=builder->getOutput();
  
  SparseRowMatrix *tdcs=output->sparse();
  
  bool found_reference_node=false;

  for(index_type i = 0; i<tdcs->nrows()-1; i++) //find reference node 
  {
    index_type ps = (*tdcs).get_row(i);
    index_type pe = (*tdcs).get_row(i+1);
    if((pe-ps)==1)
       {
        if (static_cast<double>((*tdcs).get_value(ps))==1.0)
        found_reference_node=true;
	break;
       }
  }
   
  if(!found_reference_node) 
      remark("The TDCS output matrix is not referenced yet !!! Please set the Potential of at least one node to 0 ! You can do that by setting one row and the corresponding column explicitely to 0 except for the diagonal element which should be 1 !");

 algo_end(); 
  
 return (true);
}

BuildTDCSMatrix::BuildTDCSMatrix() {}
BuildTDCSMatrix::~BuildTDCSMatrix() {}