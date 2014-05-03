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

#include <Modules/Math/ConvertMatrixType.h>
#include <Core/Algorithms/Math/ConvertMatrixType.h>
#include <Core/Datatypes/Matrix.h>

using namespace SCIRun::Modules::Math;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Core::Algorithms::Math;
using namespace SCIRun::Dataflow::Networks;
using namespace SCIRun::Core::Datatypes;

ConvertMatrixTypeModule::ConvertMatrixTypeModule() : Module(ModuleLookupInfo("ConvertMatrixType", "Math", "SCIRun")) 
{
  INITIALIZE_PORT(InputMatrix);
  INITIALIZE_PORT(ResultMatrix);
}

void ConvertMatrixTypeModule::setStateDefaults()
{
 setStateBoolFromAlgo(ConvertMatrixTypeAlgorithm::PassThrough);
 setStateBoolFromAlgo(ConvertMatrixTypeAlgorithm::Convert2ColumnMatrix);
 setStateBoolFromAlgo(ConvertMatrixTypeAlgorithm::Convert2DenseMatrix);
 setStateBoolFromAlgo(ConvertMatrixTypeAlgorithm::Convert2SparseRowMatrix);
}



void ConvertMatrixTypeModule::execute()
{
 
  auto input_matrix = getRequiredInput(InputMatrix);
  algo().set(ConvertMatrixTypeAlgorithm::PassThrough,get_state()->getValue(ConvertMatrixTypeAlgorithm::PassThrough).getBool());
  algo().set(ConvertMatrixTypeAlgorithm::Convert2ColumnMatrix,get_state()->getValue(ConvertMatrixTypeAlgorithm::Convert2ColumnMatrix).getBool());  
  algo().set(ConvertMatrixTypeAlgorithm::Convert2DenseMatrix,get_state()->getValue(ConvertMatrixTypeAlgorithm::Convert2DenseMatrix).getBool());  
  algo().set(ConvertMatrixTypeAlgorithm::Convert2SparseRowMatrix,get_state()->getValue(ConvertMatrixTypeAlgorithm::Convert2SparseRowMatrix).getBool());  
  auto output = algo().run_generic(make_input((InputMatrix, input_matrix)));
 
  sendOutputFromAlgorithm(ResultMatrix, output);
}