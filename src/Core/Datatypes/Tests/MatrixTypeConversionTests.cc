/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2012 Scientific Computing and Imaging Institute,
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
   
   author: Moritz Dannhauer
   last change: 04/21/14
*/

#include <gtest/gtest.h>
#include <Core/Datatypes/MatrixTypeConversions.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixComparison.h>
#include <Core/Datatypes/MatrixTypeConversions.h>

using namespace SCIRun;
using namespace SCIRun::Core::Datatypes;

DenseColumnMatrixHandle CreateColumnMatrix()  
{
    DenseColumnMatrixHandle m(boost::make_shared<DenseColumnMatrix>(3));

    (*m)(0) = 1;
    (*m)(1) = 2;
    (*m)(2) = 3;

    return m;
}

SparseRowMatrixHandle CreateSparseMatrixWithOneColumn()  
{ 
    SparseRowMatrixHandle m(boost::make_shared<SparseRowMatrix>(3,1));
    m->insert(0,0) = 1;
    m->insert(1,0) = 2;
    m->insert(2,0) = 3;
    m->makeCompressed();
    return m;
}

DenseColumnMatrixHandle CreateColumnMatrix_2()  
{
    DenseColumnMatrixHandle m(boost::make_shared<DenseColumnMatrix>(3));

    (*m)(0) = 1;
    (*m)(1) = 0;
    (*m)(2) = 0;
    return m;
}

DenseMatrixHandle CreateDenseMatrix_2()  
{
    DenseMatrixHandle m(boost::make_shared<DenseMatrix>(3,3));

    (*m)(0,0) = 1;
    (*m)(0,1) = 0;
    (*m)(0,2) = 0;
    (*m)(1,0) = 0;
    (*m)(1,1) = 2;
    (*m)(1,2) = 0;
    (*m)(2,0) = 0;
    (*m)(2,1) = 0;
    (*m)(2,2) = 3;

    return m;
}

SparseRowMatrixHandle CreateSparseMatrix()  
{ 
    SparseRowMatrixHandle m(boost::make_shared<SparseRowMatrix>(3,3));
    m->insert(0,0) = 1;
    m->insert(1,1) = 2;
    m->insert(2,2) = 3;
    m->makeCompressed();
    return m;
}

DenseMatrixHandle CreateDenseMatrix()  
{
    DenseMatrixHandle m(boost::make_shared<DenseMatrix>(3,1));

    (*m)(0,0) = 1;
    (*m)(1,0) = 2;
    (*m)(2,0) = 3;

    return m;
}
 
TEST(MatrixTypeConversionTests, denseTocolumn)
{
   DenseMatrixHandle from(CreateDenseMatrix());
    
   DenseColumnMatrixHandle result =  matrix_convert::to_column(from);
   
   EXPECT_EQ(from->ncols(), result->ncols());
   EXPECT_EQ(from->nrows(), result->nrows());
   
   for (int i = 0; i < result->nrows(); i++)
    for (int j = 0; j < result->ncols(); j++)
        EXPECT_EQ((*result)(i, j),(*from)(i, j));
	   
}

TEST(MatrixTypeConversionTests, denseTocolumn2)
{
   DenseMatrixHandle from(CreateDenseMatrix_2());
    
   DenseColumnMatrixHandle result =  matrix_convert::to_column(from);
   
   DenseColumnMatrixHandle expected_result(CreateColumnMatrix_2());
   
   EXPECT_EQ(expected_result->ncols(), result->ncols());
   EXPECT_EQ(expected_result->nrows(), result->nrows());
  
   for (int i = 0; i < result->nrows(); i++)
    for (int j = 0; j < result->ncols(); j++)
        EXPECT_EQ((*result)(i, j),(*expected_result)(i, j));
    
}

TEST(MatrixTypeConversionTests, denseTosparse)
{
   DenseMatrixHandle from(CreateDenseMatrix_2());
    
   SparseRowMatrixHandle result =  matrix_convert::to_sparse(from);
   
   if (!result)
   {
    std::cout << " Error: resulting null pointer !" << std::endl; 
   }
   
   SparseRowMatrixHandle expected_result  = CreateSparseMatrix();
   
   EXPECT_EQ(expected_result->ncols(), result->ncols());
   EXPECT_EQ(expected_result->nrows(), result->nrows());      
   
  for (int i = 0; i < result->nrows(); i++)
    for (int j = 0; j < result->ncols(); j++)
        EXPECT_EQ(expected_result->coeff(i, j),result->coeff(i, j));
	
}

TEST(MatrixTypeConversionTests, colvectorTosparse)
{  
  DenseColumnMatrixHandle from(CreateColumnMatrix());

  SparseRowMatrixHandle result =  matrix_convert::to_sparse(from);
  
  if (!result)
   {
    std::cout << " Error: resulting null pointer !" << std::endl; 
   }
  
  SparseRowMatrixHandle expected_result  = CreateSparseMatrixWithOneColumn();
  
   if (!result)
   {
    std::cout << " Error: resulting null pointer !" << std::endl; 
   }
   
   EXPECT_EQ(expected_result->ncols(), result->ncols());
   EXPECT_EQ(expected_result->nrows(), result->nrows());      
   
  for (int i = 0; i < result->nrows(); i++)
    for (int j = 0; j < result->ncols(); j++)
        EXPECT_EQ(expected_result->coeff(i, j),result->coeff(i, j));
}

TEST(MatrixTypeConversionTests, sparseTocol)
{  
  SparseRowMatrixHandle from = CreateSparseMatrix(); 
  
  DenseColumnMatrixHandle result =  matrix_convert::to_column(from);
  
  if (!result)
  {
    std::cout << " Error: resulting null pointer !" << std::endl; 
  } 
   
  DenseColumnMatrixHandle  expected_result  = CreateColumnMatrix_2(); 
   
  EXPECT_EQ(expected_result->ncols(), result->ncols());
  EXPECT_EQ(expected_result->nrows(), result->nrows());   
  
  for (int i = 0; i < result->nrows(); i++)
    for (int j = 0; j < result->ncols(); j++)
       EXPECT_EQ(expected_result->coeff(i, j),result->coeff(i, j)); 
	
}

TEST(MatrixTypeConversionTests, sparseTodense)
{  
  SparseRowMatrixHandle from = CreateSparseMatrix(); 
  
  DenseMatrixHandle result =  matrix_convert::to_dense(from);
  
  if (!result)
  {
    std::cout << " Error: resulting null pointer !" << std::endl; 
  } 
   
  DenseMatrixHandle  expected_result  = CreateDenseMatrix_2(); 
   
  EXPECT_EQ(expected_result->ncols(), result->ncols());
  EXPECT_EQ(expected_result->nrows(), result->nrows());   
  
  for (int i = 0; i < result->nrows(); i++)
    for (int j = 0; j < result->ncols(); j++)
        EXPECT_EQ(expected_result->coeff(i, j),result->coeff(i, j));
	
}


TEST(MatrixTypeConversionTests, colTodense)
{
 DenseColumnMatrixHandle from(CreateColumnMatrix()); 
 
 DenseMatrixHandle result =  matrix_convert::to_dense(from);

 if (!result)
 {
    std::cout << " Error: resulting null pointer !" << std::endl; 
 } 
 
 DenseMatrixHandle  expected_result  = CreateDenseMatrix();
 
 EXPECT_EQ(expected_result->ncols(), result->ncols());
 EXPECT_EQ(expected_result->nrows(), result->nrows());
 
  for (int i = 0; i < result->nrows(); i++)
    for (int j = 0; j < result->ncols(); j++)
        EXPECT_EQ(expected_result->coeff(i, j),result->coeff(i, j)); 
 
}
