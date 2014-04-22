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

#include <Core/Algorithms/Math/SelectSubMatrix.h>
#include <Core/Algorithms/Base/AlgorithmVariableNames.h>
#include <Core/Algorithms/Base/AlgorithmPreconditions.h>
#include <Core/Datatypes/SparseRowMatrix.h>
#include <Core/Datatypes/SparseRowMatrixFromMap.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixTypeConversions.h>

using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Core::Algorithms::Math;
using namespace SCIRun::Core::Datatypes;

SelectSubMatrixAlgorithm::SelectSubMatrixAlgorithm()
{
  addParameter(rowCheckBox, false);
  addParameter(columnCheckBox, false);
  addParameter(rowStartSpinBox, 0);
  addParameter(columnStartSpinBox, 0);
  addParameter(columnEndSpinBox, 0);
  addParameter(rowEndSpinBox, 0);
}


MatrixHandle SelectSubMatrixAlgorithm::get_sub_matrix(MatrixHandle& input_matrix, DenseMatrixHandle rows, DenseMatrixHandle cols) const
{

  std::vector<index_type> sel_rows;
  std::vector<index_type> sel_cols;

  try
  {
    size_type count=0;
    if (rows)
     if (rows->nrows()>0 && rows->ncols()>0)
     {
      index_type nr=rows->nrows(), nc=rows->ncols();
      if (nr==0) nr=1; if (nc==0) nc=1;
      sel_rows.resize(nr*nc);

      for (index_type i=0; i < nr; i++) 
        for (index_type j=0; j < nc; j++) 
	  {
           sel_rows[count++] = static_cast<index_type>((* rows)(i,j)); 
	  }
  
     }
    
    count=0;
    
    if (cols)
     if (cols->nrows()>0 && cols->ncols()>0)
     {
      index_type nr=cols->nrows(), nc=cols->ncols();
      if (nr==0) nr=1; if (nc==0) nc=1;
      sel_cols.resize(nr*nc);
      for (index_type i=0; i < nr; i++) 
        for (index_type j=0; j < nc; j++) 
           sel_cols[count++] = static_cast<index_type>((* cols)(i,j));	   
     } 
  }
  catch (...)
  {	  
    THROW_ALGORITHM_INPUT_ERROR(" Could not allocate enough memory ");
  }
 
 return run(input_matrix,sel_rows,sel_cols);
}

MatrixHandle SelectSubMatrixAlgorithm::run(MatrixHandle input_matrix, DenseMatrixHandle row_indices, DenseMatrixHandle col_indices) const
{  

  if (!input_matrix || input_matrix->nrows()==0 || input_matrix->ncols()==0)
  {
    remark(" No valid inputs: input matrix or row,column matrix contain NULL pointer ");
    return MatrixHandle();
  }
   
  MatrixHandle sub_matrix;
  bool row_select = get(rowCheckBox).getBool();
  bool col_select = get(columnCheckBox).getBool();
  index_type row_start = get(rowStartSpinBox).getInt();
  index_type row_end = get(rowEndSpinBox).getInt();
  index_type col_start = get(columnStartSpinBox).getInt();
  index_type col_end = get(columnEndSpinBox).getInt();
  
  if ( !row_select && !col_select )  //pipe input through
   if ( !row_indices && !col_indices)
    { 
      return input_matrix;
    }
    
  if (row_indices || col_indices)
    {
      if (row_select || col_select)
      {
        remark("Index matrices detected on inputs (indexing starts from 0), ignoring UI settings");
      }
    }
 
 if  (! (row_indices || col_indices))   
 {
   if (row_select && col_select)
   {
     row_end=row_end-row_start+1;
     col_end=col_end-col_start+1;
   } else
   if (row_select && !col_select)
    {
      col_start=0;
      row_end=row_end-row_start+1;
      col_end=input_matrix->ncols();
    } else
    if ( !row_select && col_select )
    {
      row_start=0;
      row_end=input_matrix->nrows();
      col_end=col_end-col_start+1;
    } else
    {  //pass it through
     row_start=0;
     col_start=0;
     row_end=input_matrix->nrows();
     col_end=input_matrix->ncols();
    }
 }   
  
  if ( row_start<0 || col_start<0 ||  row_end>input_matrix->nrows() || col_end>input_matrix->ncols()) 
    {
     remark(" Specified matrix indices from UI settings exceed matrix dimensions ");
     return MatrixHandle();
    }
    
  if (row_indices || col_indices) 
  {   
    sub_matrix=get_sub_matrix(input_matrix,row_indices,col_indices);
    if (!sub_matrix) return MatrixHandle(); 
  } 
  else  //GUI input only
  {     
    if(matrix_is::sparse(input_matrix))
     {
       SparseRowMatrixHandle mat (new  SparseRowMatrix(matrix_cast::as_sparse(input_matrix)->block(row_start,col_start,row_end,col_end)));
       return mat;
     } else
     if(matrix_is::dense(input_matrix))
      {
        DenseMatrixHandle mat (new  DenseMatrix(matrix_cast::as_dense(input_matrix)->block(row_start,col_start,row_end,col_end)));
	return mat;
      }  else
      {
        remark("This module needs row indices, or column indices or both from UI or input matrices");
        remark("Copying input matrix to output");
        sub_matrix=input_matrix;       
     }
  } 
    
 return sub_matrix;
}

MatrixHandle SelectSubMatrixAlgorithm::run(MatrixHandle input_matrix, std::vector<SCIRun::index_type> rows, std::vector<SCIRun::index_type> cols) const
{

  if (!input_matrix)
  {
    THROW_ALGORITHM_INPUT_ERROR("No input matrix");    
  }
  
  if ((rows.size() == 0) && (cols.size() == 0))
  {
    THROW_ALGORITHM_INPUT_ERROR("No row and column indices given");
  }
  
  size_type m = input_matrix->nrows();
  size_type n = input_matrix->ncols();

  for (size_t r=0; r<rows.size(); r++)
  {
    if (rows[r] >= static_cast<index_type>(m) || rows[r] < 0 || !IsFinite(rows[r]))
    {
      THROW_ALGORITHM_INPUT_ERROR("Selected row exceeds matrix dimensions");
    }
  }

  for (size_t r=0; r<cols.size(); r++)
  {
    if (cols[r] >= static_cast<index_type>(n) || cols[r] < 0 || !IsFinite(cols[r]))
    {
      THROW_ALGORITHM_INPUT_ERROR("Selected column exceeds matrix dimensions");
    }
  }

  auto sparse_matrix = matrix_cast::as_sparse(input_matrix);

  if (sparse_matrix)
  {

   SparseRowMatrixFromMap::Values additionalData;
   if (rows.size()>0 && cols.size()>0) //get only the indices intersection
   {
     m=rows.size();
     n=cols.size();
     for (index_type i=0; i< rows.size(); i++) 
      for (index_type j=0; j < cols.size(); j++) 
      {
       auto tmp = sparse_matrix->coeff(rows[i],cols[j]);
       if (tmp) additionalData[rows[i]][cols[j]]=tmp;
	 
      }
   }

   if (rows.size()>0 && cols.size()==0)
   {
     m=rows.size();
     for (size_t i=0; i<rows.size(); i++)
     {
      auto  matrix_row=sparse_matrix->row(rows[i]);
      for (Eigen::SparseVector<double>::InnerIterator it(matrix_row); it; ++it)
      {
       additionalData[rows[i]][it.index()]=it.value(); 
      }
     }
   }
   
   if (rows.size()==0 && cols.size()>0)
   {
     n=cols.size();
     for (size_t j=0; j<cols.size(); j++)
     {
      auto  matrix_col=sparse_matrix->col(cols[j]);
      for (Eigen::SparseVector<double>::InnerIterator it(matrix_col); it; ++it)
      {
       additionalData[it.index()][cols[j]]=it.value();
      }
     }
   }
   
   auto output_matrix = SparseRowMatrixFromMap::make(m, n, additionalData);

  } else
  {
   auto dense_input_matrix = matrix_cast::as_dense(input_matrix);

   if (!dense_input_matrix)
   {
     THROW_ALGORITHM_INPUT_ERROR("Could not convert matrix into dense matrix");
   }
   
   if (rows.size()>0 && cols.size()>0)
   {
      DenseMatrix mat(rows.size(),cols.size());

      for (index_type i = 0; i < rows.size(); i++)
        for (index_type j = 0; j < cols.size(); j++)
         (mat)(i,j) = dense_input_matrix->coeff(rows[i],cols[j]);
      
      DenseMatrixHandle  output(boost::make_shared<DenseMatrix>(mat));   
   
      return output;
   }
   
   
   if (rows.size()>0 && cols.size()==0)
   {
      DenseMatrix mat(rows.size(),n);

      for (index_type i = 0; i < rows.size(); i++)
        mat.row(i) = dense_input_matrix->row(rows[i]);
            
      DenseMatrixHandle  output(boost::make_shared<DenseMatrix>(mat));       
      
      return output;
   }
   
   if (rows.size()==0 && cols.size()>0)
   {
      DenseMatrix mat(m,cols.size());

      for (index_type i = 0; i < cols.size(); i++)
         if (cols[i]>0) mat.col(i) = dense_input_matrix->col(cols[i]);
      
      DenseMatrixHandle  output(boost::make_shared<DenseMatrix>(mat));   
   
      return output;
   }
      
  }   
  
  remark("Could not deal with input data."); 
  remark("Copying input matrix to output"); 
  return input_matrix;
}

AlgorithmInputName SelectSubMatrixAlgorithm::InputMatrix("InputMatrix");
AlgorithmOutputName SelectSubMatrixAlgorithm::RowIndicies("RowIndicies");
AlgorithmOutputName SelectSubMatrixAlgorithm::ColumnIndicies("ColumnIndicies");
AlgorithmParameterName SelectSubMatrixAlgorithm::rowCheckBox("rowCheckBox");
AlgorithmParameterName SelectSubMatrixAlgorithm::columnCheckBox("columnCheckBox");
AlgorithmParameterName SelectSubMatrixAlgorithm::rowStartSpinBox("rowStartSpinBox");
AlgorithmParameterName SelectSubMatrixAlgorithm::columnStartSpinBox("columnStartSpinBox");
AlgorithmParameterName SelectSubMatrixAlgorithm::columnEndSpinBox("columnEndSpinBox");
AlgorithmParameterName SelectSubMatrixAlgorithm::rowEndSpinBox("rowEndSpinBox");
    
AlgorithmOutput SelectSubMatrixAlgorithm::run_generic(const AlgorithmInput& input) const
{
  auto input_matrix = input.get<Matrix>(InputMatrix);
  auto rowindicies = input.get<DenseMatrix>(RowIndicies);
  auto columnindicies = input.get<DenseMatrix>(ColumnIndicies);

  MatrixHandle output_matrix;
  output_matrix = run(input_matrix, rowindicies, columnindicies);
  
  AlgorithmOutput output;
  output[Variables::ResultMatrix] = output_matrix;
  return output;
}
