/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2015 Scientific Computing and Imaging Institute,
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

#ifndef MATRIX_TEST_CASES_H
#define MATRIX_TEST_CASES_H

#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/MatrixTypeConversions.h>

namespace SCIRun 
{
namespace TestUtils
{
  static inline Core::Datatypes::DenseMatrix matrix1()
  {
    Core::Datatypes::DenseMatrix m(3, 3);
    for (int i = 0; i < m.rows(); ++ i)
      for (int j = 0; j < m.cols(); ++ j)
        m(i, j) = 3.0 * i + j;
    return m;
  }
  static inline Core::Datatypes::SparseRowMatrixHandle matrix1sparse()
  {
    return Core::Datatypes::matrix_convert::to_sparse(boost::make_shared<Core::Datatypes::DenseMatrix>(matrix1()));
  }
  static inline Core::Datatypes::DenseColumnMatrixHandle matrix1column()
  {
    return Core::Datatypes::matrix_convert::to_column(boost::make_shared<Core::Datatypes::DenseMatrix>(matrix1()));
  }
  static inline Core::Datatypes::DenseMatrix matrixNonSquare()
  {
    Core::Datatypes::DenseMatrix m(3, 4);
    for (int i = 0; i < m.rows(); ++ i)
      for (int j = 0; j < m.cols(); ++ j)
        m(i, j) = 3.0 * i + j;
    return m;
  }
  const Core::Datatypes::DenseMatrix Zero(Core::Datatypes::DenseMatrix::Zero(3, 3));
}
}

#define EXPECT_SPARSE_EQ(m1, m2) EXPECT_EQ(\
  *SCIRun::Core::Datatypes::matrix_convert::to_dense(boost::make_shared<SCIRun::Core::Datatypes::SparseRowMatrix>(m1)), \
  *SCIRun::Core::Datatypes::matrix_convert::to_dense(boost::make_shared<SCIRun::Core::Datatypes::SparseRowMatrix>(m2)))

#endif