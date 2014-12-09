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

#ifndef INTERFACE_MODULES_FORWARD_BUILDBEMATRIXDIALOG_H
#define INTERFACE_MODULES_FORWARD_BUILDBEMATRIXDIALOG_H

#include "Interface/Modules/Forward/ui_BuildBEMatrix.h"
#include <Interface/Modules/Base/ModuleDialogGeneric.h>
#include <Dataflow/Network/ModuleStateInterface.h>
#include <Interface/Modules/Forward/share.h>

namespace SCIRun {
namespace Gui {
  
class SCISHARE BuildBEMatrixDialog : public ModuleDialogGeneric, 
  public Ui::BuildBEMatrix
{
	Q_OBJECT
	
public:
  BuildBEMatrixDialog(const std::string& name, 
    SCIRun::Dataflow::Networks::ModuleStateHandle state,
    QWidget* parent = 0);
public Q_SLOTS:
  virtual void updateFromPortChange(int numPorts) override;
private Q_SLOTS:
  void pull();
  void push();
private:
  QComboBox* makeComboBoxItem(int i) const;
};

}
}

#endif
