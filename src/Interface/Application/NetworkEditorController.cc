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

#include <iostream>

#include <Interface/Application/NetworkEditorController.h>

#include <Core/Dataflow/Network/Network.h>
#include <Core/Dataflow/Network/HardCodedModuleFactory.h>
#include <Core/Dataflow/Network/ModuleDescription.h>
#include <Core/Dataflow/Network/Module.h>

using namespace SCIRun;
using namespace SCIRun::Gui;
using namespace SCIRun::Domain::Networks;

NetworkEditorController::NetworkEditorController()
{
  ModuleFactoryHandle mf(new HardCodedModuleFactory);
  theNetwork_.reset(new Network(mf));
}

void NetworkEditorController::addModule(const QString& moduleName)
{
  ModuleLookupInfo info;
  info.module_name_ = moduleName.toStdString();
  ModuleHandle realModule = theNetwork_->add_module(info);
  emit moduleAdded(moduleName, *realModule);
  printNetwork();
}

void NetworkEditorController::removeModule(const std::string& id)
{
  theNetwork_->remove_module(id);
  printNetwork();
}

void NetworkEditorController::printNetwork() const
{
  if (theNetwork_)
    std::cout << theNetwork_->toString() << std::endl;
}

