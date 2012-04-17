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

#ifndef INTERFACE_APPLICATION_MODULEWIDGET_H
#define INTERFACE_APPLICATION_MODULEWIDGET_H

#include "ui_Module.h"
#include <boost/shared_ptr.hpp>
#include <QFrame>
#include <set>
#include <Interface/Application/PositionProvider.h>

#include <Core/Dataflow/Network/NetworkFwd.h>

class QGraphicsProxyWidget;

namespace SCIRun {
namespace Gui {

class PortWidget;
class InputPortWidget;
class OutputPortWidget;
class PositionProvider;

class ModuleWidget : public QFrame, public NeedsScenePositionProvider, public Ui::Module
{
	Q_OBJECT
	
public:
  explicit ModuleWidget(const QString& name, const SCIRun::Domain::Networks::ModuleInfoProvider& moduleInfoProvider, QWidget* parent = 0);
  ~ModuleWidget();

  void trackConnections();
  QPointF inputPortPosition() const;
  QPointF outputPortPosition() const;

  double percentComplete() const;
  void setPercentComplete(double p);

  size_t numInputPorts() const;
  size_t numOutputPorts() const;
  //TODO abstract
  const std::vector<PortWidget*>& getInputPorts() const { return inputPorts_; }
  const std::vector<PortWidget*>& getOutputPorts() const { return outputPorts_; }
  
public slots:
  //for testing signal/slot of Execute
  void incrementProgressFake();
signals:
  void removeModule(const std::string& moduleId);
  void addConnection(const std::string& id1, size_t port1, const std::string& id2, size_t port2);
private:
  std::vector<PortWidget*> inputPorts_;
  std::vector<PortWidget*> outputPorts_;

  void addPorts(const SCIRun::Domain::Networks::ModuleInfoProvider& moduleInfoProvider);
  void addPort(InputPortWidget* port);
  void addPort(OutputPortWidget* port);
  std::string moduleId_;
  //
  void addPortLayouts();
  QHBoxLayout* outputPortLayout_;
  QHBoxLayout* inputPortLayout_;
};

}
}

#endif
