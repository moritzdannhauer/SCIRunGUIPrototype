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

#include <Interface/Modules/Visualization/ShowFieldDialog.h>
#include <Modules/Visualization/ShowField.h>
#include <Dataflow/Network/ModuleStateInterface.h>  //TODO: extract into intermediate
#include <Core/Datatypes/Color.h>
#include <QColorDialog>

using namespace SCIRun::Gui;
using namespace SCIRun::Dataflow::Networks;
using namespace SCIRun::Modules::Visualization;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Algorithms::Visualization;

ShowFieldDialog::ShowFieldDialog(const std::string& name, ModuleStateHandle state,
  QWidget* parent /* = 0 */)
  : ModuleDialogGeneric(state, parent),
  defaultMeshColor_(Qt::gray)
{
  setupUi(this);
  setWindowTitle(QString::fromStdString(name));
  fixSize();
  WidgetStyleMixin::tabStyle(this->displayOptionsTabs_);
  addCheckBoxManager(showNodesCheckBox_, ShowFieldModule::ShowNodes);
  addCheckBoxManager(showEdgesCheckBox_, ShowFieldModule::ShowEdges);
  addCheckBoxManager(showFacesCheckBox_, ShowFieldModule::ShowFaces);
  addCheckBoxManager(enableTransparencyNodesCheckBox_, ShowFieldModule::NodeTransparency);
  addCheckBoxManager(enableTransparencyEdgesCheckBox_, ShowFieldModule::EdgeTransparency);
  addCheckBoxManager(enableTransparencyFacesCheckBox_, ShowFieldModule::FaceTransparency);
  addCheckBoxManager(invertNormalsCheckBox, ShowFieldModule::FaceInvertNormals);
  addDoubleSpinBoxManager(transparencyDoubleSpinBox_, ShowFieldModule::FaceTransparencyValue);
  addDoubleSpinBoxManager(nodeTransparencyDoubleSpinBox_, ShowFieldModule::NodeTransparencyValue);
  addDoubleSpinBoxManager(edgeTransparencyDoubleSpinBox_, ShowFieldModule::EdgeTransparencyValue);
  addDoubleSpinBoxManager(scaleSphereDoubleSpinBox_, ShowFieldModule::SphereScaleValue);
  addDoubleSpinBoxManager(cylinder_rad_spin, Parameters::CylinderRadius);
  addSpinBoxManager(cylinder_res_spin, ShowFieldModule::CylinderResolution);
  addSpinBoxManager(sphereResolutionSpinBox, ShowFieldModule::SphereResolution);
  addRadioButtonGroupManager({ edgesAsLinesButton_, edgesAsCylindersButton_ }, ShowFieldModule::EdgesAsCylinders);
  addRadioButtonGroupManager({ nodesAsPointsButton_, nodesAsSpheresButton_ }, ShowFieldModule::NodeAsSpheres);

  connectButtonsToExecuteSignal({ showNodesCheckBox_, showEdgesCheckBox_, showFacesCheckBox_, enableTransparencyNodesCheckBox_,
    enableTransparencyEdgesCheckBox_, enableTransparencyFacesCheckBox_, invertNormalsCheckBox, edgesAsLinesButton_,
    edgesAsCylindersButton_, nodesAsPointsButton_, nodesAsSpheresButton_ });

  createExecuteInteractivelyToggleAction();

  connect(defaultMeshColorButton_, SIGNAL(clicked()), this, SLOT(assignDefaultMeshColor()));

  /////Set unused widgets to be not visible
  //Nodes Tab
  //label_4->setVisible(false); // Sphere scale lable
  //scaleSphereDoubleSpinBox_->setVisible(false); // Sphere scale spin box
  //resolutionSpinBox->setVisible(false); //resolution spin box
  //label_5->setVisible(false); //resolution label
  //groupBox_3->setVisible(false); //Node coloring
  //groupBox_4->setVisible(false); //Node Display Type Group Box

  //Edges Tab
  //groupBox_7->setVisible(false);//Edge Display Type Group Box
  //label_9->setVisible(false); //resolution label
  //cylinder_res_spin->setVisible(false); //resolution spinbox
  //label_8->setVisible(false); //scale label
  //cylinder_rad_spin->setVisible(false); //cylinder scale spinbox
  //groupBox_6->setVisible(false); //edge coloring

  //Faces Tab
  //groupBox_5->setVisible(false); //face coloring
  //checkBox->setVisible(false); //Use Face Normal box
  //checkBox_2->setVisible(false); //Images as texture box
}

void ShowFieldDialog::push()
{
  if (!pulling_)
  {
    pushColor();
  }
}

void ShowFieldDialog::createStartupNote()
{
  auto showFieldId = windowTitle().split(':')[1];
  setStartupNote("ID: " + showFieldId);
}

void ShowFieldDialog::pullSpecial()
{
  ColorRGB color(state_->getValue(ShowFieldModule::DefaultMeshColor).toString());
  //std::cout << "pull color: " << color.r() << " " << color.g() << " " << color.b() << std::endl;
  // check for old saved color format: integers 0-255.
  defaultMeshColor_ = QColor(
    static_cast<int>(color.r() > 1 ? color.r() : color.r() * 255.0),
    static_cast<int>(color.g() > 1 ? color.g() : color.g() * 255.0),
    static_cast<int>(color.b() > 1 ? color.b() : color.b() * 255.0));
}

void ShowFieldDialog::assignDefaultMeshColor()
{
  auto newColor = QColorDialog::getColor(defaultMeshColor_, this, "Choose default mesh color");
  if (newColor.isValid())
  {
    defaultMeshColor_ = newColor;
    //TODO: set color of button to this color
    //defaultMeshColorButton_->set
    pushColor();
  }
}

void ShowFieldDialog::pushColor()
{
  //std::cout << "push color: " << defaultMeshColor_.redF() << " " << defaultMeshColor_.greenF() << " " << defaultMeshColor_.blueF() << std::endl;
  state_->setValue(ShowFieldModule::DefaultMeshColor, ColorRGB(defaultMeshColor_.redF(), defaultMeshColor_.greenF(), defaultMeshColor_.blueF()).toString());
  Q_EMIT executeActionTriggered();
}
