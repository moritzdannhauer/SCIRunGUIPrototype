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

#include <Interface/Modules/Visualization/CreateStandardColorMapDialog.h>
#include <Core/Algorithms/Base/AlgorithmVariableNames.h>
#include <Modules/Visualization/CreateStandardColorMap.h>

using namespace SCIRun::Gui;
using namespace SCIRun::Dataflow::Networks;
using namespace SCIRun::Core::Algorithms::Visualization;
using namespace SCIRun::Core::Datatypes;

typedef SCIRun::Modules::Visualization::CreateStandardColorMap CreateStandardColorMapModule;


CreateStandardColorMapDialog::CreateStandardColorMapDialog(const std::string& name, ModuleStateHandle state,
  QWidget* parent /* = 0 */)
  : ModuleDialogGeneric(state, parent)
{
  setupUi(this);
  setWindowTitle(QString::fromStdString(name));
  ColorMap cm("Rainbow");
  previewColorMap_->setStyleSheet(buildGradientString(cm));

  addComboBoxManager(colorMapNameComboBox_, Parameters::ColorMapName);
  addSpinBoxManager(resolutionSpin_, Parameters::ColorMapResolution);
  addDoubleSpinBoxManager(shiftSpin_, Parameters::ColorMapShift);
  addCheckBoxManager(invertCheck_, Parameters::ColorMapInvert);

  connect(colorMapNameComboBox_, SIGNAL(currentIndexChanged(const QString&)), this, SLOT(updateColorMapPreview(const QString&)));
  connect(shiftSpin_, SIGNAL(valueChanged(double)), this, SLOT(setShiftSlider(double)));
  connect(resolutionSpin_, SIGNAL(valueChanged(int)), this, SLOT(setResolutionSlider(int)));
  connect(shiftSpin_, SIGNAL(valueChanged(double)), this, SLOT(updateColorMapPreview()));
  connect(resolutionSpin_, SIGNAL(valueChanged(int)), this, SLOT(updateColorMapPreview()));

  connect(shiftSlider_, SIGNAL(valueChanged(int)), this, SLOT(setShiftSpinner(int)));
  connect(resolutionSlider_, SIGNAL(valueChanged(int)), resolutionSpin_, SLOT(setValue(int)));
  connect(invertCheck_, SIGNAL(toggled(bool)), this, SLOT(onInvertCheck(bool)));
}

void CreateStandardColorMapDialog::updateColorMapPreview(const QString& s)
{
  ColorMap cm(s.toStdString(), resolutionSlider_->value(),
    static_cast<double>(shiftSlider_->value()) / 100.,
    invertCheck_->isChecked());
 // qDebug() << "updating color map: " << s << " " << resolutionSlider_->value() << " " << shiftSlider_->value();
  previewColorMap_->setStyleSheet(buildGradientString(cm));
}

void CreateStandardColorMapDialog::updateColorMapPreview()
{
  updateColorMapPreview(colorMapNameComboBox_->currentText());
}

const QString CreateStandardColorMapDialog::buildGradientString(const ColorMap& cm)
{
  //TODO: cache these values, GUI is slow to update.
  std::stringstream ss;
  ss << "background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0,";
  for (double i = 0.001; i < 1.0; i += 0.001) { //styling values need to be in the range [0,1]
    ss << " stop:" << i;
    ss << " rgba(";
    ColorRGB c = cm.valueToColor(i * 2. - 1.); //need to match default ColorMap data range [-1,1]
    ss << int(255.*c.r()) << ", " << int(255.*c.g()) << ", " << int(255.*c.b()) << ", 255),";
  }
  ss << ");";
  std::string str = ss.str();
  return QString::fromStdString(ss.str());
}

void CreateStandardColorMapDialog::setShiftSlider(double d)
{
  shiftSlider_->setValue(static_cast<int>(d * 100.));
}

void CreateStandardColorMapDialog::setResolutionSlider(int i)
{
  resolutionSlider_->setValue(i);
}

void CreateStandardColorMapDialog::setShiftSpinner(int i)
{
  shiftSpin_->setValue(static_cast<double>(i) / 100.);
}

void CreateStandardColorMapDialog::onInvertCheck(bool b)
{
  updateColorMapPreview();
}
