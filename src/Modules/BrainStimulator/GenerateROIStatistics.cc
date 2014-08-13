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
*/
#include <iostream>
#include <Core/Datatypes/String.h>
#include <Core/Datatypes/Scalar.h>
#include <Modules/BrainStimulator/GenerateROIStatistics.h>
#include <Core/Algorithms/BrainStimulator/GenerateROIStatisticsAlgorithm.h>
#include <Core/Datatypes/Legacy/Field/Field.h>
#include <Core/Datatypes/DenseMatrix.h>

//////////////////////////////////////////////////////////////////////////
/// @todo MORITZ
//////////////////////////////////////////////////////////////////////////
using namespace SCIRun::Modules::BrainStimulator;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Core::Algorithms::BrainStimulator;
using namespace SCIRun::Dataflow::Networks;

GenerateROIStatisticsModule::GenerateROIStatisticsModule() : Module(ModuleLookupInfo("GenerateROIStatistics", "BrainStimulator", "SCIRun"))
{
 INITIALIZE_PORT(MESH_DATA_ON_ELEMENTS);
 INITIALIZE_PORT(PHYSICAL_UNIT);
 INITIALIZE_PORT(ATLAS_MESH);
 INITIALIZE_PORT(ATLAS_MESH_LABELS);
 INITIALIZE_PORT(COORDINATE_SPACE);
 INITIALIZE_PORT(STATISTICAL_RESULTS);
 INITIALIZE_PORT(COORDINATE_SPACE_LABEL);
}

void GenerateROIStatisticsModule::setStateDefaults()
{
  /// @todo
}

void GenerateROIStatisticsModule::execute()
{
  auto mesh_data = getRequiredInput(MESH_DATA_ON_ELEMENTS);
  auto physical_unit = getOptionalInput(PHYSICAL_UNIT);
  auto atlas_mesh = getRequiredInput(ATLAS_MESH);
  auto atlas_mesh_labels = getOptionalInput(ATLAS_MESH_LABELS);
  auto coordinate_space = getOptionalInput(COORDINATE_SPACE);
  auto coordinate_space_label = getOptionalInput(COORDINATE_SPACE);
  
  //auto elc_vals_from_state = get_state()->getValue(Parameters::ElectrodeTableValues).getList();
  //algo().set(Parameters::ELECTRODE_VALUES, elc_vals_from_state);
  
  //algorithm input and run
  auto output = algo().run_generic(make_input((MESH_DATA_ON_ELEMENTS, mesh_data)(PHYSICAL_UNIT, optionalAlgoInput(physical_unit))(ATLAS_MESH, atlas_mesh)(ATLAS_MESH_LABELS, optionalAlgoInput(atlas_mesh_labels))(COORDINATE_SPACE, optionalAlgoInput(coordinate_space))(COORDINATE_SPACE_LABEL, optionalAlgoInput(coordinate_space_label))));

  //algorithm output
  sendOutputFromAlgorithm(STATISTICAL_RESULTS, output);
}
