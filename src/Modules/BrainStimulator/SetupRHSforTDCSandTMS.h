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

#ifndef MODULES_BRAINSTIMULATOR_SetupRHSforTDCSandTMS_H
#define MODULES_BRAINSTIMULATOR_SetupRHSforTDCSandTMS_H

#include <Dataflow/Network/Module.h>
#include <Modules/BrainStimulator/share.h>

namespace SCIRun {
  namespace Modules {
    namespace BrainStimulator {

class SCISHARE SetupRHSforTDCSandTMSModule : public Dataflow::Networks::Module,
  public Has5InputPorts<FieldPortTag, MatrixPortTag, FieldPortTag, FieldPortTag, MatrixPortTag>,
  public Has5OutputPorts<MatrixPortTag, MatrixPortTag, MatrixPortTag, MatrixPortTag, MatrixPortTag>
{
  public:
    SetupRHSforTDCSandTMSModule();
    virtual void execute();
    virtual void setStateDefaults();

    INPUT_PORT(0, MESH, LegacyField);
    INPUT_PORT(1, ELECTRODE_COUNT, Matrix);
    INPUT_PORT(2, SCALP_TRI_SURF_MESH, LegacyField);
    INPUT_PORT(3, ELECTRODE_TRI_SURF_MESH, LegacyField);
    INPUT_PORT(4, ELECTRODE_SPONGE_LOCATION_AVR, Matrix);

    OUTPUT_PORT(0, ELECTRODE_ELEMENT, Matrix);
    OUTPUT_PORT(1, ELECTRODE_ELEMENT_TYPE, Matrix);
    OUTPUT_PORT(2, ELECTRODE_ELEMENT_DEFINITION, Matrix); 
    OUTPUT_PORT(3, ELECTRODE_CONTACT_IMPEDANCE, Matrix);
    OUTPUT_PORT(4, RHS, Matrix);
    
};

}}}

#endif
