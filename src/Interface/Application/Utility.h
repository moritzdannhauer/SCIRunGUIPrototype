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

#ifndef UTILITY_H
#define UTILITY_H

#include <sstream>
#include <QAction>

namespace SCIRun {

template <class Point>
std::string to_string(const Point& p)
{
  std::ostringstream ostr;
  ostr << "QPoint(" << p.x() << "," << p.y() << ")";
  return ostr.str();
}

namespace Gui
{
  QColor to_color(const std::string& str, int alpha = 255);

  QColor defaultTagColor(int tag);
  typedef std::function<QColor(int)> TagColorFunc;

  QString colorToString(const QColor& color);

  inline QAction* separatorAction(QWidget* parent)
  {
    auto sep = new QAction(parent);
    sep->setSeparator(true);
    return sep;
  }

  inline QAction* disabled(QAction* action)
  {
    action->setEnabled(false);
    return action;
  }

  inline std::ostream& operator<<(std::ostream& o, const QPointF& p)
  {
    return o << "[" << p.x() << "," << p.y() << "]";
  }

  typedef boost::function<bool(const Dataflow::Networks::ModuleDescription&)> ModulePredicate;
  typedef boost::function<void(QAction*)> QActionHookup;
  QList<QAction*> fillMenuWithFilteredModuleActions(QMenu* menu, const Dataflow::Networks::ModuleDescriptionMap& moduleMap, ModulePredicate modulePred, QActionHookup hookup);
  QPointF findCenterOfNetwork(const Dataflow::Networks::ModulePositions& positions);
}

}

#endif
