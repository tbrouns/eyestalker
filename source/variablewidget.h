//  EyeStalker: robust video-based eye tracking
//  Copyright (C) 2016  Terence Brouns, t.s.n.brouns@gmail.com

//  EyeStalker is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  EyeStalker is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>

#ifndef VARIABLEWIDGET_H
#define VARIABLEWIDGET_H

#include <QGridLayout>
#include <QLabel>
#include <QSlider>
#include <QObject>
#include <QWidget>

#include "structures.h"
#include "sliderdouble.h"
#include "parameters.h"

class VariableWidget : public QWidget
{
    Q_OBJECT

public:

    void setWidgets(const dataVariables&);

    explicit VariableWidget(QWidget *parent = 0);
    ~VariableWidget();

private:

    SliderDouble *CircumferenceSlider;
    SliderDouble *AspectRatioSlider;
    QLabel *CircumferenceLabel;
    QLabel *AspectRatioLabel;

};

#endif // VARIABLEWIDGET_H
