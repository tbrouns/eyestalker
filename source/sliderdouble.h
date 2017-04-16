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

#ifndef SLIDERDOUBLE_H
#define SLIDERDOUBLE_H

// Libraries

// Qt

#include <QSlider>
#include <QWidget>

class SliderDouble : public QSlider
{
    Q_OBJECT

public:

    SliderDouble(QWidget *parent = 0);

    void setDoubleRange(double min, double max);
    void setDoubleMaximum(double max);
    void setDoubleMinimum(double min);
    void setPrecision(int val);

private:

    double precision;

signals:

    void doubleValueChanged(double value);

public slots:

    void setDoubleValue(double val);

    void notifyValueChanged(int value)
    {
        double doubleValue = value / precision;
        emit doubleValueChanged(doubleValue);
    }
};

#endif // SLIDERDOUBLE_H
