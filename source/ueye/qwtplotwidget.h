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

#ifndef QWTPLOTWIDGET_H
#define QWTPLOTWIDGET_H

#include "../parameters.h"

// Qt

#include <QObject>
#include <QVBoxLayout>
#include <QWidget>

// Qwt

#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_scaleitem.h>

class QwtPlotWidget : public QWidget
{
    Q_OBJECT

public:
    explicit QwtPlotWidget(QWidget *parent = 0);

    QSize sizeHint() const;

    void plotTimeSeries(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, double);
    void plotTrajectory(const std::vector<double>&, const std::vector<double>&);

    void plotData(std::vector<double>, std::vector<double>);
    void setWidth(int);

private:

    QwtPlot* mQwtPlot;
    int widgetWdth;

signals:

public slots:
};

#endif // QWTPLOTWIDGET_H
