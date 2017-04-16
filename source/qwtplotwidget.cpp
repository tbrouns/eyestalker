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

#include "qwtplotwidget.h"

QwtPlotWidget::QwtPlotWidget(QWidget *parent) : QWidget(parent)
{
    std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);

    mQwtPlot = new QwtPlot;
    mQwtPlot->setAxisScale(QwtPlot::xBottom, 0, Parameters::camAOI.wdth, 0);
    mQwtPlot->setAxisScale(QwtPlot::yLeft,   0, Parameters::camAOI.hght, 0);
    mQwtPlot->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

    QVBoxLayout *MainLayout = new QVBoxLayout;
    MainLayout->addWidget(mQwtPlot);
    setLayout(MainLayout);
}

void QwtPlotWidget::setWidth(int W)
{
    widgetWdth = W;
}

QSize QwtPlotWidget::sizeHint() const
{
    std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);
    double aspectRatio = (double) Parameters::camAOI.hght / Parameters::camAOI.wdth;
    double height = aspectRatio * widgetWdth;
    return QSize(widgetWdth, height);
}

void QwtPlotWidget::plotTrajectory(const std::vector<double>& x, const std::vector<double>& y)
{
    std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);

    mQwtPlot->detachItems();

    QwtPlotCurve *curve = new QwtPlotCurve();
    curve->setSamples( x.data(), y.data(), (int) x.size() );
    curve->setPen(* new QPen(Qt::red));
    curve->attach(mQwtPlot);

    mQwtPlot->setAxisScale(QwtPlot::xBottom, 0, Parameters::camAOI.wdth, 0);
    mQwtPlot->setAxisScale(QwtPlot::yLeft,   0, Parameters::camAOI.hght, 0);

    mQwtPlot->replot();
}

void QwtPlotWidget::plotTimeSeries(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& t, double trialLength)
{
    std::lock_guard<std::mutex> AOICamLock(Parameters::AOICamMutex);

    mQwtPlot->detachItems();

    QwtPlotCurve *xCurve = new QwtPlotCurve();
    xCurve->setSamples( t.data(), x.data(), (int) x.size() );
    xCurve->setPen(* new QPen(Qt::blue));
    xCurve->attach(mQwtPlot);

    QwtPlotCurve *yCurve = new QwtPlotCurve();
    yCurve->setSamples( t.data(), y.data(), (int) x.size() );
    yCurve->setPen(* new QPen(Qt::red));
    yCurve->attach(mQwtPlot);

    mQwtPlot->setAxisScale(QwtPlot::xBottom, 0,             trialLength, 0);
    mQwtPlot->setAxisScale(QwtPlot::yLeft,   0, Parameters::camAOI.wdth, 0);

    mQwtPlot->replot();
}


