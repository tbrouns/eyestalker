#include "qwtplotwidget.h"

QwtPlotWidget::QwtPlotWidget(QWidget *parent) : QWidget(parent)
{
    mQwtPlot = new QwtPlot;
    mQwtPlot->setAxisScale(QwtPlot::xBottom, 0, Parameters::eyeAOI.wdth, 0);
    mQwtPlot->setAxisScale(QwtPlot::yLeft,   0, Parameters::eyeAOI.hght, 0);
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
    double aspectRatio = (double) Parameters::eyeAOI.hght / Parameters::eyeAOI.wdth;
    double height = aspectRatio * widgetWdth;
    return QSize(widgetWdth, height);
}

void QwtPlotWidget::plotTrajectory(const std::vector<double>& x, const std::vector<double>& y)
{
    mQwtPlot->detachItems();

    QwtPlotCurve *curve = new QwtPlotCurve();
    curve->setSamples( x.data(), y.data(), (int) x.size() );
    curve->setPen(* new QPen(Qt::red));
    curve->attach(mQwtPlot);

    mQwtPlot->setAxisScale(QwtPlot::xBottom, 0, Parameters::eyeAOI.wdth, 0);
    mQwtPlot->setAxisScale(QwtPlot::yLeft,   0, Parameters::eyeAOI.hght, 0);

    mQwtPlot->replot();
}

void QwtPlotWidget::plotTimeSeries(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& t, double trialLength)
{
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
    mQwtPlot->setAxisScale(QwtPlot::yLeft,   0, Parameters::eyeAOI.wdth, 0);

    mQwtPlot->replot();
}


