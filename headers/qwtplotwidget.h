#ifndef QWTPLOTWIDGET_H
#define QWTPLOTWIDGET_H

#include <headers/parameters.h>

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
