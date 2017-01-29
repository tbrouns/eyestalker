#ifndef VARIABLEWIDGET_H
#define VARIABLEWIDGET_H

#include <QGridLayout>
#include <QLabel>
#include <QSlider>
#include <QObject>
#include <QWidget>

#include "headers/structures.h"
#include "headers/sliderdouble.h"
#include "headers/parameters.h"

class VariableWidget : public QWidget
{
    Q_OBJECT

public:

    void setWidgets(const detectionVariables&);
    detectionVariables getStructure();
    void resetStructure(const detectionParameters&);

    explicit VariableWidget(QWidget *parent = 0);
    ~VariableWidget();

private:

    SliderDouble *CircumferenceSlider;
    SliderDouble *AspectRatioSlider;
    QLabel *CircumferenceLabel;
    QLabel *AspectRatioLabel;

    detectionVariables mDetectionVariables;

signals:

public slots:

private slots:

    void setCircumference(double);
    void setAspectRatio(double);

};

#endif // VARIABLEWIDGET_H
