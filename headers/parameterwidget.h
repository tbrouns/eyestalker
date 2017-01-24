#ifndef PARAMETERWIDGET_H
#define PARAMETERWIDGET_H

#include <QObject>
#include <QWidget>

class ParameterWidget : public QWidget
{
    Q_OBJECT
public:
    explicit ParameterWidget(QWidget *parent = 0);
    ~ParameterWidget();

signals:

public slots:
};

#endif // PARAMETERWIDGET_H
