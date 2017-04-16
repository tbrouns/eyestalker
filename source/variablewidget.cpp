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

#include "variablewidget.h"

VariableWidget::VariableWidget(QWidget *parent) : QWidget(parent)
{
    // Circumference

    QLabel *CircumferenceTextBox = new QLabel;
    CircumferenceTextBox->setText("<b> Circumference:</b>");

    CircumferenceLabel  = new QLabel();
    CircumferenceSlider = new SliderDouble();
    CircumferenceSlider->setPrecision(1);
    CircumferenceSlider->setDoubleRange(0, 500);
    CircumferenceSlider->setOrientation(Qt::Horizontal);

    // Aspect ratio

    QLabel *AspectRatioTextBox = new QLabel;
    AspectRatioTextBox->setText("<b> Aspect ratio:</b>");

    AspectRatioLabel  = new QLabel();
    AspectRatioSlider = new SliderDouble();
    AspectRatioSlider->setPrecision(2);
    AspectRatioSlider->setDoubleRange(0, 1.0);
    AspectRatioSlider->setOrientation(Qt::Horizontal);

    // Title

    QLabel *TitleWidget = new QLabel;
    TitleWidget->setText("<b>Real-time variables</b>");

    // Set-up layout

    QGridLayout *MainLayout = new QGridLayout;

    MainLayout->addWidget(TitleWidget,          0, 1, 1, 1, Qt::AlignCenter);
    MainLayout->addWidget(CircumferenceTextBox, 1, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CircumferenceSlider,  1, 1);
    MainLayout->addWidget(CircumferenceLabel,   1, 2);
    MainLayout->addWidget(AspectRatioTextBox,   2, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AspectRatioSlider,    2, 1);
    MainLayout->addWidget(AspectRatioLabel,     2 ,2);

    MainLayout->setColumnStretch(0,1);
    MainLayout->setColumnStretch(1,3);
    MainLayout->setColumnStretch(2,1);

    setLayout(MainLayout);
}

VariableWidget::~VariableWidget()
{

}

void VariableWidget::setWidgets(const dataVariables& mDataVariables)
{
    CircumferenceSlider->setDoubleValue(mDataVariables.exactCircumference);
    CircumferenceLabel->setText(QString::number(mDataVariables.exactCircumference, 'f', 1));

    AspectRatioSlider->setDoubleValue(mDataVariables.exactAspectRatio);
    AspectRatioLabel->setText(QString::number(mDataVariables.exactAspectRatio, 'f', 2));
}
