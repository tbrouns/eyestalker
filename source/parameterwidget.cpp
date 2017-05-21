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

#include "parameterwidget.h"

ParameterWidget::ParameterWidget(QWidget *parent) : QWidget(parent)
{
    // Averages and Limits/Thresholds

    QLabel *ParametersTextBox = new QLabel;
    ParametersTextBox->setText("<b>Tracking parameters</b>");
    ParametersTextBox->setAlignment(Qt::AlignCenter);

    // Pupil circumference

    QLabel *CircumferenceMinTextBox = new QLabel;
    CircumferenceMinTextBox->setText("<b>Circumference min:</b>");

    CircumferenceMinLabel  = new QLabel;
    CircumferenceMinSlider = new SliderDouble;
    CircumferenceMinSlider->setPrecision(1);
    CircumferenceMinSlider->setDoubleRange(1, 500);
    CircumferenceMinSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CircumferenceMinSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCircumferenceMin(double)));

    QLabel *CircumferenceMaxTextBox = new QLabel;
    CircumferenceMaxTextBox->setText("<b>Circumference max:</b>");

    CircumferenceMaxLabel  = new QLabel;
    CircumferenceMaxSlider = new SliderDouble;
    CircumferenceMaxSlider->setPrecision(1);
    CircumferenceMaxSlider->setDoubleRange(1, 500);
    CircumferenceMaxSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CircumferenceMaxSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCircumferenceMax(double)));

    // Pupil aspect ratio

    QLabel *AspectRatioMinTextBox = new QLabel;
    AspectRatioMinTextBox->setText("<b>Aspect ratio min:</b>");

    AspectRatioMinLabel  = new QLabel;
    AspectRatioMinSlider = new SliderDouble;
    AspectRatioMinSlider->setPrecision(2);
    AspectRatioMinSlider->setDoubleRange(0.0, 1.0);
    AspectRatioMinSlider->setOrientation(Qt::Horizontal);
    QObject::connect(AspectRatioMinSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setAspectRatioMin(double)));

    // Sliders for canny edge parameters

    QLabel *CannyEdgeTextBox = new QLabel;
    CannyEdgeTextBox->setText("<b>Canny edge detection</b>");
    CannyEdgeTextBox->setAlignment(Qt::AlignCenter);

    QLabel *CannyThresholdLowTextBox = new QLabel;
    CannyThresholdLowTextBox->setText("<b>Low threshold:</b>");

    CannyThresholdLowLabel  = new QLabel;
    CannyThresholdLowSlider = new SliderDouble;
    CannyThresholdLowSlider->setPrecision(1);
    CannyThresholdLowSlider->setDoubleRange(0, 1000);
    CannyThresholdLowSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CannyThresholdLowSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCannyThresholdLow(double)));

    QLabel *CannyThresholdHighTextBox = new QLabel;
    CannyThresholdHighTextBox->setText("<b>High threshold:</b>");

    CannyThresholdHighLabel  = new QLabel;
    CannyThresholdHighSlider = new SliderDouble;
    CannyThresholdHighSlider->setPrecision(1);
    CannyThresholdHighSlider->setDoubleRange(0, 1000);
    CannyThresholdHighSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CannyThresholdHighSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCannyThresholdHigh(double)));

    QLabel *CannyBlurLevelTextBox = new QLabel;
    CannyBlurLevelTextBox->setText("<b>Blur level:</b>");

    CannyBlurLevelLabel  = new QLabel;
    CannyBlurLevelSlider = new QSlider;
    CannyBlurLevelSlider->setRange(0, 10);
    CannyBlurLevelSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CannyBlurLevelSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyBlurLevel(int)));

    QLabel *CannyKernelSizeTextBox = new QLabel;
    CannyKernelSizeTextBox->setText("<b>Kernel size:</b>");

    CannyKernelSizeLabel  = new QLabel;
    CannyKernelSizeSlider = new QSlider;
    CannyKernelSizeSlider->setRange(2, 4);
    CannyKernelSizeSlider->setOrientation(Qt::Horizontal);
    CannyKernelSizeSlider->setSingleStep(1);
    QObject::connect(CannyKernelSizeSlider, SIGNAL(valueChanged(int)), this, SLOT(setCannyKernelSize(int)));

    // Gain parameters

    QLabel *GainTextBox = new QLabel;
    GainTextBox->setText("<b>Gains</b>");
    GainTextBox->setAlignment(Qt::AlignCenter);

    QLabel *GainAveragesTextBox  = new QLabel;
    GainAveragesTextBox->setText("<b>Averages:</b>");

    GainAveragesLabel  = new QLabel;
    GainAveragesSlider = new SliderDouble;
    GainAveragesSlider->setPrecision(3);
    GainAveragesSlider->setDoubleRange(0, 0.1);
    GainAveragesSlider->setOrientation(Qt::Horizontal);
    QObject::connect(GainAveragesSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setGainAverages(double)));

    QLabel *GainPositionTextBox = new QLabel;
    GainPositionTextBox->setText("<b>Position:</b>");

    GainPositionLabel  = new QLabel;
    GainPositionSlider = new SliderDouble;
    GainPositionSlider->setPrecision(2);
    GainPositionSlider->setDoubleRange(0, 1.0);
    GainPositionSlider->setOrientation(Qt::Horizontal);
    QObject::connect(GainPositionSlider,  SIGNAL(doubleValueChanged(double)), this, SLOT(setGainPosition(double)));

    QLabel *GainAppearanceTextBox = new QLabel;
    GainAppearanceTextBox->setText("<b>Appearance:</b>");

    GainAppearanceLabel  = new QLabel;
    GainAppearanceSlider = new SliderDouble;
    GainAppearanceSlider->setPrecision(2);
    GainAppearanceSlider->setDoubleRange(0, 1.0);
    GainAppearanceSlider->setOrientation(Qt::Horizontal);
    QObject::connect(GainAppearanceSlider,  SIGNAL(doubleValueChanged(double)), this, SLOT(setGainAppearance(double)));

    QLabel *GainCertaintyTextBox = new QLabel;
    GainCertaintyTextBox->setText("<b>Certainty:</b>");

    GainCertaintyLabel  = new QLabel;
    GainCertaintySlider = new SliderDouble;
    GainCertaintySlider->setPrecision(2);
    GainCertaintySlider->setDoubleRange(0, 1.0);
    GainCertaintySlider->setOrientation(Qt::Horizontal);
    QObject::connect(GainCertaintySlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setGainCertainty(double)));

    // Threshold parameters

    QLabel *ThresholdCircumferenceUpperTextBox = new QLabel;
    QLabel *ThresholdCircumferenceLowerTextBox = new QLabel;

    QLabel *ThresholdAspectRatioUpperTextBox = new QLabel;
    QLabel *ThresholdAspectRatioLowerTextBox = new QLabel;

    QLabel *ThresholdDisplacementUpperTextBox = new QLabel;
    QLabel *ThresholdDisplacementLowerTextBox = new QLabel;

    ThresholdCircumferenceUpperTextBox->setText("<b>Circumference upper change:</b>");
    ThresholdCircumferenceLowerTextBox->setText("<b>Circumference lower change:</b>");

    ThresholdAspectRatioUpperTextBox->setText("<b>Aspect ratio upper change:</b>");
    ThresholdAspectRatioLowerTextBox->setText("<b>Aspect ratio lower change:</b>");

    ThresholdDisplacementUpperTextBox->setText("<b>Displacement upper change:</b>");
    ThresholdDisplacementLowerTextBox->setText("<b>Displacement lower change:</b>");

    ThresholdAspectRatioUpperLabel  = new QLabel;
    ThresholdAspectRatioUpperSlider = new SliderDouble;
    ThresholdAspectRatioUpperSlider->setPrecision(3);
    ThresholdAspectRatioUpperSlider->setDoubleRange(0.0, 0.2);
    ThresholdAspectRatioUpperSlider->setOrientation(Qt::Horizontal);

    ThresholdCircumferenceUpperLabel  = new QLabel;
    ThresholdCircumferenceUpperSlider = new SliderDouble;
    ThresholdCircumferenceUpperSlider->setPrecision(3);
    ThresholdCircumferenceUpperSlider->setDoubleRange(0.0, 0.2);
    ThresholdCircumferenceUpperSlider->setOrientation(Qt::Horizontal);

    ThresholdDisplacementUpperLabel  = new QLabel;
    ThresholdDisplacementUpperSlider = new SliderDouble;
    ThresholdDisplacementUpperSlider->setPrecision(2);
    ThresholdDisplacementUpperSlider->setDoubleRange(0, 20);
    ThresholdDisplacementUpperSlider->setOrientation(Qt::Horizontal);

    ThresholdAspectRatioLowerLabel  = new QLabel;
    ThresholdAspectRatioLowerSlider = new SliderDouble;
    ThresholdAspectRatioLowerSlider->setPrecision(3);
    ThresholdAspectRatioLowerSlider->setDoubleRange(0.0, 0.2);
    ThresholdAspectRatioLowerSlider->setOrientation(Qt::Horizontal);

    ThresholdCircumferenceLowerLabel  = new QLabel;
    ThresholdCircumferenceLowerSlider = new SliderDouble;
    ThresholdCircumferenceLowerSlider->setPrecision(3);
    ThresholdCircumferenceLowerSlider->setDoubleRange(0.0, 0.2);
    ThresholdCircumferenceLowerSlider->setOrientation(Qt::Horizontal);

    ThresholdDisplacementLowerLabel  = new QLabel;
    ThresholdDisplacementLowerSlider = new SliderDouble;
    ThresholdDisplacementLowerSlider->setPrecision(2);
    ThresholdDisplacementLowerSlider->setDoubleRange(0, 20);
    ThresholdDisplacementLowerSlider->setOrientation(Qt::Horizontal);

    QObject::connect(ThresholdAspectRatioUpperSlider,   SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdAspectRatioUpper(double)));
    QObject::connect(ThresholdAspectRatioLowerSlider,   SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdAspectRatioLower(double)));

    QObject::connect(ThresholdCircumferenceUpperSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdCircumferenceUpper(double)));
    QObject::connect(ThresholdCircumferenceLowerSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdCircumferenceLower(double)));

    QObject::connect(ThresholdDisplacementUpperSlider,  SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdDisplacementUpper(double)));
    QObject::connect(ThresholdDisplacementLowerSlider,  SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdDisplacementLower(double)));

    QLabel *ThresholdScoreEdgeTextBox = new QLabel;
    ThresholdScoreEdgeTextBox->setText("<b>Edge score:</b>");

    ThresholdScoreEdgeLabel  = new QLabel;
    ThresholdScoreEdgeSlider = new SliderDouble;
    ThresholdScoreEdgeSlider->setPrecision(2);
    ThresholdScoreEdgeSlider->setDoubleRange(0, 0.5);
    ThresholdScoreEdgeSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdScoreEdgeSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdScoreEdge(double)));

    QLabel *ThresholdScoreFitTextBox = new QLabel;
    ThresholdScoreFitTextBox->setText("<b>Fit score:</b>");

    ThresholdScoreFitLabel  = new QLabel;
    ThresholdScoreFitSlider = new SliderDouble;
    ThresholdScoreFitSlider->setPrecision(2);
    ThresholdScoreFitSlider->setDoubleRange(0, 0.5);
    ThresholdScoreFitSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdScoreFitSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdScoreFit(double)));

    QLabel *ThresholdScoreDiffEdgeTextBox = new QLabel;
    ThresholdScoreDiffEdgeTextBox->setText("<b>Edge score difference:</b>");

    ThresholdScoreDiffEdgeLabel  = new QLabel;
    ThresholdScoreDiffEdgeSlider = new SliderDouble;
    ThresholdScoreDiffEdgeSlider->setPrecision(2);
    ThresholdScoreDiffEdgeSlider->setDoubleRange(0, 1.0);
    ThresholdScoreDiffEdgeSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdScoreDiffEdgeSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdScoreDiffEdge(double)));

    QLabel *ThresholdScoreDiffFitTextBox = new QLabel;
    ThresholdScoreDiffFitTextBox->setText("<b>Fit score difference:</b>");

    ThresholdScoreDiffFitLabel  = new QLabel;
    ThresholdScoreDiffFitSlider = new SliderDouble;
    ThresholdScoreDiffFitSlider->setPrecision(2);
    ThresholdScoreDiffFitSlider->setDoubleRange(0, 0.5);
    ThresholdScoreDiffFitSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdScoreDiffFitSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdScoreDiffFit(double)));

    QLabel *ThresholdFitErrorTextBox = new QLabel;
    ThresholdFitErrorTextBox->setText("<b>Fit error:</b>");

    ThresholdFitErrorLabel  = new QLabel;
    ThresholdFitErrorSlider = new SliderDouble;
    ThresholdFitErrorSlider->setPrecision(2);
    ThresholdFitErrorSlider->setDoubleRange(0, 1.0);
    ThresholdFitErrorSlider->setOrientation(Qt::Horizontal);
    QObject::connect(ThresholdFitErrorSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setThresholdFitError(double)));

    // Miscellaneous parameters

    QLabel *MiscParametersTextBox = new QLabel;
    MiscParametersTextBox->setText("<b>Miscellaneous</b>");
    MiscParametersTextBox->setAlignment(Qt::AlignCenter);

    QLabel *GlintSizeTextBox = new QLabel;
    GlintSizeTextBox->setText("<b>Glint size:</b>");

    GlintSizeLabel  = new QLabel;
    GlintSizeSlider = new QSlider;
    GlintSizeSlider->setRange(0, 10);
    GlintSizeSlider->setOrientation(Qt::Horizontal);
    QObject::connect(GlintSizeSlider, SIGNAL(valueChanged(int)), this, SLOT(setGlintSize(int)));

    QLabel *CurvatureOffsetTextBox = new QLabel;
    CurvatureOffsetTextBox->setText("<b>Curvature offset:</b>");

    CurvatureOffsetLabel  = new QLabel;
    CurvatureOffsetSlider = new SliderDouble;
    CurvatureOffsetSlider->setPrecision(1);
    CurvatureOffsetSlider->setDoubleRange(0, 10);
    CurvatureOffsetSlider->setOrientation(Qt::Horizontal);
    QObject::connect(CurvatureOffsetSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setCurvatureOffset(double)));

    QLabel *WindowLengthEdgeTextBox = new QLabel;
    WindowLengthEdgeTextBox->setText("<b>Edge window length:</b>");

    WindowLengthEdgeLabel  = new QLabel;
    WindowLengthEdgeSlider = new QSlider;
    WindowLengthEdgeSlider->setRange(windowLengthMin, 11);
    WindowLengthEdgeSlider->setOrientation(Qt::Horizontal);
    QObject::connect(WindowLengthEdgeSlider, SIGNAL(valueChanged(int)), this, SLOT(setWindowLengthEdge(int)));

    QLabel *FitEdgeFractionTextBox = new QLabel;
    FitEdgeFractionTextBox->setText("<b>Fit error edge fraction:");

    FitEdgeFractionLabel  = new QLabel;
    FitEdgeFractionSlider = new SliderDouble;
    FitEdgeFractionSlider->setPrecision(2);
    FitEdgeFractionSlider->setDoubleRange(0.0, 0.2);
    FitEdgeFractionSlider->setOrientation(Qt::Horizontal);
    QObject::connect(FitEdgeFractionSlider, SIGNAL(doubleValueChanged(double)), this, SLOT(setFitEdgeFraction(double)));

    QLabel *FitEdgeMaximumTextBox = new QLabel;
    FitEdgeMaximumTextBox->setText("<b>Maximum number of edges:</b>");

    FitEdgeMaximumLabel  = new QLabel;
    FitEdgeMaximumSlider = new QSlider;
    FitEdgeMaximumSlider->setRange(1, 7);
    FitEdgeMaximumSlider->setOrientation(Qt::Horizontal);
    QObject::connect(FitEdgeMaximumSlider, SIGNAL(valueChanged(int)), this, SLOT(setFitEdgeMaximum(int)));

    QLabel *FitMaximumTextBox = new QLabel;
    FitMaximumTextBox->setText("<b>Maximum number of fits:</b>");

    FitMaximumLabel  = new QLabel;
    FitMaximumSlider = new QSlider;
    FitMaximumSlider->setRange(1, 10);
    FitMaximumSlider->setOrientation(Qt::Horizontal);
    QObject::connect(FitMaximumSlider, SIGNAL(valueChanged(int)), this, SLOT(setFitMaximum(int)));

    QLabel *TitleLimitTextBox  = new QLabel;
    QLabel *TitleCannyTextBox  = new QLabel;
    QLabel *TitleLearnTextBox  = new QLabel;
    QLabel *TitleChangeTextBox = new QLabel;
    QLabel *TitleMiscTextBox   = new QLabel;

    TitleLimitTextBox ->setText("<b>Variable limits</b>");
    TitleCannyTextBox ->setText("<b>Canny edge detection</b>");
    TitleLearnTextBox ->setText("<b>Gains</b>");
    TitleChangeTextBox->setText("<b>Thresholds</b>");
    TitleMiscTextBox  ->setText("<b>Miscellaneous</b>");

    QWidget *MainWidget = new QWidget;
    QGridLayout *MainLayout = new QGridLayout(MainWidget);

    // Text labels

    MainLayout->addWidget(CircumferenceMaxTextBox,              1, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CircumferenceMinTextBox,              2, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(AspectRatioMinTextBox,                3, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(CannyThresholdHighTextBox,            5, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyThresholdLowTextBox,             6, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyKernelSizeTextBox,               7, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CannyBlurLevelTextBox,                8, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(GainAveragesTextBox,                  10, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(GainPositionTextBox,                  11, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(GainAppearanceTextBox,                  12, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(GainCertaintyTextBox,                 13, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(ThresholdAspectRatioLowerTextBox,     15, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdAspectRatioUpperTextBox,     16, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdCircumferenceLowerTextBox,   17, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdCircumferenceUpperTextBox,   18, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdDisplacementLowerTextBox,    19, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdDisplacementUpperTextBox,    20, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdScoreEdgeTextBox,            21, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdScoreFitTextBox,             22, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdScoreDiffEdgeTextBox,        23, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdScoreDiffFitTextBox,         24, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(ThresholdFitErrorTextBox,             25, 0, 1, 1, Qt::AlignRight);

    MainLayout->addWidget(GlintSizeTextBox,                     27, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(CurvatureOffsetTextBox,               28, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(WindowLengthEdgeTextBox,              29, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(FitEdgeFractionTextBox,               30, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(FitEdgeMaximumTextBox,                31, 0, 1, 1, Qt::AlignRight);
    MainLayout->addWidget(FitMaximumTextBox,                    32, 0, 1, 1, Qt::AlignRight);

    // Sliders and titles

    MainLayout->addWidget(TitleLimitTextBox,                    0, 1, 1, 1, Qt::AlignCenter);

    MainLayout->addWidget(CircumferenceMaxSlider,               1, 1);
    MainLayout->addWidget(CircumferenceMinSlider,               2, 1);
    MainLayout->addWidget(AspectRatioMinSlider,                 3, 1);

    MainLayout->addWidget(TitleCannyTextBox,                    4, 1, 1, 1, Qt::AlignCenter);

    MainLayout->addWidget(CannyThresholdHighSlider,             5, 1);
    MainLayout->addWidget(CannyThresholdLowSlider,              6, 1);
    MainLayout->addWidget(CannyKernelSizeSlider,                7, 1);
    MainLayout->addWidget(CannyBlurLevelSlider,                 8, 1);

    MainLayout->addWidget(TitleLearnTextBox,                    9, 1, 1, 1, Qt::AlignCenter);

    MainLayout->addWidget(GainAveragesSlider,                   10, 1);
    MainLayout->addWidget(GainPositionSlider,                   11, 1);
    MainLayout->addWidget(GainAppearanceSlider,                   12, 1);
    MainLayout->addWidget(GainCertaintySlider,                  13, 1);

    MainLayout->addWidget(TitleChangeTextBox,                   14, 1, 1, 1, Qt::AlignCenter);

    MainLayout->addWidget(ThresholdAspectRatioLowerSlider,      15, 1);
    MainLayout->addWidget(ThresholdAspectRatioUpperSlider,      16, 1);
    MainLayout->addWidget(ThresholdCircumferenceLowerSlider,    17, 1);
    MainLayout->addWidget(ThresholdCircumferenceUpperSlider,    18, 1);
    MainLayout->addWidget(ThresholdDisplacementLowerSlider,     19, 1);
    MainLayout->addWidget(ThresholdDisplacementUpperSlider,     20, 1);
    MainLayout->addWidget(ThresholdScoreEdgeSlider,             21, 1);
    MainLayout->addWidget(ThresholdScoreFitSlider,              22, 1);
    MainLayout->addWidget(ThresholdScoreDiffEdgeSlider,         23, 1);
    MainLayout->addWidget(ThresholdScoreDiffFitSlider,          24, 1);
    MainLayout->addWidget(ThresholdFitErrorSlider,              25, 1);

    MainLayout->addWidget(TitleMiscTextBox,                     26, 1, 1, 1, Qt::AlignCenter);

    MainLayout->addWidget(GlintSizeSlider,                      27, 1);
    MainLayout->addWidget(CurvatureOffsetSlider,                28, 1);
    MainLayout->addWidget(WindowLengthEdgeSlider,               29, 1);
    MainLayout->addWidget(FitEdgeFractionSlider,                30, 1);
    MainLayout->addWidget(FitEdgeMaximumSlider,                 31, 1);
    MainLayout->addWidget(FitMaximumSlider,                     32, 1);

    // Value labels

    MainLayout->addWidget(CircumferenceMaxLabel,                1, 2);
    MainLayout->addWidget(CircumferenceMinLabel,                2, 2);
    MainLayout->addWidget(  AspectRatioMinLabel,                3, 2);

    MainLayout->addWidget(CannyThresholdHighLabel,              5, 2);
    MainLayout->addWidget(CannyThresholdLowLabel,               6, 2);
    MainLayout->addWidget(CannyKernelSizeLabel,                 7, 2);
    MainLayout->addWidget(CannyBlurLevelLabel,                  8, 2);

    MainLayout->addWidget(GainAveragesLabel,                    10, 2);
    MainLayout->addWidget(GainPositionLabel,                    11, 2);
    MainLayout->addWidget(GainAppearanceLabel,                    12, 2);
    MainLayout->addWidget(GainCertaintyLabel,                   13, 2);

    MainLayout->addWidget(ThresholdAspectRatioLowerLabel,       15, 2);
    MainLayout->addWidget(ThresholdAspectRatioUpperLabel,       16, 2);
    MainLayout->addWidget(ThresholdCircumferenceLowerLabel,     17, 2);
    MainLayout->addWidget(ThresholdCircumferenceUpperLabel,     18, 2);
    MainLayout->addWidget(ThresholdDisplacementLowerLabel,      19, 2);
    MainLayout->addWidget(ThresholdDisplacementUpperLabel,      20, 2);
    MainLayout->addWidget(ThresholdScoreEdgeLabel,              21, 2);
    MainLayout->addWidget(ThresholdScoreFitLabel,               22, 2);
    MainLayout->addWidget(ThresholdScoreDiffEdgeLabel,          23, 2);
    MainLayout->addWidget(ThresholdScoreDiffFitLabel,           24, 2);
    MainLayout->addWidget(ThresholdFitErrorLabel,               25, 2);

    MainLayout->addWidget(GlintSizeLabel,                       27, 2);
    MainLayout->addWidget(CurvatureOffsetLabel,                 28, 2);
    MainLayout->addWidget(WindowLengthEdgeLabel,                29, 2);
    MainLayout->addWidget(FitEdgeFractionLabel,                 30, 2);
    MainLayout->addWidget(FitEdgeMaximumLabel,                  31, 2);
    MainLayout->addWidget(FitMaximumLabel,                      32, 2);

    MainLayout->setColumnStretch(0,1);
    MainLayout->setColumnStretch(1,3);
    MainLayout->setColumnStretch(2,1);

    setLayout(MainLayout);
}

ParameterWidget::~ParameterWidget()
{

}

detectionParameters ParameterWidget::getStructure()
{
    return mDetectionParameters;
}

bool ParameterWidget::getState()
{
  return mDetectionParameters.DETECTION_ON;
}

void ParameterWidget::setState(bool state)
{
    mDetectionParameters.DETECTION_ON = state;
}

void ParameterWidget::setStructure(detectionParameters mEyePropertiesParametersNew)
{
    mDetectionParameters = mEyePropertiesParametersNew;
    this->reset();
}

void ParameterWidget::reset()
{
    // Size and shape limits

    CircumferenceMaxSlider->setDoubleValue(mDetectionParameters.thresholdCircumferenceMax);
    CircumferenceMaxLabel ->setText(QString::number(mDetectionParameters.thresholdCircumferenceMax, 'f', 1));

    CircumferenceMinSlider->setDoubleValue(mDetectionParameters.thresholdCircumferenceMin);
    CircumferenceMinLabel ->setText(QString::number(mDetectionParameters.thresholdCircumferenceMin, 'f', 1));

    AspectRatioMinSlider->setDoubleValue(mDetectionParameters.thresholdAspectRatioMin);
    AspectRatioMinLabel ->setText(QString::number(mDetectionParameters.thresholdAspectRatioMin, 'f', 2));

    // Canny edge detection

    CannyThresholdHighSlider->setDoubleValue(mDetectionParameters.cannyThresholdHigh);
    CannyThresholdHighLabel ->setText(QString::number(mDetectionParameters.cannyThresholdHigh, 'f', 1));

    CannyThresholdLowSlider->setDoubleValue(mDetectionParameters.cannyThresholdLow);
    CannyThresholdLowLabel ->setText(QString::number(mDetectionParameters.cannyThresholdLow, 'f', 1));

    CannyBlurLevelSlider->setValue(mDetectionParameters.cannyBlurLevel);
    CannyBlurLevelLabel ->setText(QString::number(mDetectionParameters.cannyBlurLevel));

    CannyKernelSizeSlider->setValue(ceil(0.5 * mDetectionParameters.cannyKernelSize));
    CannyKernelSizeLabel ->setText(QString::number(mDetectionParameters.cannyKernelSize));

    // Gains

    GainAveragesSlider  ->setDoubleValue(mDetectionParameters.gainAverages);
    GainPositionSlider  ->setDoubleValue(mDetectionParameters.gainPosition);
    GainAppearanceSlider->setDoubleValue(mDetectionParameters.gainAppearance);
    GainCertaintySlider ->setDoubleValue(mDetectionParameters.gainCertainty);

    GainAveragesLabel  ->setText(QString::number(mDetectionParameters.gainAverages,   'f', 3));
    GainPositionLabel  ->setText(QString::number(mDetectionParameters.gainPosition,   'f', 2));
    GainAppearanceLabel->setText(QString::number(mDetectionParameters.gainAppearance, 'f', 2));
    GainCertaintyLabel ->setText(QString::number(mDetectionParameters.gainCertainty,  'f', 2));

    // Thresholds

    ThresholdCircumferenceLowerSlider->setDoubleValue(mDetectionParameters.thresholdChangeCircumferenceLower);
    ThresholdCircumferenceUpperSlider->setDoubleValue(mDetectionParameters.thresholdChangeCircumferenceUpper);

    ThresholdCircumferenceLowerLabel ->setText(QString::number(mDetectionParameters.thresholdChangeCircumferenceLower, 'f', 3));
    ThresholdCircumferenceUpperLabel ->setText(QString::number(mDetectionParameters.thresholdChangeCircumferenceUpper, 'f', 3));

    ThresholdAspectRatioLowerSlider->setDoubleValue(mDetectionParameters.thresholdChangeAspectRatioLower);
    ThresholdAspectRatioUpperSlider->setDoubleValue(mDetectionParameters.thresholdChangeAspectRatioUpper);

    ThresholdAspectRatioLowerLabel ->setText(QString::number(mDetectionParameters.thresholdChangeAspectRatioLower, 'f', 3));
    ThresholdAspectRatioUpperLabel ->setText(QString::number(mDetectionParameters.thresholdChangeAspectRatioUpper, 'f', 3));

    ThresholdDisplacementLowerSlider->setDoubleValue(mDetectionParameters.thresholdChangePositionLower);
    ThresholdDisplacementUpperSlider->setDoubleValue(mDetectionParameters.thresholdChangePositionUpper);

    ThresholdDisplacementLowerLabel ->setText(QString::number(mDetectionParameters.thresholdChangePositionLower, 'f', 2));
    ThresholdDisplacementUpperLabel ->setText(QString::number(mDetectionParameters.thresholdChangePositionUpper, 'f', 2));

    ThresholdScoreEdgeSlider->setDoubleValue(mDetectionParameters.thresholdScoreEdge);
    ThresholdScoreEdgeLabel ->setText(QString::number(mDetectionParameters.thresholdScoreEdge, 'f', 2));

    ThresholdScoreFitSlider->setDoubleValue(mDetectionParameters.thresholdScoreFit);
    ThresholdScoreFitLabel ->setText(QString::number(mDetectionParameters.thresholdScoreFit, 'f', 2));

    ThresholdScoreDiffEdgeSlider->setDoubleValue(mDetectionParameters.thresholdScoreDiffEdge);
    ThresholdScoreDiffEdgeLabel ->setText(QString::number(mDetectionParameters.thresholdScoreDiffEdge, 'f', 2));

    ThresholdScoreDiffFitSlider->setDoubleValue(mDetectionParameters.thresholdScoreDiffFit);
    ThresholdScoreDiffFitLabel ->setText(QString::number(mDetectionParameters.thresholdScoreDiffFit, 'f', 2));

    ThresholdFitErrorSlider->setDoubleValue(mDetectionParameters.thresholdFitError);
    ThresholdFitErrorLabel ->setText(QString::number(mDetectionParameters.thresholdFitError, 'f', 2));

    // Misc

    GlintSizeSlider->setValue(round(0.5 * mDetectionParameters.glintWdth));
    GlintSizeLabel->setText(QString::number(mDetectionParameters.glintWdth));

    CurvatureOffsetSlider->setDoubleValue(mDetectionParameters.curvatureOffset);
    CurvatureOffsetLabel ->setText(QString::number(mDetectionParameters.curvatureOffset, 'f', 1));

    WindowLengthEdgeSlider->setValue(mDetectionParameters.windowLengthEdge);
    WindowLengthEdgeLabel ->setText(QString::number(mDetectionParameters.windowLengthEdge));

    FitEdgeFractionSlider->setDoubleValue(mDetectionParameters.fitEdgeFraction);
    FitEdgeFractionLabel ->setText(QString::number(mDetectionParameters.fitEdgeFraction, 'f', 2));

    FitEdgeMaximumSlider->setValue(mDetectionParameters.fitEdgeMaximum);
    FitEdgeMaximumLabel ->setText(QString::number(mDetectionParameters.fitEdgeMaximum));

    FitMaximumSlider->setValue(mDetectionParameters.fitMaximum);
    FitMaximumLabel ->setText(QString::number(mDetectionParameters.fitMaximum));
}

void ParameterWidget::setCircumferenceMin(double value)
{
    if (mDetectionParameters.thresholdCircumferenceMax < value)
    {   CircumferenceMaxSlider->setDoubleValue(value); }

    mDetectionParameters.thresholdCircumferenceMin = value;
    CircumferenceMinLabel->setText(QString::number(mDetectionParameters.thresholdCircumferenceMin, 'f', 1));
}

void ParameterWidget::setCircumferenceMax(double value)
{
    if (mDetectionParameters.thresholdCircumferenceMin > value)
    {   CircumferenceMinSlider->setDoubleValue(value); }

    mDetectionParameters.thresholdCircumferenceMax = value;
    CircumferenceMaxLabel->setText(QString::number(mDetectionParameters.thresholdCircumferenceMax, 'f', 1));
}

void ParameterWidget::setAspectRatioMin(double value)
{
    mDetectionParameters.thresholdAspectRatioMin = value;
    AspectRatioMinLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setFitEdgeFraction(double value)
{
    mDetectionParameters.fitEdgeFraction = value;
    FitEdgeFractionLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setCannyThresholdLow(double value)
{
    if (value > mDetectionParameters.cannyThresholdHigh)
    {   CannyThresholdHighSlider->setDoubleValue(value); }

    mDetectionParameters.cannyThresholdLow = value;
    CannyThresholdLowLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setCannyThresholdHigh(double value)
{
    if (value < mDetectionParameters.cannyThresholdLow)
    {   CannyThresholdLowSlider->setDoubleValue(value); }

    mDetectionParameters.cannyThresholdHigh = value;
    CannyThresholdHighLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setCannyKernelSize(int value)
{
    int newValue = 2 * value - 1;
    mDetectionParameters.cannyKernelSize = newValue;
    CannyKernelSizeLabel->setText(QString::number(newValue));
}

void ParameterWidget::setCannyBlurLevel(int value)
{
    mDetectionParameters.cannyBlurLevel = value;
    CannyBlurLevelLabel->setText(QString::number(value));
}

void ParameterWidget::setGainAverages(double value)
{
    mDetectionParameters.gainAverages = value;
    GainAveragesLabel->setText(QString::number(value, 'f', 3));
}

void ParameterWidget::setGainPosition(double value)
{
    mDetectionParameters.gainPosition = value;
    GainPositionLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setGainAppearance(double value)
{
    mDetectionParameters.gainAppearance = value;
    GainAppearanceLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setGainCertainty(double value)
{
    mDetectionParameters.gainCertainty = value;
    GainCertaintyLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdCircumferenceLower(double value)
{
    if (value > mDetectionParameters.thresholdChangeCircumferenceUpper)
    {   ThresholdCircumferenceUpperSlider->setDoubleValue(value); }

    mDetectionParameters.thresholdChangeCircumferenceLower = value;
    ThresholdCircumferenceLowerLabel->setText(QString::number(value, 'f', 3));
}

void ParameterWidget::setThresholdCircumferenceUpper(double value)
{
    if (value < mDetectionParameters.thresholdChangeCircumferenceLower)
    {   ThresholdCircumferenceLowerSlider->setDoubleValue(value); }

    mDetectionParameters.thresholdChangeCircumferenceUpper = value;
    ThresholdCircumferenceUpperLabel->setText(QString::number(value, 'f', 3));
}

void ParameterWidget::setThresholdAspectRatioLower(double value)
{
    if (value > mDetectionParameters.thresholdChangeAspectRatioUpper)
    {   ThresholdAspectRatioUpperSlider->setDoubleValue(value); }

    mDetectionParameters.thresholdChangeAspectRatioLower = value;
    ThresholdAspectRatioLowerLabel->setText(QString::number(value, 'f', 3));
}

void ParameterWidget::setThresholdAspectRatioUpper(double value)
{
    if (value < mDetectionParameters.thresholdChangeAspectRatioLower)
    {   ThresholdAspectRatioLowerSlider->setDoubleValue(value); }

    mDetectionParameters.thresholdChangeAspectRatioUpper = value;
    ThresholdAspectRatioUpperLabel->setText(QString::number(value, 'f', 3));
}

void ParameterWidget::setThresholdDisplacementLower(double value)
{
    if (value > mDetectionParameters.thresholdChangePositionUpper)
    {   ThresholdDisplacementUpperSlider->setDoubleValue(value); }

    mDetectionParameters.thresholdChangePositionLower = value;
    ThresholdDisplacementLowerLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdDisplacementUpper(double value)
{
    if (value < mDetectionParameters.thresholdChangePositionLower)
    {   ThresholdDisplacementLowerSlider->setDoubleValue(value); }

    mDetectionParameters.thresholdChangePositionUpper = value;
    ThresholdDisplacementUpperLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdScoreEdge(double value)
{
    mDetectionParameters.thresholdScoreEdge = value;
    ThresholdScoreEdgeLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdScoreFit(double value)
{
    mDetectionParameters.thresholdScoreFit = value;
    ThresholdScoreFitLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdScoreDiffEdge(double value)
{
    mDetectionParameters.thresholdScoreDiffEdge = value;
    ThresholdScoreDiffEdgeLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setThresholdScoreDiffFit(double value)
{
    mDetectionParameters.thresholdScoreDiffFit = value;
    ThresholdScoreDiffFitLabel->setText(QString::number(value, 'f', 2));
}

void ParameterWidget::setGlintSize(int value)
{
    int newValue = 2 * value;
    mDetectionParameters.glintWdth = newValue;
    GlintSizeLabel->setText(QString::number(newValue));
}

void ParameterWidget::setCurvatureOffset(double value)
{
    mDetectionParameters.curvatureOffset = value;
    CurvatureOffsetLabel->setText(QString::number(value, 'f', 1));
}

void ParameterWidget::setWindowLengthEdge(int value)
{
    mDetectionParameters.windowLengthEdge = value;
    WindowLengthEdgeLabel->setText(QString::number(value));
}

void ParameterWidget::setFitEdgeMaximum(int value)
{
    mDetectionParameters.fitEdgeMaximum = value;
    FitEdgeMaximumLabel->setText(QString::number(value));
}

void ParameterWidget::setFitMaximum(int value)
{
    mDetectionParameters.fitMaximum = value;
    FitMaximumLabel->setText(QString::number(value));
}

void ParameterWidget::setThresholdFitError(double value)
{
    mDetectionParameters.thresholdFitError = value;
    ThresholdFitErrorLabel->setText(QString::number(value, 'f', 2));
}
