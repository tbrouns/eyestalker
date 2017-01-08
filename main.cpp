//  Copyright (C) 2016  Terence Brouns

//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>

#include "headers/mainwindow.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <QApplication>

#include <iostream>     // for console output
#include <stdio.h>      // for sprintf()
#include <string>       // for std::string

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    MainWindow mMainWindow;
    mMainWindow.setWindowTitle("EyeStalker");

    QString styleSheetLocation = QString(":qdarkstyle/style.qss");
    QFile styleSheetFile(styleSheetLocation);

    if (!styleSheetFile.exists()) { printf("Unable to set stylesheet, file not found\n"); }
    else
    {
        QTextStream ts(&styleSheetFile);
        styleSheetFile.open(QFile::ReadOnly | QFile::Text);
        mMainWindow.setStyleSheet(ts.readAll());
    }

    mMainWindow.show();

    return app.exec();
}
