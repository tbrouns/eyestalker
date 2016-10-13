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


#include "mainwindow.h"
#include "source/startupwindow.h"

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
    mMainWindow.setWindowTitle("Eye-tracking application");

    StartUpWindow *mStartUpWindow = new StartUpWindow;

    QString styleSheetLocation = QString(":qdarkstyle/style.qss");
    QFile styleSheetFile(styleSheetLocation);

    if (!styleSheetFile.exists())
    {
        printf("Unable to set stylesheet, file not found\n");
    }
    else
    {
        {
            QTextStream ts(&styleSheetFile);
            styleSheetFile.open(QFile::ReadOnly | QFile::Text);
            mMainWindow.setStyleSheet(ts.readAll());
        }

        {
            QFile newStyleSheetFile(styleSheetLocation);
            QTextStream ts(&newStyleSheetFile);
            newStyleSheetFile.open(QFile::ReadOnly | QFile::Text);
            mStartUpWindow->setStyleSheet(ts.readAll());
        }
    }

    // First-time start-up window

    if (!boost::filesystem::exists("config.ini")) // Check if global config.ini file is present in folder
    {
        mStartUpWindow->setWindowTitle("Start-up window");
        mStartUpWindow->exec();

        if(!mStartUpWindow->getStatus())
        {
            return 0;
        }
    }
    else
    {
        delete mStartUpWindow;
    }

    mMainWindow.show();

    return app.exec();
}
