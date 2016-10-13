#include "mainwindow.h"

void MainWindow::CalibrationGridChange(int value)
{
    gridSize = value;
    CalibrationPlotSpinBox->setMaximum(gridSize);
    setCalibrationDataLabel(1);
}

void MainWindow::CalibrationStartAuto()
{
    if (Parameters::CAMERA_RUNNING)
    {
        // Create folders

        {
            std::stringstream directory_name;
            directory_name << dataDirectory << "/" << currentDate << "/" << (NameInputLineEdit->text()).toStdString();
            mkdir(directory_name.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }

        {
            std::stringstream directory_name;
            directory_name << dataDirectory << "/" << currentDate << "/" << (NameInputLineEdit->text()).toStdString() << "/calibration/";
            mkdir(directory_name.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        }

        Parameters::CALIBRATION_MOUSE_PRESS_LOCK = true;
        CALIBRATION_RUNNING = true;

        for (int i = 0; i < gridSize && CALIBRATION_RUNNING; i++)
        {
            // create folders for calibration point

            {
                std::stringstream directory_name;
                directory_name << dataDirectory << "/" << currentDate << "/" << (NameInputLineEdit->text()).toStdString() << "/calibration/point_" << i;
                mkdir(directory_name.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            }

            {
                std::stringstream directory_name;
                directory_name << dataDirectory << "/" << currentDate << "/" << (NameInputLineEdit->text()).toStdString() << "/calibration/point_" << i << "/images/";
                mkdir(directory_name.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            }

            {
                std::stringstream directory_name;
                directory_name << dataDirectory << "/" << currentDate << "/" << (NameInputLineEdit->text()).toStdString() << "/calibration/point_" << i << "/images/raw/";
                mkdir(directory_name.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            }

            ExpWindow->CalibrationPointShow(i, Parameters::calibrationPointColour);
            CalWidget->PointColourFill(i, Qt::red); // Red colour fill

            msWait(Parameters::calibrationTimerLength);

            CalWidget->PointColourFill(i, Qt::green); // Green colour fill
            CalibrationRecordPoint(i);

            std::unique_lock<std::mutex> lck(Parameters::calibration_mtx);
            while (Parameters::CALIBRATION_POINT_RECORDING) Parameters::calibration_cv.wait(lck);

            ExpWindow->CalibrationPointShow(i, Qt::black);
            CalWidget->PointColourFill(i); // No colour fill

            msWait(round(0.5 * Parameters::calibrationTimerLength));

            setCalibrationDataLabel(i + 1); // Set point index and show plot
        }

        Parameters::CALIBRATION_MOUSE_PRESS_LOCK = false;
    }
}

void MainWindow::CalibrationTerminateAuto()
{
    CALIBRATION_RUNNING = false;
}

void MainWindow::CalibrationRecordPoint(int index)
{
    calibrationIndex = index;

    frameCount = 0;

    Parameters::CALIBRATION_POINT_RECORDING = true;
}

std::vector<double> MainWindow::CalibrationGetScreenPositions(int coordinateIndex, int numberOfCalibrationPoints)
{
    std::vector<double> screenCoordinates(numberOfCalibrationPoints); // positions of calibration points in screen coordinates

    switch(numberOfCalibrationPoints)
    {

    case 4:

        switch(coordinateIndex)
        {

        case 0: // x

            screenCoordinates[0] = -1;
            screenCoordinates[1] =  1;
            screenCoordinates[2] = -1;
            screenCoordinates[3] =  1;

            break;

        case 1: // y

            screenCoordinates[0] =  1;
            screenCoordinates[1] =  1;
            screenCoordinates[2] = -1;
            screenCoordinates[3] = -1;

            break;
        }

        break;

    case 9:

        switch(coordinateIndex)
        {

        case 0: // x

            screenCoordinates[0] = -1;
            screenCoordinates[1] =  0;
            screenCoordinates[2] =  1;
            screenCoordinates[3] = -1;
            screenCoordinates[4] =  0;
            screenCoordinates[5] =  1;
            screenCoordinates[6] = -1;
            screenCoordinates[7] =  0;
            screenCoordinates[8] =  1;

            break;

        case 1: // y

            screenCoordinates[0] =  1;
            screenCoordinates[1] =  1;
            screenCoordinates[2] =  1;
            screenCoordinates[3] =  0;
            screenCoordinates[4] =  0;
            screenCoordinates[5] =  0;
            screenCoordinates[6] = -1;
            screenCoordinates[7] = -1;
            screenCoordinates[8] = -1;

            break;
        }

        break;

    case 16:

        switch(coordinateIndex)
        {

        case 0: // x

            screenCoordinates[0]  = -(3/3);
            screenCoordinates[1]  = -(1/3);
            screenCoordinates[2]  =  (1/3);
            screenCoordinates[3]  =  (3/3);
            screenCoordinates[4]  = -(3/3);
            screenCoordinates[5]  = -(1/3);
            screenCoordinates[6]  =  (1/3);
            screenCoordinates[7]  =  (3/3);
            screenCoordinates[8]  = -(3/3);
            screenCoordinates[9]  = -(1/3);
            screenCoordinates[10] =  (1/3);
            screenCoordinates[11] =  (3/3);
            screenCoordinates[12] = -(3/3);
            screenCoordinates[13] = -(1/3);
            screenCoordinates[14] =  (1/3);
            screenCoordinates[15] =  (3/3);

            break;

        case 1: // y

            screenCoordinates[0]  =  (3/3);
            screenCoordinates[1]  =  (3/3);
            screenCoordinates[2]  =  (3/3);
            screenCoordinates[3]  =  (3/3);
            screenCoordinates[4]  =  (1/3);
            screenCoordinates[5]  =  (1/3);
            screenCoordinates[6]  =  (1/3);
            screenCoordinates[7]  =  (1/3);
            screenCoordinates[8]  = -(1/3);
            screenCoordinates[9]  = -(1/3);
            screenCoordinates[10] = -(1/3);
            screenCoordinates[11] = -(1/3);
            screenCoordinates[12] = -(3/3);
            screenCoordinates[13] = -(3/3);
            screenCoordinates[14] = -(3/3);
            screenCoordinates[15] = -(3/3);

            break;
        }

        break;

    case 25:

        switch(coordinateIndex)
        {

        case 0: // x

            screenCoordinates[0]  = -1.0;
            screenCoordinates[1]  = -0.5;
            screenCoordinates[2]  =  0.0;
            screenCoordinates[3]  =  0.5;
            screenCoordinates[4]  =  1.0;
            screenCoordinates[5]  = -1.0;
            screenCoordinates[6]  = -0.5;
            screenCoordinates[7]  =  0.0;
            screenCoordinates[8]  =  0.5;
            screenCoordinates[9]  =  1.0;
            screenCoordinates[10] = -1.0;
            screenCoordinates[11] = -0.5;
            screenCoordinates[12] =  0.0;
            screenCoordinates[13] =  0.5;
            screenCoordinates[14] =  1.0;
            screenCoordinates[15] = -1.0;
            screenCoordinates[16] = -0.5;
            screenCoordinates[17] =  0.0;
            screenCoordinates[18] =  0.5;
            screenCoordinates[19] =  1.0;
            screenCoordinates[20] = -1.0;
            screenCoordinates[21] = -0.5;
            screenCoordinates[22] =  0.0;
            screenCoordinates[23] =  0.5;
            screenCoordinates[24] =  1.0;

            break;

        case 1: // y

            screenCoordinates[0]  =  1.0;
            screenCoordinates[1]  =  1.0;
            screenCoordinates[2]  =  1.0;
            screenCoordinates[3]  =  1.0;
            screenCoordinates[4]  =  1.0;
            screenCoordinates[5]  =  0.5;
            screenCoordinates[6]  =  0.5;
            screenCoordinates[7]  =  0.5;
            screenCoordinates[8]  =  0.5;
            screenCoordinates[9]  =  0.5;
            screenCoordinates[10] =  0.0;
            screenCoordinates[11] =  0.0;
            screenCoordinates[12] =  0.0;
            screenCoordinates[13] =  0.0;
            screenCoordinates[14] =  0.0;
            screenCoordinates[15] = -0.5;
            screenCoordinates[16] = -0.5;
            screenCoordinates[17] = -0.5;
            screenCoordinates[18] = -0.5;
            screenCoordinates[19] = -0.5;
            screenCoordinates[20] = -1.0;
            screenCoordinates[21] = -1.0;
            screenCoordinates[22] = -1.0;
            screenCoordinates[23] = -1.0;
            screenCoordinates[24] = -1.0;

            break;

        }

        break;
    }

    return screenCoordinates;
}

std::vector<double> MainWindow::CalibrationGetTransformVector(const std::vector<double>& calibrationXPos, const std::vector<double>& calibrationYPos)
{
    int numberOfCalibrationPoints = calibrationXPos.size();

    std::vector<double> transformVector;

    std::vector<double> screenXPos = CalibrationGetScreenPositions(0, numberOfCalibrationPoints);
    std::vector<double> screenYPos = CalibrationGetScreenPositions(1, numberOfCalibrationPoints);

    Eigen::Matrix4d calibrationMatrix;
    Eigen::Matrix4d screenMatrix;
    Eigen::Matrix4d transformMatrix;

    switch(numberOfCalibrationPoints)
    {

    case 4:

        // pupil matrix

        calibrationMatrix(0,0) = calibrationXPos[0];
        calibrationMatrix(0,1) = calibrationXPos[1];
        calibrationMatrix(0,2) = calibrationXPos[2];
        calibrationMatrix(0,3) = calibrationXPos[3];

        calibrationMatrix(1,0) = calibrationYPos[0];
        calibrationMatrix(1,1) = calibrationYPos[1];
        calibrationMatrix(1,2) = calibrationYPos[2];
        calibrationMatrix(1,3) = calibrationYPos[3];

        calibrationMatrix(2,0) = calibrationXPos[0] * calibrationYPos[0];
        calibrationMatrix(2,1) = calibrationXPos[1] * calibrationYPos[1];
        calibrationMatrix(2,2) = calibrationXPos[2] * calibrationYPos[2];
        calibrationMatrix(2,3) = calibrationXPos[3] * calibrationYPos[3];

        calibrationMatrix(3,0) = 1;
        calibrationMatrix(3,1) = 1;
        calibrationMatrix(3,2) = 1;
        calibrationMatrix(3,3) = 1;

        // screen matrix

        screenMatrix = Eigen::Matrix4d::Zero();

        for (int i = 0; i < 4; i++)
        {
            screenMatrix(0, i) = screenXPos[i];
            screenMatrix(1, i) = screenYPos[i];
        }

        // calculate transform matrix

        transformMatrix = screenMatrix * calibrationMatrix.inverse();

        // convert to vector

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                transformVector.push_back(transformMatrix(i,j));
            }
        }

        break;

    case 9:

        for (int iScreenArea = 0; iScreenArea < 4; iScreenArea++) // iterate over all areas
        {
            int iScreenAreaTopLeft = iScreenArea + floor(iScreenArea / 2);
            int iScreenAreaTopRght = iScreenAreaTopLeft + 1;
            int iScreenAreaBtmLeft = iScreenAreaTopLeft + 3;
            int iScreenAreaBtmRght = iScreenAreaTopLeft + 4;

            // pupil matrix

            calibrationMatrix(0,0) = calibrationXPos[iScreenAreaTopLeft];
            calibrationMatrix(0,1) = calibrationXPos[iScreenAreaTopRght];
            calibrationMatrix(0,2) = calibrationXPos[iScreenAreaBtmLeft];
            calibrationMatrix(0,3) = calibrationXPos[iScreenAreaBtmRght];

            calibrationMatrix(1,0) = calibrationYPos[iScreenAreaTopLeft];
            calibrationMatrix(1,1) = calibrationYPos[iScreenAreaTopRght];
            calibrationMatrix(1,2) = calibrationYPos[iScreenAreaBtmLeft];
            calibrationMatrix(1,3) = calibrationYPos[iScreenAreaBtmRght];

            calibrationMatrix(2,0) = calibrationXPos[iScreenAreaTopLeft] * calibrationYPos[iScreenAreaTopLeft];
            calibrationMatrix(2,1) = calibrationXPos[iScreenAreaTopRght] * calibrationYPos[iScreenAreaTopRght];
            calibrationMatrix(2,2) = calibrationXPos[iScreenAreaBtmLeft] * calibrationYPos[iScreenAreaBtmLeft];
            calibrationMatrix(2,3) = calibrationXPos[iScreenAreaBtmRght] * calibrationYPos[iScreenAreaBtmRght];

            calibrationMatrix(3,0) = 1;
            calibrationMatrix(3,1) = 1;
            calibrationMatrix(3,2) = 1;
            calibrationMatrix(3,3) = 1;

            // screen matrix

            screenMatrix = Eigen::Matrix4d::Zero();

            screenMatrix(0, 0) = screenXPos[iScreenAreaTopLeft];
            screenMatrix(0, 1) = screenXPos[iScreenAreaTopRght];
            screenMatrix(0, 2) = screenXPos[iScreenAreaBtmLeft];
            screenMatrix(0, 3) = screenXPos[iScreenAreaBtmRght];

            screenMatrix(1, 0) = screenYPos[iScreenAreaTopLeft];
            screenMatrix(1, 1) = screenYPos[iScreenAreaTopRght];
            screenMatrix(1, 2) = screenYPos[iScreenAreaBtmLeft];
            screenMatrix(1, 3) = screenYPos[iScreenAreaBtmRght];

            // calculate transform matrix

            transformMatrix = screenMatrix * calibrationMatrix.inverse();

            // convert to vector

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    transformVector.push_back(transformMatrix(i,j));
                }
            }
        }

        break;

    case 16:

        for (int iScreenArea = 0; iScreenArea < 9; iScreenArea++) // iterate over all areas
        {
            int iScreenAreaTopLeft = iScreenArea + floor(iScreenArea / 3);
            int iScreenAreaTopRght = iScreenAreaTopLeft + 1;
            int iScreenAreaBtmLeft = iScreenAreaTopLeft + 4;
            int iScreenAreaBtmRght = iScreenAreaTopLeft + 5;

            // pupil matrix

            calibrationMatrix(0,0) = calibrationXPos[iScreenAreaTopLeft];
            calibrationMatrix(0,1) = calibrationXPos[iScreenAreaTopRght];
            calibrationMatrix(0,2) = calibrationXPos[iScreenAreaBtmLeft];
            calibrationMatrix(0,3) = calibrationXPos[iScreenAreaBtmRght];

            calibrationMatrix(1,0) = calibrationYPos[iScreenAreaTopLeft];
            calibrationMatrix(1,1) = calibrationYPos[iScreenAreaTopRght];
            calibrationMatrix(1,2) = calibrationYPos[iScreenAreaBtmLeft];
            calibrationMatrix(1,3) = calibrationYPos[iScreenAreaBtmRght];

            calibrationMatrix(2,0) = calibrationXPos[iScreenAreaTopLeft] * calibrationYPos[iScreenAreaTopLeft];
            calibrationMatrix(2,1) = calibrationXPos[iScreenAreaTopRght] * calibrationYPos[iScreenAreaTopRght];
            calibrationMatrix(2,2) = calibrationXPos[iScreenAreaBtmLeft] * calibrationYPos[iScreenAreaBtmLeft];
            calibrationMatrix(2,3) = calibrationXPos[iScreenAreaBtmRght] * calibrationYPos[iScreenAreaBtmRght];

            calibrationMatrix(3,0) = 1;
            calibrationMatrix(3,1) = 1;
            calibrationMatrix(3,2) = 1;
            calibrationMatrix(3,3) = 1;

            // screen matrix

            screenMatrix = Eigen::Matrix4d::Zero();

            screenMatrix(0, 0) = screenXPos[iScreenAreaTopLeft];
            screenMatrix(0, 1) = screenXPos[iScreenAreaTopRght];
            screenMatrix(0, 2) = screenXPos[iScreenAreaBtmLeft];
            screenMatrix(0, 3) = screenXPos[iScreenAreaBtmRght];

            screenMatrix(1, 0) = screenYPos[iScreenAreaTopLeft];
            screenMatrix(1, 1) = screenYPos[iScreenAreaTopRght];
            screenMatrix(1, 2) = screenYPos[iScreenAreaBtmLeft];
            screenMatrix(1, 3) = screenYPos[iScreenAreaBtmRght];

            // calculate transform matrix

            transformMatrix = screenMatrix * calibrationMatrix.inverse();

            // convert to vector

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    transformVector.push_back(transformMatrix(i, j));
                }
            }
        }

        break;

    case 25:

        for (int iScreenArea = 0; iScreenArea < 16; iScreenArea++) // iterate over all areas
        {
            int iScreenAreaTopLeft = iScreenArea + floor(iScreenArea / 4);
            int iScreenAreaTopRght = iScreenAreaTopLeft + 1;
            int iScreenAreaBtmLeft = iScreenAreaTopLeft + 5;
            int iScreenAreaBtmRght = iScreenAreaTopLeft + 6;

            // pupil matrix

            calibrationMatrix(0,0) = calibrationXPos[iScreenAreaTopLeft];
            calibrationMatrix(0,1) = calibrationXPos[iScreenAreaTopRght];
            calibrationMatrix(0,2) = calibrationXPos[iScreenAreaBtmLeft];
            calibrationMatrix(0,3) = calibrationXPos[iScreenAreaBtmRght];

            calibrationMatrix(1,0) = calibrationYPos[iScreenAreaTopLeft];
            calibrationMatrix(1,1) = calibrationYPos[iScreenAreaTopRght];
            calibrationMatrix(1,2) = calibrationYPos[iScreenAreaBtmLeft];
            calibrationMatrix(1,3) = calibrationYPos[iScreenAreaBtmRght];

            calibrationMatrix(2,0) = calibrationXPos[iScreenAreaTopLeft] * calibrationYPos[iScreenAreaTopLeft];
            calibrationMatrix(2,1) = calibrationXPos[iScreenAreaTopRght] * calibrationYPos[iScreenAreaTopRght];
            calibrationMatrix(2,2) = calibrationXPos[iScreenAreaBtmLeft] * calibrationYPos[iScreenAreaBtmLeft];
            calibrationMatrix(2,3) = calibrationXPos[iScreenAreaBtmRght] * calibrationYPos[iScreenAreaBtmRght];

            calibrationMatrix(3,0) = 1;
            calibrationMatrix(3,1) = 1;
            calibrationMatrix(3,2) = 1;
            calibrationMatrix(3,3) = 1;

            // screen matrix

            screenMatrix = Eigen::Matrix4d::Zero();

            screenMatrix(0, 0) = screenXPos[iScreenAreaTopLeft];
            screenMatrix(0, 1) = screenXPos[iScreenAreaTopRght];
            screenMatrix(0, 2) = screenXPos[iScreenAreaBtmLeft];
            screenMatrix(0, 3) = screenXPos[iScreenAreaBtmRght];

            screenMatrix(1, 0) = screenYPos[iScreenAreaTopLeft];
            screenMatrix(1, 1) = screenYPos[iScreenAreaTopRght];
            screenMatrix(1, 2) = screenYPos[iScreenAreaBtmLeft];
            screenMatrix(1, 3) = screenYPos[iScreenAreaBtmRght];

            // calculate transform matrix

            transformMatrix = screenMatrix * calibrationMatrix.inverse();

            // convert to vector

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    transformVector.push_back(transformMatrix(i, j));
                }
            }
        }

        break;
    }

    return transformVector;
}

std::vector<double> MainWindow::CalibrationGeometricTransform(const std::vector<double>& transformVector, const std::vector<double>& calibrationXPos, const std::vector<double>& calibrationYPos, double rawTrackingXPos, double rawTrackingYPos)
{
    int numberOfCalibrationPoints = calibrationXPos.size();

    std::vector<double> newTrackingCoordinates(2); // coordinates

    int screenAreaIndex = 0;

    switch(numberOfCalibrationPoints)
    {

    case 4:

        newTrackingCoordinates[0] = transformVector[0] * rawTrackingXPos + transformVector[1] * rawTrackingYPos + transformVector[2] * rawTrackingXPos * rawTrackingYPos + transformVector[3];
        newTrackingCoordinates[1] = transformVector[4] * rawTrackingXPos + transformVector[5] * rawTrackingYPos + transformVector[6] * rawTrackingXPos * rawTrackingYPos + transformVector[7];

        break;

    case 9:

        if (rawTrackingYPos <= calibrationYPos[4]) // first row
        {
            if (rawTrackingXPos <= calibrationXPos[4])
            {
                screenAreaIndex = 0;
            }
            else
            {
                screenAreaIndex = 1;
            }
        }
        else if (rawTrackingYPos > calibrationYPos[4]) // second row
        {
            if (rawTrackingXPos <= calibrationXPos[4])
            {
                screenAreaIndex = 2;
            }
            else
            {
                screenAreaIndex = 3;
            }
        }

        screenAreaIndex = screenAreaIndex * 8;

        newTrackingCoordinates[0] = transformVector[0 + screenAreaIndex] * rawTrackingXPos + transformVector[1 + screenAreaIndex] * rawTrackingYPos + transformVector[2 + screenAreaIndex] * rawTrackingXPos * rawTrackingYPos + transformVector[3 + screenAreaIndex];
        newTrackingCoordinates[1] = transformVector[4 + screenAreaIndex] * rawTrackingXPos + transformVector[5 + screenAreaIndex] * rawTrackingYPos + transformVector[6 + screenAreaIndex] * rawTrackingXPos * rawTrackingYPos + transformVector[7 + screenAreaIndex];

        break;

    case 16:

        if (rawTrackingYPos <= calibrationYPos[5]) // first row
        {
            if (rawTrackingXPos <= calibrationXPos[5])
            {
                screenAreaIndex = 0;
            }
            else if (rawTrackingXPos > calibrationXPos[5] && rawTrackingXPos <= calibrationXPos[10])
            {
                screenAreaIndex = 1;
            }
            else
            {
                screenAreaIndex = 2;
            }
        }
        else if (rawTrackingYPos > calibrationYPos[5] && rawTrackingYPos <= calibrationYPos[10]) // second row
        {
            if (rawTrackingXPos <= calibrationXPos[5])
            {
                screenAreaIndex = 3;
            }
            else if (rawTrackingXPos > calibrationXPos[5] && rawTrackingXPos <= calibrationXPos[10])
            {
                screenAreaIndex = 4;
            }
            else
            {
                screenAreaIndex = 5;
            }
        }
        else if (rawTrackingYPos > calibrationYPos[10]) // third row
        {
            if (rawTrackingXPos <= calibrationXPos[5])
            {
                screenAreaIndex = 6;
            }
            else if (rawTrackingXPos > calibrationXPos[5] && rawTrackingXPos <= calibrationXPos[10])
            {
                screenAreaIndex = 7;
            }
            else
            {
                screenAreaIndex = 8;
            }
        }

        screenAreaIndex = screenAreaIndex * 8;

        newTrackingCoordinates[0] = transformVector[0 + screenAreaIndex] * rawTrackingXPos + transformVector[1 + screenAreaIndex] * rawTrackingYPos + transformVector[2 + screenAreaIndex] * rawTrackingXPos * rawTrackingYPos + transformVector[3 + screenAreaIndex];
        newTrackingCoordinates[1] = transformVector[4 + screenAreaIndex] * rawTrackingXPos + transformVector[5 + screenAreaIndex] * rawTrackingYPos + transformVector[6 + screenAreaIndex] * rawTrackingXPos * rawTrackingYPos + transformVector[7 + screenAreaIndex];

        break;

    case 25:

        if (rawTrackingYPos <= calibrationYPos[6]) // first row
        {
            if (rawTrackingXPos <= calibrationXPos[6])
            {
                screenAreaIndex = 0;
            }
            else if (rawTrackingXPos > calibrationXPos[6] && rawTrackingXPos <= calibrationXPos[12])
            {
                screenAreaIndex = 1;
            }
            else if (rawTrackingXPos > calibrationXPos[12] && rawTrackingXPos <= calibrationXPos[18])
            {
                screenAreaIndex = 2;
            }
            else
            {
                screenAreaIndex = 3;
            }
        }
        else if (rawTrackingYPos > calibrationYPos[6] && rawTrackingYPos <= calibrationYPos[12]) // second row
        {
            if (rawTrackingXPos <= calibrationXPos[6])
            {
                screenAreaIndex = 4;
            }
            else if (rawTrackingXPos > calibrationXPos[6] && rawTrackingXPos <= calibrationXPos[12])
            {
                screenAreaIndex = 5;
            }
            else if (rawTrackingXPos > calibrationXPos[12] && rawTrackingXPos <= calibrationXPos[18])
            {
                screenAreaIndex = 6;
            }
            else
            {
                screenAreaIndex = 7;
            }
        }
        else if (rawTrackingYPos > calibrationYPos[12] && rawTrackingYPos <= calibrationYPos[18]) // third row
        {
            if (rawTrackingXPos <= calibrationXPos[6])
            {
                screenAreaIndex = 8;
            }
            else if (rawTrackingXPos > calibrationXPos[6] && rawTrackingXPos <= calibrationXPos[12])
            {
                screenAreaIndex = 9;
            }
            else if (rawTrackingXPos > calibrationXPos[12] && rawTrackingXPos <= calibrationXPos[18])
            {
                screenAreaIndex = 10;
            }
            else
            {
                screenAreaIndex = 11;
            }
        }
        else if (rawTrackingYPos > calibrationYPos[18])  // fourth row
        {
            if (rawTrackingXPos <= calibrationXPos[6])
            {
                screenAreaIndex = 12;
            }
            else if (rawTrackingXPos > calibrationXPos[6] && rawTrackingXPos <= calibrationXPos[12])
            {
                screenAreaIndex = 13;
            }
            else if (rawTrackingXPos > calibrationXPos[12] && rawTrackingXPos <= calibrationXPos[18])
            {
                screenAreaIndex = 14;
            }
            else
            {
                screenAreaIndex = 15;
            }
        }

        screenAreaIndex = screenAreaIndex * 8;

        newTrackingCoordinates[0] = transformVector[0 + screenAreaIndex] * rawTrackingXPos + transformVector[1 + screenAreaIndex] * rawTrackingYPos + transformVector[2 + screenAreaIndex] * rawTrackingXPos * rawTrackingYPos + transformVector[3 + screenAreaIndex];
        newTrackingCoordinates[1] = transformVector[4 + screenAreaIndex] * rawTrackingXPos + transformVector[5 + screenAreaIndex] * rawTrackingYPos + transformVector[6 + screenAreaIndex] * rawTrackingXPos * rawTrackingYPos + transformVector[7 + screenAreaIndex];

        break;
    }

    return newTrackingCoordinates;
}


bool MainWindow::CalibrateData(int trackingDataIndex)
{
    std::vector<double> calibrationLeftXPosAvg(gridSize);
    std::vector<double> calibrationLeftYPosAvg(gridSize);
    std::vector<double> calibrationRghtXPosAvg(gridSize);
    std::vector<double> calibrationRghtYPosAvg(gridSize);

    for (int iCalibrationPoint = 0; iCalibrationPoint < gridSize; iCalibrationPoint++)
    {
        std::vector<double> calibrationLeftXPos;
        std::vector<double> calibrationLeftYPos;
        std::vector<double> calibrationRghtXPos;
        std::vector<double> calibrationRghtYPos;

        {
            // Grab calibration data

            std::stringstream filename;
            filename << dataDirectory
                     << "/" << currentDate
                     << "/" << (NameInputLineEdit->text()).toStdString()
                     << "/calibration/point_" << iCalibrationPoint
                     << "/" << (NameInputLineEdit->text()).toStdString()
                     << "_calibration_data_" << iCalibrationPoint << ".dat";

            std::ifstream data(filename.str().c_str(), std::ios::in);

            if (data.is_open())
            {
                std::string line;
                while (std::getline(data, line))
                {
                    std::istringstream iss(line);
                    double leftXPos, leftYPos, rghtXPos, rghtYPos;
                    bool leftFlag, rghtFlag;

                    if (!(iss >> leftXPos >> leftYPos >> leftFlag >> rghtXPos >> rghtYPos >> rghtFlag))
                    {
                        break; // error
                    }

                    if (leftFlag)
                    {
                        calibrationLeftXPos.push_back(leftXPos);
                        calibrationLeftYPos.push_back(leftYPos);
                    }

                    if (rghtFlag)
                    {
                        calibrationRghtXPos.push_back(rghtXPos);
                        calibrationRghtYPos.push_back(rghtYPos);
                    }
                }
            }
            else
            {
                std::cout << "Calibration data of point " << iCalibrationPoint << " is missing" << std::endl;
                return false;
            }
        }

        calibrationLeftXPosAvg[iCalibrationPoint] = removeOutliers(calibrationLeftXPos, Parameters::outlierFactor);
        calibrationLeftYPosAvg[iCalibrationPoint] = removeOutliers(calibrationLeftYPos, Parameters::outlierFactor);
        calibrationRghtXPosAvg[iCalibrationPoint] = removeOutliers(calibrationRghtXPos, Parameters::outlierFactor);
        calibrationRghtYPosAvg[iCalibrationPoint] = removeOutliers(calibrationRghtYPos, Parameters::outlierFactor);
    }

    // tracking data

    std::vector<double> trackingLeftXPos;
    std::vector<double> trackingLeftYPos;
    std::vector<double> trackingRghtXPos;
    std::vector<double> trackingRghtYPos;
    std::vector<bool> trackingLeftDetections;
    std::vector<bool> trackingRghtDetections;

    // time stamps

    std::vector<double> timeStamps;

    std::vector<double> transformVectorLeft = CalibrationGetTransformVector(calibrationLeftXPosAvg, calibrationLeftYPosAvg);
    std::vector<double> transformVectorRght = CalibrationGetTransformVector(calibrationRghtXPosAvg, calibrationRghtYPosAvg);

    { // grab raw tracking data

        std::stringstream filename;
        filename << dataDirectory
                 << "/" << currentDate
                 << "/" << (NameInputLineEdit->text()).toStdString()
                 << "/experiment/trial_" << trackingDataIndex
                 << "/" << (NameInputLineEdit->text()).toStdString()
                 << "_tracking_data_raw_" << trackingDataIndex << ".dat";

        std::ifstream data(filename.str().c_str(), std::ios::in);

        if (data.is_open())
        {
            std::string line;
            while (std::getline(data, line))
            {
                std::istringstream iss(line);

                double leftXPos, leftYPos, rghtXPos, rghtYPos, timeStamp;
                bool leftFlag, rghtFlag;

                if (!(iss >> leftXPos >> leftYPos >> leftFlag >> rghtXPos >> rghtYPos >> rghtFlag >> timeStamp))
                {
                    break;
                }

                trackingLeftXPos.push_back(leftXPos);
                trackingLeftYPos.push_back(leftYPos);
                trackingRghtXPos.push_back(rghtXPos);
                trackingRghtYPos.push_back(rghtYPos);

                trackingLeftDetections.push_back(leftFlag);
                trackingRghtDetections.push_back(rghtFlag);

                timeStamps.push_back(timeStamp);
            }
        }
        else
        {
            return false;
        }
    }

    // perform geometric transformation and save data

    std::stringstream filename;
    filename << dataDirectory
             << "/" << currentDate
             << "/" << (NameInputLineEdit->text()).toStdString()
             << "/experiment/trial_" << trackingDataIndex
             << "/" << (NameInputLineEdit->text()).toStdString()
             << "_tracking_data_" << trackingDataIndex << ".dat";

    std::ofstream data;
    data.open(filename.str(), std::ios::out | std::ios::ate);

    for (int i = 0, lineCount = trackingLeftXPos.size(); i < lineCount; i++)
    {
        std::vector<double> leftCoords = CalibrationGeometricTransform(transformVectorLeft, calibrationLeftXPosAvg, calibrationLeftYPosAvg, trackingLeftXPos[i], trackingLeftYPos[i]);
        std::vector<double> rghtCoords = CalibrationGeometricTransform(transformVectorRght, calibrationRghtXPosAvg, calibrationRghtYPosAvg, trackingRghtXPos[i], trackingRghtYPos[i]);

        data << leftCoords[0] << " " << leftCoords[1] << " " << trackingLeftDetections[i] << " " << rghtCoords[0] << " " << rghtCoords[1] << " " << trackingRghtDetections[i] << " " << timeStamps[i] << std::endl;
    }

    data.close();

    return true;
}

void MainWindow::CalibrateDataAll()
{
    int iTrial = 0; while(CalibrateData(iTrial)) { iTrial++; }
}

void MainWindow::CalibrationSave()
{
    // save info in .ini file

    std::stringstream saveFileNameSS;
    saveFileNameSS << dataDirectory
                   << "/" << currentDate
                   << "/" << (NameInputLineEdit->text()).toStdString()
                   << "/calibration/point_" << calibrationIndex;

    editDataDirectory = QString::fromStdString(saveFileNameSS.str());

    saveFileNameSS << "/info.ini";

    editDataIndex = calibrationIndex;
    editDataType  = 0;

    editImageTotal = calibrationFrames;

    saveSettings(QString::fromStdString(saveFileNameSS.str()));
}
