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

#include "startupwindow.h"

StartUpWindow::StartUpWindow(QWidget *parent) : QDialog(parent)
{
    SAVE_SUCCESS = false;

    QLabel *InfoTextBox = new QLabel;
    InfoTextBox->setText("First time set-up. Please specify where data should be saved.");

    QLabel *DirectoryTextBox = new QLabel;
    DirectoryTextBox->setText("<b>Save directory</b>");

    DirectoryLineEdit = new QLineEdit;
    DirectoryLineEdit->setAlignment(Qt::AlignLeft);

    QPushButton *DirectoryButton = new QPushButton("&Browse");
    QObject::connect(DirectoryButton, SIGNAL(clicked()), this, SLOT(SelectDirectory()));

    QPushButton *SaveButton = new QPushButton("&OK");
    QObject::connect(SaveButton, SIGNAL(clicked()), this, SLOT(SaveDirectory()));

    QHBoxLayout *DirectoryLayout = new QHBoxLayout;
    DirectoryLayout->addWidget(DirectoryTextBox);
    DirectoryLayout->addWidget(DirectoryLineEdit);
    DirectoryLayout->addWidget(DirectoryButton);

    QVBoxLayout *MainLayout = new QVBoxLayout;
    MainLayout->addWidget(InfoTextBox);
    MainLayout->addLayout(DirectoryLayout);
    MainLayout->addWidget(SaveButton, 0, Qt::AlignRight);
    setLayout(MainLayout);

    DirectoryLineEdit->text();
}

StartUpWindow::~StartUpWindow()
{

}

void StartUpWindow::SelectDirectory()
{
    directory = QFileDialog::getExistingDirectory(this, tr("Open Directory"), "/home", QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    DirectoryLineEdit->setText(directory);
}

void StartUpWindow::SaveDirectory()
{
    QSettings settings("config.ini", QSettings::IniFormat);
    settings.setValue("dataDirectory", (DirectoryLineEdit->text()));

    SAVE_SUCCESS = true;
    close();
}

bool StartUpWindow::getStatus()
{
    return SAVE_SUCCESS;
}
