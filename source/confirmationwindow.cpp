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

#include "headers/confirmationwindow.h"

ConfirmationWindow::ConfirmationWindow(QString text, QWidget *parent) : QDialog(parent)
{
    RETURN_VALUE = false;

    QLabel *InfoTextBox = new QLabel;
    InfoTextBox->setText(text);

    QPushButton *ApproveButton = new QPushButton("&OK");
    QObject::connect(ApproveButton, SIGNAL(clicked()), this, SLOT(setApprove()));

    QPushButton *DisapproveButton = new QPushButton("&Cancel");
    QObject::connect(DisapproveButton, SIGNAL(clicked()), this, SLOT(setDisapprove()));

    QHBoxLayout *ButtonLayout = new QHBoxLayout;
    ButtonLayout->addStretch();
    ButtonLayout->addWidget(ApproveButton);
    ButtonLayout->addStretch();
    ButtonLayout->addWidget(DisapproveButton);
    ButtonLayout->addStretch();

    QVBoxLayout *MainLayout = new QVBoxLayout;
    MainLayout->addWidget(InfoTextBox);
    MainLayout->addLayout(ButtonLayout);
    setLayout(MainLayout);
}

void ConfirmationWindow::setApprove()
{
    accept();
    close();
}

void ConfirmationWindow::setDisapprove()
{
    reject();
    close();
}

bool ConfirmationWindow::getStatus()
{
    return RETURN_VALUE;
}

