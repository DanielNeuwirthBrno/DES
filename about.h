/*******************************************************************************
 Copyright 2017 Daniel Neuwirth
 This program is distributed under the terms of the GNU General Public License.
*******************************************************************************/

#ifndef ABOUT_H
#define ABOUT_H

// This header file defines so-called About window which displays basic program
// information regarding its version, author, programming environment etc.


#include "ui/ui_about.h"

#include <QDialog>
#include <QWidget>

class AboutWindow: public QDialog {

    Q_OBJECT

    public:
        explicit AboutWindow(QWidget * parent = 0);
        ~AboutWindow();

    private:
        Ui_AboutWindow * ui;

    private slots:
};

#endif // ABOUT_H
