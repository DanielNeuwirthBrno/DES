/*******************************************************************************
 Copyright 2017 Daniel Neuwirth
 This program is distributed under the terms of the GNU General Public License.
*******************************************************************************/

#ifndef STEPS_H
#define STEPS_H

// This header file defines behaviour (and reactions on user's input) of
// window which displays individual steps made during computation of equation.


#include <token.h>
#include "ui/ui_steps.h"

#include <utility>
#include <QChar>
#include <QDialog>
#include <QString>
#include <QWidget>

class StepsWindow: public QDialog {

    Q_OBJECT

    public:
        explicit StepsWindow(const QChar, const QString &,
            const std::pair<bool, Matrix<double_t>> &, QWidget * parent = 0);
        ~StepsWindow();

    private:
        QChar diffVariable;
        Matrix<double_t> fundamentalMatrix;
        Ui_StepsWindow * ui;

    private slots:
        int fundamentalMatrix_clicked();
};

#endif // STEPS_H
