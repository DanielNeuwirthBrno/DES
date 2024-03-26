/*******************************************************************************
 Copyright 2017 Daniel Neuwirth
 This program is distributed under the terms of the GNU General Public License.
*******************************************************************************/

#ifndef CONDITIONS_H
#define CONDITIONS_H

// This header file defines window designated for insertion/editing of initial
// conditions; geometry of this window is adjusted according to current
// number of initial conditions (depends on the highest order of the VWADWRT)


#include "token.h"
#include "ui/ui_conditions.h"

#include <QChar>
#include <QDialog>
#include <QVector>
#include <QWidget>

class ConditionsWindow: public QDialog {

    Q_OBJECT

    public:
        explicit ConditionsWindow(const QChar, const uint8_t,
                                  QVector<InitConditions> &, QWidget * = 0);
        ~ConditionsWindow();

    private:
        QVector<InitConditions> & conditions;
        Ui_ConditionsWindow * ui;

    private slots:
        void saveConditions();
};

#endif // CONDITIONS_H
