/*******************************************************************************
 Copyright 2017 Daniel Neuwirth
 This program is distributed under the terms of the GNU General Public License.
*******************************************************************************/

#ifndef MATRIX_H
#define MATRIX_H

// This header file defines window used for displaying fundamental matrix
// ie. square matrix whose determinant (wronskian) serves the purpose of:
// a] determining if computed roots form a fundamental set of solutions
// b] computation of numerical values of constants


#include "token.h"
#include "ui/ui_matrix.h"

#include <QDialog>
#include <QWidget>

class MatrixWindow: public QDialog {

    Q_OBJECT

    public:
        explicit MatrixWindow(const QChar, const Matrix<double_t> &, QWidget * = 0);
        ~MatrixWindow();

    private:
        Ui_MatrixWindow * ui;

    private slots:
};

#endif // MATRIX_H
