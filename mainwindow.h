/*******************************************************************************
 Copyright 2016-17 Daniel Neuwirth
 This program is distributed under the terms of the GNU General Public License.
*******************************************************************************/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

// This header file defines MainWindow's behaviour, how it reacts on user's
// input and preprocesses (checks) some elements of inserted diff.equation.


#include "equation.h"
#include "token.h"
#include "ui/ui_mainwindow.h"

#include <memory>
#include <utility>
#include <QApplication>
#include <QChar>
#include <QMainWindow>
#include <QString>
#include <QVector>

class MainWindow: public QMainWindow {

    Q_OBJECT

    public:
        // enumeration for (catchable) syntax errors
        enum syntaxError { UNSPECIFIED = 0, OK = 1, UNRECOGNIZED_CHAR = 10,
                           NOT_MATCHING_BRACES = 20, DIFF_SIGN_MISPLACED = 30,
                           NO_VARIABLE_TO_DIFF = 40};

        explicit MainWindow(QWidget * parent = 0);
        ~MainWindow();

        // put through only allowed characters
        inline bool isCharAllowed(const QChar inChar) const
            { if (inChar >= '1' && inChar <= '9')
                return (true); else return (false); }
        // display warning about type of detected error
        const QString errorDescription(const syntaxError = UNSPECIFIED,
                                       const QChar = '\0') const;
        // first-wave check for syntax errors
        syntaxError checkForSyntaxErrors
            (const QString &, const QString &, const QChar = '\0') const;

    private:
        std::unique_ptr<Equation> diffEquationToSolve;
        QVector<InitConditions> setOfInitConditions;
        std::pair<bool, Matrix<double_t>> fundamentalMatrix;
        int cursorPosition;
        const static char multipleDiffVars = '_';
        const static QString findDiffVarRegex;
        const static QString allowedSpecialCharacters;
        const static QString containsErrorBackground;
        const static QString noSolution;
        Ui_MainWindow * ui;

    private slots:
        // extract VWADWRT from given expression
        bool extractDiffVariable(QString) const;
        void preprocessing_l(const QString &) const;
        void preprocessing_r(const QString &) const;
        void insertDiffAtCaret() const;
        void enableSubmitButton() const;
        void enableAddConditionsButton() const;
        void checkIfOrderEqualsNoOfConditions();
        int addInitialConditions();
        void setInitialConditions();
        void applyConditionsIcon(const bool) const;
        void enableShowStepsButton() const;
        void submitEquation_clicked();
        void clearEquation_clicked();
        int showSteps_clicked();

        // actions
        int action_About();
        void action_Exit() { qApp->quit(); }
};

#endif // MAINWINDOW_H
