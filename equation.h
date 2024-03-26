/*******************************************************************************
 Copyright 2016-17 Daniel Neuwirth
 This program is distributed under the terms of the GNU General Public License.
*******************************************************************************/

#ifndef EQUATION_H
#define EQUATION_H

// This header file defines Equation class, this application's core part
// containing definitions of variables and corresponding functions needed
// for successful solution of a given differential equation.


#include "general.h"
#include "token.h"

#include <bitset>
#include <utility>
#include <QChar>
#include <QString>
#include <QVector>

class Equation {

    public:
        // enumeration for goniometric functions
        enum gonFunction { UNKNOWN, SIN, COS };
        // enumeration for equation types
        enum typeOfEquation { LINEAR, HOMOGENEOUS, ORDINARY, CONST_COEFF };

        Equation(const bool, const QVector<InitConditions>, const QChar,
                 const uint8_t, const QString &, const QString &);
        ~Equation() {}

        inline QString getEquationLeftSide() const
            { return equationLeftSide; }
        inline QString getEquationRightSide() const
            { return equationRightSide; }
        inline QString getEquationComplete() const
            { return (equationLeftSide + " = " + equationRightSide); }
        inline QChar getDiffVariable() const
            { return diffVariable; }
        inline bool applyInitialConditions() const
            { return initialConditions.first; }
        inline const QVector<InitConditions> & getInitialConditions() const
            { return initialConditions.second; }
        inline QString getFinalResult() const
            { return finalResult; }
        inline const QVector<QString> & getSolution() const
            { return completeSolution; }
        inline const QVector<Token> & getTokens() const
            { return tokens; }

        // store information about particular solution step
        inline void saveStep(const QString & inText)
            { completeSolution.push_back(inText); return; }

        // determine if the equation is linear or not
        bool isLinear() const;
        // determine if the equation is homogeneous or not
        bool isHomogeneous() const;
        // determine if the equation is ordinary or not
        bool isOrdinary() const;

        // inject delimiters into the equation to simplify tokenization process
        QString tokenDelimiters(QString, const int16_t = -1) const;
        // move certain elements from the left to the right side of the equation
        bool moveToRightSide();
        // parse equation ~ break it into tokens
        bool parseEquation(const QString &);

        // express order of the equation as a textual representation
        const QString convertToOrdinalNumber(const uint8_t) const;
        // insert appropriate symbols/subscripts for constants
        inline const QString constantsInResult(const uint8_t inNumber) const
            { return ("c<sub>" + QString::number(inNumber) + "</sub>"); }
        // don't display elements raised to the 0th power (= 1)
        inline const QString notZeroPowerElems(const uint8_t inNumber) const
            { return ((inNumber) ? ((inNumber != 1) ?
                "x<sup>" + QString::number(inNumber) + "</sup>" : "x") : ""); }

        // build string representation of derivative
        const QString buildDerivative(const QString &, const double_t,
            const uint8_t, const uint8_t = 0, QString = "");
        // evaluate numerical values of derivatives at initial point
        double_t evalDerivativeRealPart(
            const double_t, const uint8_t, const double_t) const;
        double_t evalDerivativeImagPart(
            const double_t, const uint8_t, double_t (*)(double_t)) const;
        // determine type of the equation (ordinary, linear, etc.)
        const QString determineType(bool &);
        // calculate general solution
        bool generalSolution(std::pair<bool, Matrix<double_t>> &);
        // change (description of) goniometric functions during differentation
        gonFunction changeGonFunction(QString &);
        // (symbolically) express all derivatives up to order-1
        bool computeDerivatives(const uint8_t);
        // calculate wronskian - to determine linear independence of roots
        bool computeWronskian(double_t &, Matrix<double_t> &);
        // infer exact numerical values of constants
        bool determineValuesOfConstants(const double_t, QVector<double_t> &);
        // build solution according to principle of superposition
        const QString principleOfSuperposition();

    private:
        // could be declared as initConditionsStruct (see using declarative)
        // left in this form for clear understanding what's under the hood
        std::pair<bool,QVector<InitConditions>> initialConditions;
        QString finalResult;
        QVector<QString> completeSolution;
        CharEq characteristicEquation;
        QVector<Token> tokens;

        QChar diffVariable;
        uint8_t order;
        std::bitset<4> type;
        QString equationLeftSide;
        QString equationRightSide;
};

#endif // EQUATION_H
