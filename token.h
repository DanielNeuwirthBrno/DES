/*******************************************************************************
 Copyright 2016-17 Daniel Neuwirth
 This program is distributed under the terms of the GNU General Public License.
*******************************************************************************/

#ifndef TOKEN_H
#define TOKEN_H

// This header file contains definitions of several classes.
// Token class is used for storing individual tokens which is submitted
// differential equation broken into.
// Roots class is used for storing roots of characteristic equation and
// its textual description.
// InitConditions class stores entered initial conditions.


#include "compute.h"

#include <QChar>
#include <QMatrix4x4>
#include <QString>
#include <QVector>

class Token {

    public:
        Token();
        Token(const double_t, const QChar, const uint8_t, const int16_t);
        ~Token() {}

        inline double_t getCoefficient() const
            { return coefficient; }
        inline complex getComplexCoef() const
            { return compute::toComplex(coefficient); }
        inline const QChar getVariable() const
            { return variable; }
        inline bool isDerived() const
            { return (static_cast<bool>(order)); }
        inline int8_t getPower() const
            { return power; }
        const QString getToken() const;

    private:
        double_t coefficient;
        QChar variable;
        uint8_t order;
        int8_t power;
};

class Derivatives {

    public:
        Derivatives() = default;
        Derivatives(const QString, const double_t);
        ~Derivatives() {}

        inline const QString getDerivative() const
            { return derivative; }
        inline double_t evalDerivative() const
            { return numValueAtInit; }

    private:
        QString derivative;
        double_t numValueAtInit;
};

class Roots {

    public:
        // enumeration for solution output
        enum returnPartOfSolution { COMPLETE_SOLUTION, ONLY_PART_A, ONLY_PART_B,
                                    TEST_REAL_ROOT_FIRST };

        Roots() = default;
        Roots(const complex, const QString &, const QString &, const double_t = 0,
              const QVector<Derivatives> & = QVector<Derivatives>());
        ~Roots() {}

        inline double_t getRealPartOfRoot() const
            { return root.real(); }
        inline double_t getImagPartOfRoot() const
            { return root.imag(); }
        const QString getSolution(returnPartOfSolution = COMPLETE_SOLUTION) const;
        inline double_t getNumValue() const
            { return rootNumValueAtInit; }
        inline const QString getDerivative(const uint8_t row) const
            { return derivatives.at(row).getDerivative(); }
        inline double_t evalDerivativeAtInit(const uint8_t row) const
            { return derivatives.at(row).evalDerivative(); }
        inline void addDerivative(const Derivatives & inDerivatives)
            { derivatives.push_back(inDerivatives); return; }

    private:
        complex root;
        QString solutionPartA;
        QString solutionPartB;
        double_t rootNumValueAtInit;
        QVector<Derivatives> derivatives;
};

class InitConditions {

    public:
        InitConditions() {}
        ~InitConditions() {}

        inline double_t getXValue() const { return xValue; }
        inline double_t getFunctionValue() const { return functionValue; }
        inline void setValues(const double_t inXValue, const double_t inYValue)
            { xValue = inXValue; functionValue = inYValue; return; }

        const QString display(const uint8_t) const;

    private:
        double_t xValue;
        double_t functionValue;
};

template <class T>
class Matrix {

    public:
        Matrix() = default;
        Matrix(const uint8_t inOrder) { setSize(inOrder); }
        ~Matrix() {}

        inline T getValue(const uint8_t row, const uint8_t column) const
            { return matrix[row][column]; }
        inline T & at(const uint8_t row, const uint8_t column)
            { return matrix[row][column]; }

        inline uint8_t size() const { return matrix.size(); }
        inline void empty() { matrix.resize(0); return; }
        void setSize(const uint8_t inOrder) {
            QVector<QVector<T>> matrixWithSize(inOrder, QVector<T>(inOrder));
            matrix = matrixWithSize;
        }

        // compute determinant of 2x2 and 3x3 matrix
        const T determinant3() const;
        // compute determinant of 4x4 matrix
        const T determinant4(const float * const inValues) const
          { const QMatrix4x4 qMatrix(inValues); return qMatrix.determinant(); }
\
        void operator=(QVector<QVector<T>> & inMatrix)
            { matrix = inMatrix; }

    private:
        QVector<QVector<T>> matrix;
};

#endif // TOKEN_H
