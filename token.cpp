#include "token.h"


Token::Token(): coefficient(0.0), variable('\0'), order(0), power(1) {}

Token::Token(const double_t inCoef, const QChar inVariable,
             const uint8_t inOrder, const int16_t inPower):
    coefficient(inCoef), variable(inVariable), order(inOrder), power(inPower) {}

const QString Token::getToken() const {

    const QString token = QString::number(coefficient, 'f') + variable + "{" +
                          QString::number(order) + "}^" + QString::number(power);
    return (token);
}

Derivatives::Derivatives(const QString inDerivative,
                         const double_t inNumValueAtInit = 0):
    derivative(inDerivative), numValueAtInit(inNumValueAtInit) {}

Roots::Roots(const complex inNumber, const QString & inPartA,
             const QString & inPartB, const double_t inNumValueOfRoot,
             const QVector<Derivatives> & inDerivatives):
             root(inNumber), solutionPartA(inPartA), solutionPartB(inPartB),
             rootNumValueAtInit(inNumValueOfRoot), derivatives(inDerivatives) {}

const QString Roots::getSolution(returnPartOfSolution inWhichPart) const {

    switch (inWhichPart) {

        case ONLY_PART_A: return (solutionPartA);
                          break;
        case ONLY_PART_B: return (solutionPartB);
                          break;
        case TEST_REAL_ROOT_FIRST:
                          if (root.real() == 0)
                              return (solutionPartB);
        case COMPLETE_SOLUTION:
        default: return (solutionPartA+solutionPartB);
    }
}

const QString InitConditions::display(const uint8_t inOrder) const {

    QString displayCondition = format::expressDerivativesAsSymbols(inOrder) +
        "(" + format::numberForOutput(getXValue()) +
        ") = " + format::numberForOutput(getFunctionValue());

    return (displayCondition);
}

template <class T>
const T Matrix<T>::determinant3() const {

    T determinant = 0;

    // run only once for 2x2 matrix (determinant = ad - bc)
    // run three times for 3x3 matrix (determinant computed using Sarrus' rule)
    for (uint8_t i = 0; i < size()-(size()%3)/2; ++i) {

        T partToAdd = 1, partToSubtract = 1;
        for (int j = 0; j < size(); ++j) {

            partToAdd *= matrix[j % size()][(j + i) % size()];
            partToSubtract *= matrix[j % size()][(size()-j-1 + i) % size()];
        }
        determinant += partToAdd - partToSubtract;
    }
    return (determinant);
}

template class Matrix<double_t>;
template class Matrix<int>;
