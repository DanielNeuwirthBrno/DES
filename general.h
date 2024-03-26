/*******************************************************************************
 Copyright 2016-17 Daniel Neuwirth
 This program is distributed under the terms of the GNU General Public License.
*******************************************************************************/

#ifndef GENERAL_H
#define GENERAL_H

// This header file defines CharEq class used for computing and storing roots
// of characteristic equation which is used as a main building block during
// general solution processing.


#include "compute.h"
#include "token.h"

#include <utility>
#include <QString>
#include <QVector>

using initConditionsStruct = std::pair<bool, QVector<InitConditions>>;

// sort using a custom function objects
struct sortRealRootsOverComplex {

    bool operator()(const Roots inFirst, const Roots inSecond)
        { return (std::abs(inFirst.getImagPartOfRoot())) <
                  std::abs(inSecond.getImagPartOfRoot()); }
};

struct sortRealRootsInAscOrder {

    bool operator()(const Roots inFirst, const Roots inSecond)
       { return (inFirst.getRealPartOfRoot() < inSecond.getRealPartOfRoot()); }
};

class CharEq {

    public:
        enum equationDegree { QUADRATIC = 1, BIQUADRATIC = 2 };

        CharEq(const uint8_t noOfElements);
        ~CharEq() {}

        inline const QVector<Roots> & getRoots() const
            { return rootsOfCharEquation; }
        // for sorting purposes and storing of derivatives
        inline QVector<Roots> & getRootsModifiable()
            { return rootsOfCharEquation; }

        uint8_t noOfRealRoots() const
            { return (rootsOfCharEquation.size()-noOfComplexRoots()); }
        uint8_t noOfComplexRoots() const;
        // is discriminant greater than, less than or equal to zero?
        const QString stateOfDiscriminant(const double_t) const;
        // express solution in text form
        const Roots solution(const complex, const initConditionsStruct &);
        // text description of roots' types used as part of complete solution
        const QString rootsDescription();
        // convert complex number to real and back
        inline void cutOffImagPart(complex & inNumber)
            { inNumber = static_cast<complex>(inNumber.real()); return; }

        // compute roots of quadratic equation
        const QString rootsOfQuadraticEquation(const QVector<Token> &,
            const initConditionsStruct &, const equationDegree = QUADRATIC);
        // compute roots of cubic equation
        const QString rootsOfCubicEquation(const QVector<Token> &,
            const initConditionsStruct &);
        // compute roots of quartic equation
        const QString rootsOfQuarticEquation(const QVector<Token> &,
            const initConditionsStruct &, bool &);
        // compute roots of n-th degree polynomial equation
        const QString rootsOfPolynomialEquation(const QVector<Token> &);

    private:
        QVector<Roots> rootsOfCharEquation;
};

#endif // GENERAL
