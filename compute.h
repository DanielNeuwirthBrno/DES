/*******************************************************************************
 Copyright 2016-17 Daniel Neuwirth
 This program is distributed under the terms of the GNU General Public License.
*******************************************************************************/

#ifndef COMPUTE_H
#define COMPUTE_H

// This header file contains declarations of auxiliary mathematical functions
// in namespace compute and other (helper, formatting) functions in namespaces
// determine and format.


#include <complex>
#include <QLineEdit>
#include <QString>

using complex = std::complex<double_t>;

namespace compute {

    const uint8_t compareFloatsToDecPlaces = 6;
    const complex rootTwoOfCubic(-0.5,0.5*sqrt(3));
    const complex rootThreeOfCubic(-0.5,-0.5*sqrt(3));

    // compute cubic root of complex number
    complex cbrt(const complex inNumber);

    // round floating-point numbers
    complex round(const complex, const uint8_t = compareFloatsToDecPlaces);
    double_t round(const double_t, const uint8_t = compareFloatsToDecPlaces);

    // convert real numbers to complex numbers
    template<typename T>
    inline const complex toComplex(const T inNumber)
        { return (static_cast<complex>(inNumber)); }
}

namespace determine {

    // test if all chars in a string are the same
    inline bool allCharsTheSame(const QString & inString) {

        return (inString.size() == inString.count(inString[0]));
    }
}

namespace format {

    const QString expressDerivativesAsSymbols(const uint8_t);

    // load html file (description of particular element) from resources
    const QString loadRes(const QString &);
    // format output of floating-point numbers (represented as strings)
    const QString numberForOutput(const double_t, const bool = false);
    // set appropriate width of QLineEdit element (according to font metrics)
    void setElementWidth(QLineEdit * const);
}

#endif // COMPUTE_H
