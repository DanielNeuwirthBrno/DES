#include "compute.h"
#include "general.h"

#include <cmath>


CharEq::CharEq(const uint8_t noOfElements) {

    rootsOfCharEquation.reserve(noOfElements);
}

uint8_t CharEq::noOfComplexRoots() const {

    uint8_t noOfRoots = 0;
    for (auto it : rootsOfCharEquation)
        if (it.getImagPartOfRoot() != 0) ++noOfRoots;

    return (noOfRoots);
}

const QString CharEq::stateOfDiscriminant(const double_t inDiscriminant) const {

    QString disDescription = " = 0";
    if (inDiscriminant)
        disDescription[1] = (inDiscriminant > 0) ? '>' : '<';
    return (disDescription);
}

const Roots CharEq::solution(
    const complex root, const initConditionsStruct & inConditions) {

    // format output (truncate trailing zeros)
    const QString rootFormatted = format::numberForOutput(root.real(), true);

    // root part A (real part)
    const QString rootPartA = QStringLiteral("e<sup>") + rootFormatted +
                              QStringLiteral("x</sup>");

    // if initial conditions are supplied then compute root's numerical value
    const double_t conditionValue =
        (inConditions.first) ? inConditions.second.at(0).getXValue() : 0;
    double_t rootNumValue = std::exp(root.real()*conditionValue);

    // root part B (imaginary part)
    QString rootPartB = "";

    if (compute::round(root.imag()) != 0) {

        QString gonFunction;
        if (root.imag() > 0) {
            gonFunction = QStringLiteral("cos(");
            rootNumValue *= std::cos(std::abs(root.imag())*conditionValue);
        }
        else {
            gonFunction = QStringLiteral("sin(");
            rootNumValue *= std::sin(std::abs(root.imag())*conditionValue);
        }
        rootPartB =
            gonFunction + format::numberForOutput(std::abs(root.imag()), true) + "x)";
    }

    // check for repeated roots => numValue must be multiplied accordingly
    if (inConditions.first)
        for (auto it: rootsOfCharEquation)
            if (it.getSolution() == rootPartA+rootPartB)
                rootNumValue *= conditionValue;

    Roots newRoot(root, rootPartA, rootPartB, rootNumValue);
    return newRoot;
}

const QString CharEq::rootsDescription() {

    QString roots = "type of roots: ";
    const uint8_t noOfRealRoots = CharEq::noOfRealRoots();
    uint8_t noOfRealDistinctRoots = noOfRealRoots;

    // row of vector represents number of repetitions (double, triple etc.),
    // contents represent number of roots with this amount of repetitions
    QVector<uint8_t> repeatedRoots(noOfRealRoots, 0);
    uint8_t repeated = 1;
    auto endOfVector = rootsOfCharEquation.begin()+noOfRealRoots;

    // cycle only through real roots
    for (auto it = endOfVector-noOfRealRoots+1; it < endOfVector; ++it) {

        // repeated ...
        if ((it-1)->getRealPartOfRoot() == it->getRealPartOfRoot()) {

            ++repeated;
            --noOfRealDistinctRoots;
        }

        // not repeated ...
        if (((it-1)->getRealPartOfRoot() != it->getRealPartOfRoot() ||
            (it+1) == endOfVector) && (repeated != 1)) {

            ++(repeatedRoots[repeated-1]);
            repeated = 1;
        }
    }

    if (noOfRealRoots > 0) {

        roots += QString::number(noOfRealRoots) + " real roots";

        // all real roots are distinct
        if (noOfRealRoots == noOfRealDistinctRoots) {

            roots.insert(roots.lastIndexOf("roots"),"distinct ");
            if (noOfRealRoots == 1)
                roots.chop(1);

        }
        // some of the roots are repetitive
        else {

            QString multiples[] = { "double", "triple" };
            roots += " (";

            for (uint8_t i = 0; i < repeatedRoots.size(); ++i) {

                if (repeatedRoots.at(i) > 0) {

                    roots += QString::number(repeatedRoots.at(i))+ "x ";
                    if (i <= 2) // use textual description
                        roots += multiples[i-1] + " root, ";
                    else // use number description
                        roots += QString::number(i+1) + "-times repeated root, ";
                }
            }
            roots.chop(2);
            roots.append(')');
        }
    }

    if (noOfComplexRoots() > 0) {

        if (noOfRealRoots > 0)
            roots += ", ";
        roots += QString::number(noOfComplexRoots()) + " complex roots";
    }

    return (roots);
}

const QString CharEq::rootsOfQuadraticEquation(const QVector<Token> & inTokens,
    const initConditionsStruct & inConditions, const equationDegree inDegree) {

    const complex a = inTokens.at(0).getComplexCoef(),
                  b = inTokens.at(1*inDegree).getComplexCoef(),
                  c = inTokens.at(2*inDegree).getComplexCoef();

    // compute discriminant
    const complex discriminant = b*b - compute::toComplex(4)*a*c;

    // compute roots
    const complex denominator = (a*compute::toComplex(2));
    for (uint8_t i = 0; i < inDegree; ++i) {

        complex root1 = (-b - std::sqrt(discriminant)) / denominator;
        complex root2 = (-b + std::sqrt(discriminant)) / denominator;

        if (inDegree == BIQUADRATIC) {
            root1 = compute::toComplex(1-i*2/1)*sqrt(root1);
            root2 = compute::toComplex(1-i*2/1)*sqrt(root2);
        }
        rootsOfCharEquation.push_back(solution(root1,inConditions));
        rootsOfCharEquation.push_back(solution(root2,inConditions));
    }

    if (inDegree == BIQUADRATIC)
        return (QStringLiteral("quartic (biquadratic) equation"));

    const QString equationType =
        QStringLiteral("quadratic equation (discriminant ") +
        stateOfDiscriminant(discriminant.real()) + ')';

    return (equationType);
}

const QString CharEq::rootsOfCubicEquation(const QVector<Token> & inTokens,
    const initConditionsStruct & inConditions) {

    const complex a = inTokens.at(0).getComplexCoef(),
                  b = inTokens.at(1).getComplexCoef(),
                  c = inTokens.at(2).getComplexCoef(),
                  d = inTokens.at(3).getComplexCoef();

    // compute discriminant
    const complex temp1 = compute::toComplex(27)*a*a;
    const complex discriminant =
        compute::toComplex(18)*a*b*c*d - compute::toComplex(4)*std::pow(b,3)*d +
        b*b*c*c - compute::toComplex(4)*a*std::pow(c,3) - temp1*d*d;

    const complex disD0 = b*b - compute::toComplex(3)*a*c;
    const complex disD1 = compute::toComplex(2)*std::pow(b,3) -
                          compute::toComplex(9)*a*b*c + temp1*d;

    // compute roots
    if (compute::round(discriminant.real()) == 0) { // discriminant is zero

        if (compute::round(disD0.real()) == 0) { // one triple real root

            const complex root = -b / compute::toComplex(3)*a;
            for (uint8_t i=0; i<3; ++i)
                rootsOfCharEquation.push_back(solution(root,inConditions));
        }
        else { // one double real root + one distinct real root

            const complex temp2 = compute::toComplex(9)*a*d;
            complex root1 = (temp2 - b*c) / (compute::toComplex(2)*disD0);
            rootsOfCharEquation.push_back(solution(root1,inConditions));
            rootsOfCharEquation.push_back(solution(root1,inConditions));
            complex root2 = (compute::toComplex(4)*a*b*c -
                             temp2*a - std::pow(b,3)) / (a*disD0);
            rootsOfCharEquation.push_back(solution(root2,inConditions));
        }
    }
    else { // discriminant is non-zero

        complex C;
        bool tripleRoot = false;

        if (discriminant.real() < 0) { // one distinct real root +
                                       // two complex conjugate roots

            // use standard cbrt function => computes cbrt of real numbers

            // firstly calculate only part of C expression (that is under sqrt)
            double_t sqrtPart = std::sqrt(-temp1.real()*discriminant.real());
            // round floating point numbers before comparing them
            if (compute::round(disD1.real()) == compute::round(sqrtPart))
                // expressions in the numerator part of the fraction under cbrt
                // ie. discriminant D1 and sqrtPart must not cancel each other
                C = compute::toComplex(std::cbrt(((disD1.real()+sqrtPart)/2)));
            else
                C = compute::toComplex(std::cbrt(((disD1.real()-sqrtPart)/2)));
        }
        else { // three distinct real roots

            // use custom-made cbrt function => computes cbrt of complex numbers
            C = compute::cbrt(((disD1 + std::sqrt(-temp1*discriminant)) /
                compute::toComplex(2)));
            // although computed root is a real number it could (in theory)
            // have a non-zero imag-part because intermediate result is a
            // complex number and due to the way floating-point numbers are
            // stored in memory there could have emerged imag-part residuals
            // during performing mathematical operations => result must be
            // normalized ie. its imag-part must be zeroed out via
            // cutOffImagPart() function
            tripleRoot = true;
        }

        complex root1 = (compute::toComplex(-1.0/3)/a) * (b+C+disD0/C);
        if (tripleRoot) cutOffImagPart(root1);
        rootsOfCharEquation.push_back(solution(root1,inConditions));

        complex Cx = compute::toComplex(C)*compute::rootTwoOfCubic;
        complex root2 = (-compute::toComplex(1.0/3)/a)*(b+Cx+disD0/Cx);
        if (tripleRoot) cutOffImagPart(root2);
        rootsOfCharEquation.push_back(solution(root2,inConditions));

        Cx = compute::toComplex(C)*compute::rootThreeOfCubic;
        complex root3 = (-compute::toComplex(1.0/3)/a)*(b+Cx+disD0/Cx);
        if (tripleRoot) cutOffImagPart(root3);
        rootsOfCharEquation.push_back(solution(root3,inConditions));
    }

    QString EquationType = QStringLiteral("cubic equation (discriminant ") +
                           stateOfDiscriminant(discriminant.real()) + ')';
    return (EquationType);
}

const QString CharEq::rootsOfQuarticEquation(const QVector<Token> & inTokens,
    const initConditionsStruct & inConditions, bool & canBeSolved) {

    const complex a = inTokens.at(0).getComplexCoef(),
                  b = inTokens.at(1).getComplexCoef(),
                  c = inTokens.at(2).getComplexCoef(),
                  d = inTokens.at(3).getComplexCoef(),
                  e = inTokens.at(4).getComplexCoef();

    if (compute::round(b.real()) == 0 && compute::round(d.real()) == 0) {

        // biquadratic
        const QString equationType =
            rootsOfQuadraticEquation(inTokens, inConditions, BIQUADRATIC);
        return equationType;
    }

    // compute discriminant
    const complex temp1 = a*a*e*e, temp2 = b*b*d*d;
    const complex C4 = compute::toComplex(4);
    const complex discriminant =
        compute::toComplex(256)*temp1*a*e-compute::toComplex(192)*temp1*b*d-
        compute::toComplex(128)*temp1*c*c+compute::toComplex(144)*a*a*c*d*d*e-
        compute::toComplex(27)*a*a*pow(d,4)+compute::toComplex(144)*a*b*b*c*e*e-
        compute::toComplex(6)*a*temp2*e-compute::toComplex(80)*a*b*c*c*d*e+
    compute::toComplex(18)*a*b*c*pow(d,3)+C4*C4*a*pow(c,4)*e-C4*a*pow(c,3)*d*d-
    compute::toComplex(27)*pow(b,4)*e*e+compute::toComplex(18)*pow(b,3)*c*d*e-
    C4*temp2*b*d-C4*b*b*pow(c,3)*e+temp2*c*c;

    const complex disDelta0 =
        c*c - compute::toComplex(3)*b*d + compute::toComplex(12)*a*e;
    const complex disD =
        compute::toComplex(64)*pow(a,3)*e - C4*C4*a*a*c*c + C4*C4*a*b*b*c -
        C4*C4*a*a*b*d - compute::toComplex(3)*pow(b,4);

    if (compute::round(disDelta0.real()) == 0 && compute::round(disD.real()) == 0) {

        // all four roots are equal (quadruple root)
        const complex root = -b / (C4*a);
        for (uint8_t i = 0; i < 4; ++i) {

            // even though all four roots are equal, newRoot must be
            // re-initialized each time because its numerical value changes
            const Roots newRoot = solution(root, inConditions);
            rootsOfCharEquation.push_back(newRoot);
        }
    }
    else // other types of roots
        canBeSolved = false;

    const QString equationType =
        QStringLiteral("quartic equation (discriminant ") +
        stateOfDiscriminant(discriminant.real()) + ')';

    return equationType;
}

const QString CharEq::rootsOfPolynomialEquation(const QVector<Token> & inTokens) {

    return ("polynomial equation of " + QString::number(inTokens.size()-1) + "th degree");
}
