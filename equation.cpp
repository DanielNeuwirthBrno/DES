#include "compute.h"
#include "equation.h"
#include "general.h"

#include <algorithm>
#include <cmath>
#include <QRegularExpression>
#include <QStringList>


Equation::Equation
    (const bool inApplyConds, const QVector<InitConditions> inInitConditions,
     const QChar inDiffVar, const uint8_t inDiffOrder,
     const QString & inLeftSide, const QString & inRightSide):

    // construct tokens vector with specified number of elements
    // characteristic Equation has 'inDiffOrder' roots
    characteristicEquation(inDiffOrder), tokens(inDiffOrder+1),
    diffVariable(inDiffVar), order(inDiffOrder), type("0000"),
    equationLeftSide(inLeftSide), equationRightSide(inRightSide) {

    initialConditions.second = inInitConditions;
    initialConditions.first = inApplyConds;
}

bool Equation::isLinear() const {

    for (auto it: tokens)
        if (it.getPower() != 1)
            return (false);
    return (true);
}

bool Equation::isHomogeneous() const {

    if (equationRightSide == "0" || equationRightSide.isEmpty())
        return (true);
    return (false);
}

bool Equation::isOrdinary() const {

    QString diffVariables = "";
    for (auto it: tokens)
        if (it.isDerived())
            diffVariables += it.getVariable();
    return (determine::allCharsTheSame(diffVariables));
}

bool Equation::moveToRightSide() {

    int16_t checkPointFrom = 0, checkPointTo = -1;
    QString partToMove = "";

    // cycle through the equation as long as an element to be copied is
    // identified or the end of the equation is encountered
    while (checkPointTo < equationLeftSide.length()) {

        // determine where current element (substring) ends
        checkPointTo = equationLeftSide.indexOf(QRegularExpression("[+|-]"),
                                                checkPointFrom);
        if (checkPointTo == -1)
            checkPointTo = equationLeftSide.length();

        // copy identified element (substring) for further examination
        if (checkPointFrom != 0)
            --checkPointFrom;
        partToMove = equationLeftSide.mid(checkPointFrom,
                                          checkPointTo-checkPointFrom);

        // is it a double_t type?
        bool test = false;
        partToMove.toDouble(&test);

        // if either of these two conditions is satisfied ...
        if ((partToMove.indexOf(diffVariable) == -1 &&
             partToMove.indexOf(QRegularExpression("[a-z]")) != -1) || test ) {

            // ... move current element to the right side of the equation
            equationLeftSide.remove(checkPointFrom,
                                    checkPointTo-checkPointFrom);

            // clear contents of the right side if it's equal to "0"
            if (equationRightSide == "0")
                equationRightSide.clear();

            // the left side of the equation need not start with a '+' sign
            if (equationLeftSide.startsWith('+'))
                equationLeftSide.remove(0,1);

            // reverse sign according to mathematical rules
            if (partToMove.startsWith('-')) {

                if (equationRightSide.isEmpty())
                    // remove '+' sign at the beginning of the expression
                    partToMove.remove(0,1);
                else
                    // just reverse sign
                    partToMove.replace(0,1,'+');
            }
            else {

                if (partToMove.startsWith('+'))
                    // just reverse sign
                    partToMove.replace(0,1,'-');
                else
                    // prepend '-' sign at the beginning of the expression
                    partToMove.prepend('-');
            }

            // add the identified part of the equation to its right side
            equationRightSide.append(partToMove);

            // returns true if any part of the equation has been moved
            return (true);
        }
        else
            checkPointFrom = checkPointTo+1;
    }
    return (false);
}

QString Equation::tokenDelimiters(QString inEquation,
                                  const int16_t currentPos) const {

    // insert '+' sign before any '-' sign to simplify process of breaking
    // the equation down to individual elements during its parsing
    int16_t position = inEquation.lastIndexOf('-',currentPos);
    if (position == -1)
        return (inEquation);
    else
        return (tokenDelimiters(inEquation.insert(position, '+'), position));
}

bool Equation::parseEquation(const QString & inRegExp) {

    // split (left side of) the equation into individual elements
    QStringList tokensAsStrings =
        tokenDelimiters(equationLeftSide).split('+', QString::SkipEmptyParts);
    // pointer for manipulating data inside tokens vector
    Token * const pTokens = tokens.data();

    for (auto it = tokensAsStrings.begin();
              it != tokensAsStrings.end(); ++it) {

        // determine position of VWADWRT
        QChar variable;
        uint8_t elementOrder = 0;
        uint8_t numberOfCharsToRemove = 1;
        int16_t power = 1;

        // if VWADWRT in proper format is found
        int16_t positionOfVariable =
            it->indexOf(QRegularExpression(diffVariable+inRegExp));

        if (positionOfVariable != -1) {

            // determine order of differentation
            elementOrder = (it->at(positionOfVariable+2)).digitValue();
            numberOfCharsToRemove = 4;
        }
        else {
            // if other VWADWRT is found than indicated by diffSign variable
            positionOfVariable = it->indexOf(QRegularExpression("[a-z]"+inRegExp));

            if (positionOfVariable != -1) {
                return (false);
            }
            else {
                // if no derivative (just the VWADWRT itself) is found
                positionOfVariable = it->indexOf(QRegularExpression(diffVariable));
                if (positionOfVariable == -1)
                    return (false);
            }
        }
        variable = it->at(positionOfVariable);

        // determine position of power (of VWADWRT)
        int16_t positionOfPower = it->indexOf(QRegularExpression("\\^\\d+"));

        if (positionOfPower == positionOfVariable+numberOfCharsToRemove) {

            ++numberOfCharsToRemove; // remove ^ sign
            uint8_t numberOfDigits = 1;
            QString powerAsString;

            do {
                powerAsString += it->at(positionOfPower+numberOfDigits);
                ++numberOfCharsToRemove;
                ++numberOfDigits;

            } while (it->indexOf(QRegularExpression("\\^\\d{"+
                     QString::number(numberOfDigits)+"}")) == positionOfPower);
            power = powerAsString.toInt();
        }
        it->remove(positionOfVariable, numberOfCharsToRemove);

        // determine coefficient; basically it's (at least it should be) just
        // the remainder of the original input after performing above-mentioned
        // operations of extracting VWADWRT, its order and its power
        double_t coefficient = 1.0;
        bool hasConstCoeffs = true;

        if (!it->isEmpty()) {

            if (*it != "-") {

                // check for any remaining non-digit characters
                coefficient *= it->toDouble(&hasConstCoeffs);

                if (!hasConstCoeffs) {

                    // search for non-constant coefficients
                    // => equation could be parsed but is not solvable (by this software)
                    if (it->indexOf(QRegularExpression("[a-z]")) != 0 ||
                        it->length() != 1)
                        return (false); // parsing error
                    // else // proceed (hasConstCoeffs set to false)
                }
            }
            else
                coefficient = -1.0;
        }
        // set equation type - with const/non-const coefficients
        type.set(3, hasConstCoeffs);

        // check if any elements of the same order as is currently
        // being subject of inspection have been already processed;
        // if yes, add up their coefficients
        const uint8_t currentRow = tokens.size()-elementOrder-1;
        const double_t lastCoef = pTokens[currentRow].getCoefficient();

        // update token information inside vector
        pTokens[currentRow] = Token(coefficient+lastCoef, variable,
                                    elementOrder, power);
    }

    delete pTokens;
    return (true);
}

const QString Equation::convertToOrdinalNumber(const uint8_t inNumber) const {

    const QString ending = "stndrdth";
    const uint8_t from = (inNumber <= 4) ? (inNumber*2-2) : 6;

    return (QString::number(inNumber)+ending.mid(from,2));
}

const QString Equation::buildDerivative(
    const QString & inSolution, const double_t inCoeff, const uint8_t inConst,
    const uint8_t inPower, QString inSign) {

    // do not display expressions with null coefficient
    if (compute::round(inCoeff) == 0)
        return "";

    if (inCoeff < 0)
         inSign.resize(0);

    return (inSign + format::numberForOutput(inCoeff, true) +
          constantsInResult(inConst) + notZeroPowerElems(inPower) + inSolution);
}

double_t Equation::evalDerivativeRealPart(const double_t rootReal,
    const uint8_t derivativeOrder, const double_t inCoeff) const {

    const double_t evaluatedValue = (initialConditions.first) ? inCoeff *
    std::exp(rootReal*initialConditions.second.at(derivativeOrder).getXValue()) : 0;

    return (evaluatedValue);
}

double_t Equation::evalDerivativeImagPart(const double_t rootImag,
    const uint8_t derivativeOrder, double_t (*inFuncPointer)(const double_t)) const {

    const double_t evaluatedValue = (initialConditions.first) ? inFuncPointer(
        rootImag*initialConditions.second.at(derivativeOrder).getXValue()) : 0;

    return (evaluatedValue);
}

const QString Equation::determineType(bool & canBeSolved) {

    // 1] order
    QString description = convertToOrdinalNumber(order)+ " order ";

    // 2] is linear?
    type.set(0, isLinear());
    if (type.test(LINEAR))
        description += "linear ";
    else {
        description += "nonlinear ";
        canBeSolved = false;
    }

    // 3] is homogeneous?
    type.set(1, isHomogeneous());
    if (type.test(HOMOGENEOUS))
        description += "homogeneous ";
    else
        description += "nonhomogeneous ";

    // 4] is ordinary?
    type.set(2, isOrdinary());
    if (type.test(ORDINARY))
        description += "ordinary ";
    else
        canBeSolved = false;

    description += "differential equation ";

    // 5] determine coefficients
    // appropriate bit already set during parsing stage
    if (type.test(CONST_COEFF))
        description += "with constant coefficients";
    else {
        description += "with non-constant coefficients ";
        canBeSolved = false;
    }

    return (description);
}

bool Equation::generalSolution(std::pair<bool, Matrix<double_t>> & fundamentalMatrix) {

    // [1] calculate roots of characteristic equation

    bool canBeSolved = true;
    QString descriptionOfCharEq = "not determined";

    // call appropriate function according to the degree of the char.equation
    switch (order) {

        case 1: canBeSolved = false; break;
        case 2: // quadratic
          descriptionOfCharEq = characteristicEquation.rootsOfQuadraticEquation(
                                    tokens, initialConditions);
            break;
        case 3: // cubic
          descriptionOfCharEq = characteristicEquation.rootsOfCubicEquation(
                                    tokens, initialConditions);
            break;
        case 4: // quartic
          descriptionOfCharEq = characteristicEquation.rootsOfQuarticEquation(
                                    tokens, initialConditions, canBeSolved);
            break;
        default: // higher degrees
            descriptionOfCharEq =
                characteristicEquation.rootsOfPolynomialEquation(tokens);
            canBeSolved = false;
    }
    saveStep("characteristic equation: " + descriptionOfCharEq);
    if (!canBeSolved)
        return (false);

    // [2] sort roots
    // std::stable_sort (unlike std::sort and std::partial_sort) preserves
    // order of equal elements; this is essential in case of multiple roots

    // declare non_const iterators
    auto from = characteristicEquation.getRootsModifiable().begin();
    auto to = characteristicEquation.getRootsModifiable().end();
    // phase one: move real roots towards the beginning (top) of the vector,
    // and complex roots to the end (bottom) of the vector
    sortRealRootsOverComplex customSortPhaseOne;
    std::stable_sort(from, to, customSortPhaseOne);
    // phase two: sort real roots (ie. only part of the vector) [ASC]
    to = from + characteristicEquation.noOfRealRoots();
    sortRealRootsInAscOrder customSortPhaseTwo;
    std::stable_sort(from, to, customSortPhaseTwo);

    // save description of roots
    saveStep(characteristicEquation.rootsDescription());

    // save values of individual roots
    QString roots = QStringLiteral("roots of characteristic equation: ");
    QString solutions = QStringLiteral("solutions of differential equation: ");

    for (auto it: characteristicEquation.getRoots()) {

        roots += QString::number(it.getRealPartOfRoot(),'f');

        if (it.getImagPartOfRoot() != 0) {
            if (it.getImagPartOfRoot() > 0)
                roots += "+";
            roots += QString::number(it.getImagPartOfRoot(),'f') + "i";
        }
        roots += ", ";
        solutions += it.getSolution() + ", ";
    }
    roots.chop(2);
    saveStep(roots);
    solutions.chop(2);
    saveStep(solutions);

    // save (not yet approved) general solution
    const QString descriptionOfSolution =
        QStringLiteral("according to Principle of superposition: ") +
        principleOfSuperposition();
    saveStep(descriptionOfSolution);

    // if initial conditions are not set (or set but not applied)
    if (!applyInitialConditions())
        return (true);

    // if initial conditions are set and applied

    // [3] compute derivatives
    const uint8_t noOfDerivatives = order-1;
    computeDerivatives(noOfDerivatives);

    // save description of derivatives
    for (uint8_t row = 0; row < noOfDerivatives; ++row) {

        QString derivativeDescription =
            convertToOrdinalNumber(row+1) + QStringLiteral(" derivative: ") +
            diffVariable + format::expressDerivativesAsSymbols(row+1) +
            QStringLiteral(" = ");

        for (auto it: characteristicEquation.getRoots()) {

            if (!it.getDerivative(row).startsWith('-') &&
                !it.getDerivative(row).isEmpty() &&
                !derivativeDescription.endsWith("= "))
                derivativeDescription += '+';
            derivativeDescription += it.getDerivative(row);
        }
        saveStep(derivativeDescription);
    }

    // [4] compute wronskian
    double_t wronskian = 0;
    QString wronskianDescription = QStringLiteral("wronskian: W = ");
    if (!computeWronskian(wronskian, fundamentalMatrix.second)) {
        wronskianDescription += QStringLiteral("cannot be determined");
        saveStep(wronskianDescription);
        return (false);
    }
    fundamentalMatrix.first = true;
    wronskianDescription += format::numberForOutput(wronskian) +
                            QStringLiteral(" (see fundamental matrix)");
    saveStep(wronskianDescription);

    if (wronskian == 0) {

        saveStep(QStringLiteral("solutions don't form general solution"));
        return (true);
    }
    saveStep(QStringLiteral
        ("wronskian is non-zero => solutions form a general solution"));
    const QString generalSolutionDescription =
        QStringLiteral("general solution: ") + diffVariable +
        QStringLiteral(" = ") + finalResult;
    saveStep(generalSolutionDescription);

    // [5] plug in initial conditions and compute coefficients
    QVector<double_t> constants(order);
    determineValuesOfConstants(wronskian, constants);

    QString constantsDescription = QStringLiteral("constants: ");
    QString constantNumber, constantSymbol;
    uint8_t constantNo = 0;

    for (auto it : constants) {

        // display each individual constant
        constantsDescription += (constantSymbol =
            constantsInResult(++constantNo)) + QStringLiteral(" = ") +
            (constantNumber = format::numberForOutput(it));
        if (constantNo != constants.size())
            constantsDescription += ", ";

        // build general solution = inject numerical values of constants
        const QString sign =
            (!finalResult.startsWith(constantSymbol) && it > 0)
                ? "+" : ((it == -1) ? "-" : "");
        // don't show multiples of one
        if (std::abs(it) == 1)
            constantNumber.resize(0);
        // if constant is zero => remove whole element (until next + sign)
        if (it == 0) {
            const int16_t pos = finalResult.indexOf(constantSymbol);
            constantSymbol = finalResult.mid(pos,finalResult.indexOf('+',pos));
            constantNumber.resize(0);
        }
        if (!finalResult.startsWith(constantSymbol))
            constantSymbol.prepend('+');

        finalResult.replace(constantSymbol, sign+constantNumber);
    }
    saveStep(constantsDescription);

    // [6] save actual solution
    const QString generalSolution = QStringLiteral("actual solution: ") +
                                    diffVariable + QStringLiteral(" = ");
    if (finalResult.isEmpty())
        finalResult.append('0');
    saveStep(generalSolution+finalResult);

    return (true);
}

Equation::gonFunction Equation::changeGonFunction(QString & inResultPartB) {

    if (inResultPartB.contains("sin")) {
        inResultPartB.replace("sin","cos");
        return (SIN);
    }
    if (inResultPartB.contains("cos")) {
        inResultPartB.replace("cos","sin");
        return (COS);
    }
    return (UNKNOWN);
}

bool Equation::computeDerivatives(const uint8_t noOfDerivatives) {

    uint8_t constantNo = 1;
    QString lastRoot = "";
    uint8_t lastRootRepeats = 0;

    // cycle through all roots
    for (auto & it: characteristicEquation.getRootsModifiable()) {

            // coefficients
            double_t coeff1 = 1, coeff2 = -1;
            QVector<double_t> coeffsRepeatedRoots; // for multiple roots
            uint8_t noOfComplexRootsWithZeroRealPart = 0;
            const double_t real = it.getRealPartOfRoot();
            const double_t imag = it.getImagPartOfRoot();

            // function pointers to goniometric functions
            double_t (*gonFunctionFirstPart)(double_t) = &(std::sin);
            double_t (*gonFunctionSecondPart)(double_t) = &(std::cos);

            // compute all derivatives up to order-1
            for (uint8_t row = 0; row < noOfDerivatives; ++row) {

                // description of derivatives - its string representation and
                // numerical value (in case initial conditions have been applied)
                QString derivative = "";
                double_t derivativeAtInit = 0;

                // root is NULL
                if (compute::round(real) == 0 && compute::round(imag) == 0) {

                    // root is NULL and MULTIPLE
                    if (lastRoot == it.getSolution(Roots::ONLY_PART_A)) {

                        // first derivative - initial coefficients
                        if (row == 0) {
                            coeff1 = ++lastRootRepeats;
                            coeff2 = coeff1-1; // exponent
                        }

                        // build derivative
                        derivative = buildDerivative(
                            it.getSolution(Roots::TEST_REAL_ROOT_FIRST),
                            coeff1, constantNo, coeff2);

                        // modify coefficients
                        coeff1 *= coeff2--;
                    }
                    // root is NULL and DISTINCT
                    else {
                        // derivative is null - no need to assign any value

                        if (row+1 == noOfDerivatives)
                            lastRoot = it.getSolution(Roots::ONLY_PART_A);
                        else {
                            lastRoot = "";
                            lastRootRepeats = 0;
                        }
                    }
                }

                // root is REAL
                if (compute::round(real) != 0 && compute::round(imag) == 0) {

                    coeff1 *= real;
                    double_t currentCoeff = 0;

                    // root is REAL and MULTIPLE
                    if (lastRoot == it.getSolution(Roots::ONLY_PART_A)) {

                        // compute first derivative
                        if (row == 0) {

                            // count for how many times is this root repeated
                            currentCoeff = ++lastRootRepeats;
                            // ... and resize storage for coeffs accordingly
                            coeffsRepeatedRoots.resize(lastRootRepeats+1);
                            QString sign = "";

                            // build derivative (which comprises of two parts)
                            for (uint8_t i=0; i<2; ++i) {

                                derivative += buildDerivative(
                                    it.getSolution(Roots::ONLY_PART_A),
                                    currentCoeff, constantNo,
                                    lastRootRepeats+i-1, sign);

                                if (!derivative.isEmpty())
                                    sign = "+";
                                currentCoeff = coeff1;
                            }

                            // evaluate derivative at initial point
                            if (initialConditions.first)
                                derivativeAtInit =
                                evalDerivativeRealPart(real, row, lastRootRepeats)*
                                    std::pow(initialConditions.second.at(row).
                                             getXValue(), lastRootRepeats-1)+
                                evalDerivativeRealPart(real, row, currentCoeff)*
                                    std::pow(initialConditions.second.at(row).
                                             getXValue(), lastRootRepeats);

                            // determine coefficients essential
                            // for computing higher derivatives
                            coeffsRepeatedRoots[row] =
                                (real+coeff1)*lastRootRepeats;
                            coeffsRepeatedRoots[row+1] =
                                lastRootRepeats*(lastRootRepeats-1);
                        }

                        // compute higher derivatives
                        else {

                            // although one-character long, string variable
                            // must be used to avoid adding a null character
                            QString sign = "";

                            // compound derivative of the first through
                            // to the one-before-the-last part of the expression
                            for (uint8_t i = lastRootRepeats; i>0; --i) {

                                // call buildDerivative function with Power
                                // variable explicitly set to lastRootRepeats-i
                                derivative += buildDerivative(
                                    it.getSolution(Roots::ONLY_PART_A),
                                    coeffsRepeatedRoots.at(i-1), constantNo,
                                    lastRootRepeats-i, sign);

                                // evaluate derivative at initial point (first part)
                                if (initialConditions.first)
                                    derivativeAtInit +=
                                        evalDerivativeRealPart(real, row,
                                        coeffsRepeatedRoots.at(i-1))*
                                        std::pow(initialConditions.second.
                                        at(row).getXValue(), lastRootRepeats-i);

                                if (!derivative.isEmpty())
                                    sign = "+";

                                // determine new coefficients
                                // for computing higher derivatives
                                currentCoeff =
                                    (i > 1) ? coeffsRepeatedRoots[i-2] : coeff1;
                                coeffsRepeatedRoots[i-1] =
                                    real*coeffsRepeatedRoots.at(i-1) +
                                    (lastRootRepeats-i+1)*currentCoeff;
                            }

                            // derivative of the last part of the expression
                            derivative += buildDerivative(
                                it.getSolution(Roots::ONLY_PART_A), coeff1,
                                constantNo, lastRootRepeats, sign);

                            // evaluate derivative at initial point (add second part)
                            if (initialConditions.first)
                                derivativeAtInit +=
                                    evalDerivativeRealPart(real, row, coeff1)*
                                    std::pow(initialConditions.second.at(row).
                                             getXValue(), lastRootRepeats);
                        }
                    }
                    // root is REAL and DISTINCT
                    else {

                        // build derivative
                        derivative =
                            buildDerivative(it.getSolution(Roots::ONLY_PART_A),
                                            coeff1, constantNo);

                        // evaluate derivative at initial point
                        derivativeAtInit =
                            evalDerivativeRealPart(real, row, coeff1);

                        if (row+1 == noOfDerivatives)
                            lastRoot = it.getSolution(Roots::ONLY_PART_A);
                        else {
                            lastRoot = "";
                            lastRootRepeats = 0;
                        }
                    }
                }

                // root is COMPLEX (its real part is non-zero)
                if (compute::round(real) != 0 && compute::round(imag) != 0) {

                    // first derivative
                    if (row == 0)
                        { coeff1 *= real; coeff2 *= imag; }
                    else {
                        // higher derivatives
                        double_t lastCoeff1 = coeff1;

                        // coeff1[n-degree] = coeff2[(n-1)-degree] *
                        // imag part of root + coeff1[(n-1)-degree] * real part
                        coeff1 = coeff2 * imag + coeff1 * real;
                        // coeff2[n-degree] = coeff2[(n-1)-degree] *
                        // real part of root - coeff1[(n-1)-degree] * imag part
                        coeff2 = coeff2 * real - lastCoeff1 * imag;
                    }

                    // interchange sin/cos functions in second part of derivation
                    QString diffedResultPartB =  it.getSolution();
                    switch (changeGonFunction(diffedResultPartB)) {
                        case SIN: break;
                        case COS: std::swap(gonFunctionFirstPart,
                                            gonFunctionSecondPart);
                                  break;
                        case UNKNOWN:
                        default: return (false);
                    }

                    // build derivative
                    derivative =
                        buildDerivative(it.getSolution(), coeff1, constantNo);
                    const QString sign = (!derivative.isEmpty()) ? "+" : "";
                    derivative += buildDerivative(diffedResultPartB, coeff2,
                                                  constantNo, 0, sign);

                    // evaluate derivative at initial point
                    derivativeAtInit =
                        evalDerivativeRealPart(real, row, coeff1) *
                        evalDerivativeImagPart(std::abs(imag), row,
                                               gonFunctionFirstPart) +
                        evalDerivativeRealPart(real, row, coeff2) *
                        evalDerivativeImagPart(std::abs(imag), row,
                                               gonFunctionSecondPart);
                }

                // root is COMPLEX (its real part is zero)
                if (compute::round(real) == 0 && compute::round(imag) != 0) {

                    double_t coeffWithSign = coeff1 *= std::abs(imag);
                    QString diffedResult = it.getSolution(Roots::ONLY_PART_B);

                    // interchange functions only in case of odd derivatives
                    if (row % 2 == 0)
                        switch (changeGonFunction(diffedResult)) {
                            case COS: break;
                            case SIN: std::swap(gonFunctionFirstPart,
                                                gonFunctionSecondPart);
                                      break;
                            case UNKNOWN:
                            default: return (false);
                        }

                    // determine sign to be put in front of sin/cos functions
                    // -> sin is negative every 2nd and 3rd derivative
                    // -> cos is negative every 1st and 2nd derivative
                    if (noOfComplexRootsWithZeroRealPart % 4 == 1 ||
                       (diffedResult.contains("sin") &&
                        noOfComplexRootsWithZeroRealPart == 0) ||
                        (diffedResult.contains("cos") &&
                        noOfComplexRootsWithZeroRealPart == 2))
                        coeffWithSign = -coeffWithSign;

                    // build derivative
                    derivative = buildDerivative(diffedResult, coeffWithSign,
                                                 constantNo);
                    ++noOfComplexRootsWithZeroRealPart;

                    // evaluate derivative at initial point
                    derivativeAtInit =
                        evalDerivativeRealPart(real, row, coeff1)*
                        evalDerivativeImagPart(std::abs(imag), row,
                                               gonFunctionFirstPart);
                }

                // add to Roots (to vector of derivatives)
                Derivatives newDerivative(derivative, derivativeAtInit);
                it.addDerivative(newDerivative);
                derivative.resize(0); // erase description
            }
            ++constantNo;
    }
    return (true);
}

bool Equation::computeWronskian(double_t & wronskian, Matrix<double_t> & fMatrix) {

    // fMatrix has been empty so far
    fMatrix.setSize(order);

    // in case of a square matrix of 4th order is for computing its determinant
    // used a built-in class QMatrix4x4; unfortunately this can be initialized
    // only by plugging in an ordinary array (pointer); this array is filled in
    // in a column-major order (according to existing code it's easier to fill
    // it up this way); be it as it may it can be filled up in whatever order
    // you choose (row-major or column-major) because determinant of any square
    // matrix is equal to determinant of its transposed "cousin"
    float values[16]; float * p = values;

    for (uint8_t column = 0; column < fMatrix.size(); ++column) {

        // first row of fundamental matrix = individual roots
        *p++ = fMatrix.at(0,column) =
            characteristicEquation.getRoots().at(column).getNumValue();
        // other rows of individual matrix = derivatives of roots
        for (uint8_t row = 1; row < fMatrix.size(); ++row)
            *p++ = fMatrix.at(row,column) =
                characteristicEquation.getRoots().at(column).
                                       evalDerivativeAtInit(row-1);
    }

    switch (fMatrix.size()) {
        case 2:
        case 3: wronskian = fMatrix.determinant3(); break;
        case 4: wronskian = fMatrix.determinant4(values); break;
        default: return (false);
    }

    return (true);
}

bool Equation::determineValuesOfConstants(const double_t wronskian,
                                          QVector<double_t> & inConstants) {

    Matrix<double_t> matrix(inConstants.size());

    // square matrix of 4th order (see computeWronskian for explanation)
    float values[16];

    // repeat "number-of-constants" times
    for (uint8_t constantNo = 0; constantNo < inConstants.size(); ++constantNo) {

        float * p = values + constantNo * inConstants.size();
        // constantNo-th column contains function values at initial points
        for (uint8_t i = 0; i < matrix.size(); ++i)
            *p++ = matrix.at(i,constantNo) =
                initialConditions.second.at(i).getFunctionValue();

        p = values;
        for (uint8_t column = 0; column < matrix.size(); ++column) {

            // constantNo-th column is already filled in (see above)
            if (column == constantNo)
                { p+=4; continue; }
            // first row (from second column onwards) contains
            // function values of individual roots at initial points
            *p++ = matrix.at(0,column) =
                characteristicEquation.getRoots().at(column).getNumValue();

            // other rows (from second column onwards) contain
            // values of individual root's derivatives at initial points
            for (uint8_t row = 1; row < matrix.size(); ++row)
                *p++ = matrix.at(row,column) =
                    characteristicEquation.getRoots().at(column).
                    evalDerivativeAtInit(row-1);
        }
        double_t determinant = 0;

        switch (matrix.size()) {
            case 2:
            case 3: determinant = matrix.determinant3(); break;
            case 4: determinant = matrix.determinant4(values); break;
            default: return (false);
        }

        // Cramer's rule
        inConstants[constantNo] = determinant / wronskian;
    }
    return (true);
}

const QString Equation::principleOfSuperposition() {

    QString solution = diffVariable + QStringLiteral(" = ");
    QString repeatedRoot;
    uint8_t constantNo = 1, exponent = 0;
    double_t lastRealRoot = 0;

    for (auto it : characteristicEquation.getRoots()) {

        // check for repeated real roots
        if (lastRealRoot == it.getRealPartOfRoot() &&
            it.getImagPartOfRoot() == 0 && !repeatedRoot.isNull()) {

            repeatedRoot = "x";
            if (++exponent > 1)
                repeatedRoot += "<sup>"+QString::number(exponent)+"</sup>";
        }
        else {
            exponent = 0;
            repeatedRoot = "";
        }

        finalResult += constantsInResult(constantNo++) + repeatedRoot +
                     it.getSolution(Roots::TEST_REAL_ROOT_FIRST) + "+";
        lastRealRoot = it.getRealPartOfRoot();
    }

    finalResult.chop(1);
    solution += finalResult;

    return solution;
}
