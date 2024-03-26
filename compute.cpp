#include "compute.h"

#include <cmath>
#include <QFile>
#include <QFontMetrics>


complex compute::cbrt(const complex inNumber) {

    double_t a = inNumber.real();
    double_t b = inNumber.imag();

    double_t r = std::sqrt(a*a + b*b);
    double_t fi = std::atan(b/a);

    a = std::cbrt(r) * std::cos(fi/3);
    b = std::cbrt(r) * std::sin(fi/3);

    // quadrant adjustment
    if (inNumber.real() < 0)
        a = -a;

    complex result(a,b);

    return result;
}

complex compute::round(const complex inNumber, const uint8_t inDecPlaces) {

    const double_t realPart = compute::round(inNumber.real(),inDecPlaces);
    const double_t imagPart = compute::round(inNumber.imag(),inDecPlaces);

    return (complex(realPart,imagPart));
}

double_t compute::round(const double_t inNumber, const uint8_t inDecPlaces) {

    const long int modifier = std::pow(10,inDecPlaces);
    double_t roundedNumber = std::round(inNumber*modifier);
    roundedNumber /= modifier;
    return (roundedNumber);
}

const QString format::expressDerivativesAsSymbols(const uint8_t inOrder) {

    QString expression = "";
    switch (inOrder) {
        case 0: break;
        case 3: expression += '\'';
        case 2: expression += '\'';
        case 1: expression += '\''; break;
        default: expression = "<sup>(" + QString::number(inOrder) + ")</sup>";
    }
    return (expression);
}

const QString format::loadRes(const QString & inResourceName) {

    const QString fromResource = ":/html/" + inResourceName + ".html/";
    QFile fromRes(fromResource);
    fromRes.open(QIODevice::ReadOnly | QIODevice::Text);

    return (fromRes.readAll());
}

const QString format::numberForOutput(const double_t inNumber, const bool remove) {

    QString textFormat = QString::number(inNumber, 'f');
    while (textFormat.endsWith('0'))
        textFormat.chop(1);
    if (textFormat.endsWith('.'))
        textFormat.chop(1);
    if (remove && textFormat.endsWith('1'))
        textFormat.chop(1);
    return (textFormat);
}

void format::setElementWidth(QLineEdit * const inElement) {

    QFontMetrics metrics(inElement->font());
    inElement->setFixedWidth(metrics.width(inElement->text())+10);
    return;
}
