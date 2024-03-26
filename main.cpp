/* Application name:    Differential Equations Solver (des.exe)
 * Current version:     0.55 - work in progress
 *
 * Author:              Daniel Neuwirth (d.neuwirth@tiscali.cz)
 * Created:             2016-10-14
 * Lastly modified:     2016-04-20
 * Time spent:          204 hrs (of development)
 *
 * IDE/framework:       Qt 5.7.0
 * Compiler:            MinGW 5.3.0 32-bit
 * Language standard:   C++11
 *
 * Objective:           a simple application for solving 2nd order, 3rd order
 *                      and some types of 4th order linear homogeneous ordinary
 *                      differential equations with constant coefficients;
 *                      with the option to apply initial conditions
 *                      /project for my bachelor thesis/
 *
 * Description:         N/A
 */

#include "mainwindow.h"

#include <QApplication>


int main(int argc, char * argv[])
{
    Q_INIT_RESOURCE(resources);

    QApplication application(argc, argv);
    MainWindow mainApplicationWindow;

    mainApplicationWindow.show();

    return application.exec();
}
