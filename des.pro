#-----------------------------------------#
#                                         #
# Project created by QtCreator 2016-10-08 #
#                                         #
#-----------------------------------------#

QT        += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = des
TEMPLATE = app

SOURCES   += main.cpp \
             mainwindow.cpp \
             equation.cpp \
             token.cpp \
             general.cpp \
             compute.cpp \
             steps.cpp \
             about.cpp \
             conditions.cpp \
             matrix.cpp

HEADERS   += about.h \
             compute.h \
             conditions.h \
             equation.h \
             general.h \
             mainwindow.h \
             matrix.h \
             steps.h \
             token.h \
             # GUI
             ui/ui_about.h \
             ui/ui_conditions.h \
             ui/ui_mainwindow.h \
             ui/ui_matrix.h \
             ui/ui_steps.h

FORMS     +=

RESOURCES += resources.qrc

DISTFILES += notes.txt \
             solution.txt \
             todolist.txt
