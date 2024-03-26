#include "conditions.h"
#include "ui/ui_conditions.h"


ConditionsWindow::ConditionsWindow(
    const QChar inDiffVariable, const uint8_t inNumberOfConditions,
    QVector<InitConditions> & inConditions, QWidget * parent):
    QDialog(parent), conditions(inConditions), ui(new Ui_ConditionsWindow) {

    ui->setupUi(this, conditions, inDiffVariable, inNumberOfConditions);

    connect(ui->pushButton_OK, SIGNAL(clicked()), this, SLOT(saveConditions()));
    connect(ui->pushButton_OK, SIGNAL(clicked()), this, SLOT(close()));
    connect(ui->pushButton_Cancel, SIGNAL(clicked()), this, SLOT(close()));
}

ConditionsWindow::~ConditionsWindow() {

    delete ui;
}

void ConditionsWindow::saveConditions() {

    uint8_t rowNumber = 0;
    // if initial conditions are being assigned for the first time
    if (conditions.isEmpty())
        conditions.resize(ui->individualRowsUI.size());

    rowNumber = 0;
    // save currently assigned values of initial conditions
    for (auto it: ui->individualRowsUI) {

        conditions[rowNumber++].setValues(
            it.lineEdit_ValueAtPoint->text().toDouble(),
            it.lineEdit_ValueOfDiffVar->text().toDouble());

    }
    return;
}
