#include "matrix.h"
#include "steps.h"
#include "ui/ui_steps.h"


StepsWindow::StepsWindow(const QChar diffVariable, const QString & steps,
    const std::pair<bool, Matrix<double_t>> & matrix, QWidget * parent):
    QDialog(parent, Qt::WindowMaximizeButtonHint), diffVariable(diffVariable),
    fundamentalMatrix(matrix.second), ui(new Ui_StepsWindow) {

    ui->setupUi(this, steps, matrix.first);

    connect(ui->pushButton_OK, SIGNAL(clicked()), this, SLOT(close()));
    connect(ui->pushButton_CopyToClipboard, SIGNAL(clicked()),
            ui->textEdit_ComputationSteps, SLOT(selectAll()));
    connect(ui->pushButton_CopyToClipboard, SIGNAL(clicked()),
            ui->textEdit_ComputationSteps, SLOT(copy()));
    connect(ui->pushButton_FundamentalMatrix, SIGNAL(clicked()),
            this, SLOT(fundamentalMatrix_clicked()));
}

StepsWindow::~StepsWindow() {

    delete ui;
}

int StepsWindow::fundamentalMatrix_clicked() {

    MatrixWindow fundamentalMatrixWindow(diffVariable, fundamentalMatrix, this);
    return fundamentalMatrixWindow.exec();
}
