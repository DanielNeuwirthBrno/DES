#include "matrix.h"
#include "ui/ui_matrix.h"


MatrixWindow::MatrixWindow(const QChar inDiffVariable,
                           const Matrix<double_t> & inMatrix, QWidget * parent):
    QDialog(parent), ui(new Ui_MatrixWindow) {

    ui->setupUi(this, inDiffVariable, inMatrix);

    connect(ui->pushButton_OK, SIGNAL(clicked()), this, SLOT(close()));
}

MatrixWindow::~MatrixWindow() {

    delete ui;
}
