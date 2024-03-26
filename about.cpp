#include "about.h"
#include "ui/ui_about.h"


AboutWindow::AboutWindow(QWidget * parent):
    QDialog(parent), ui(new Ui_AboutWindow) {

    ui->setupUi(this);

    connect(ui->pushButton_OK, SIGNAL(clicked()), this, SLOT(close()));
}

AboutWindow::~AboutWindow() {

    delete ui;
}
