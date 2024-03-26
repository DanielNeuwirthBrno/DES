#include "about.h"
#include "compute.h"
#include "conditions.h"
#include "equation.h"
#include "mainwindow.h"
#include "steps.h"
#include "ui/ui_mainwindow.h"

#include <memory>
#include <QFontMetrics>
#include <QIcon>
#include <QMessageBox>
#include <QRegularExpression>


const QString MainWindow::findDiffVarRegex = "\\{[1-9]\\}";
const QString MainWindow::allowedSpecialCharacters = "+-^. ";
const QString MainWindow::containsErrorBackground = "background-color:#ffb0ff;";
const QString MainWindow::noSolution = QStringLiteral
    ("This type of equation cannot be solved using this software.");

MainWindow::MainWindow(QWidget * parent): QMainWindow(parent),
    cursorPosition(0), ui(new Ui_MainWindow) {

    ui->setupUi(this);

    connect(ui->lineEdit_EquationLeftSide, SIGNAL(textChanged(QString)),
            this, SLOT(extractDiffVariable(QString)));
    connect(ui->lineEdit_EquationLeftSide, SIGNAL(textChanged(QString)),
            this, SLOT(preprocessing_l(QString)));
    connect(ui->lineEdit_EquationLeftSide, SIGNAL(textChanged(QString)),
            this, SLOT(enableSubmitButton()));
    connect(ui->lineEdit_EquationLeftSide, &QLineEdit::editingFinished,
            this, [this]() -> void { cursorPosition =
                  ui->lineEdit_EquationLeftSide->cursorPosition(); });
    // textEdited() SIGNAL used instead of textChanged() to prevent undesirable
    // check of the right side of the equation during parsing stage
    // (which is performed automatically in case the right side is left blank)
    connect(ui->lineEdit_EquationRightSide, SIGNAL(textEdited(QString)),
            this, SLOT(preprocessing_r(QString)));
    connect(ui->lineEdit_EquationRightSide, SIGNAL(textChanged(QString)),
            this, SLOT(enableSubmitButton()));
    connect(ui->lineEdit_DiffSign, SIGNAL(textChanged(QString)),
            this, SLOT(enableSubmitButton()));
    connect(ui->lineEdit_DiffSign, SIGNAL(textChanged(QString)),
            this, SLOT(enableAddConditionsButton()));
    connect(ui->lineEdit_DiffSign, SIGNAL(textChanged(QString)),
            this, SLOT(checkIfOrderEqualsNoOfConditions()));
    connect(ui->pushButton_AddConditions, SIGNAL(clicked()),
            this, SLOT(addInitialConditions()));
    connect(ui->pushButton_AddConditions, SIGNAL(clicked()),
            this, SLOT(setInitialConditions()));
    connect(ui->pushButton_RemoveConditions, &QPushButton::clicked,
            this, [this]() -> void { setOfInitConditions.clear(); });
    connect(ui->pushButton_RemoveConditions, SIGNAL(clicked()),
            this, SLOT(setInitialConditions()));
    connect(ui->pushButton_ApplyConditions, SIGNAL(toggled(bool)),
            this, SLOT(applyConditionsIcon(bool)));
    connect(ui->textEdit_Result, SIGNAL(textChanged()),
            this, SLOT(enableShowStepsButton()));
    connect(ui->pushButton_InsertDiff, SIGNAL(clicked()),
            this, SLOT(insertDiffAtCaret()));
    connect(ui->pushButton_SubmitEquation, SIGNAL(clicked()),
            this, SLOT(submitEquation_clicked()));
    connect(ui->pushButton_ClearEquation, SIGNAL(clicked()),
            this, SLOT(clearEquation_clicked()));
    connect(ui->pushButton_ShowSteps, SIGNAL(clicked()),
            this, SLOT(showSteps_clicked()));

    connect(ui->action_Exit, SIGNAL(triggered()), this, SLOT(action_Exit()));
    connect(ui->action_About, SIGNAL(triggered()), this, SLOT(action_About()));
}

MainWindow::~MainWindow() {

    delete ui;
}

const QString MainWindow::errorDescription(const syntaxError inError,
                                           const QChar inDiffVar) const {

    QString explainError = QStringLiteral("parsing error: ");
    switch (inError) {
        case UNRECOGNIZED_CHAR:
            explainError += QStringLiteral("character not recognized"); break;
        case NOT_MATCHING_BRACES:
            explainError +=
                QStringLiteral("empty, missing or mismatched brace(s)"); break;
        case DIFF_SIGN_MISPLACED:
            explainError +=
                QStringLiteral("element containing variable (") + inDiffVar +
                QStringLiteral(") to be differentiated with respect to misplaced");
            break;
        case NO_VARIABLE_TO_DIFF:
            explainError += QStringLiteral("missing variable"); break;
        default:
            explainError += "unknown error";
    }
    return (explainError);
}

MainWindow::syntaxError MainWindow::checkForSyntaxErrors
    (const QString & inEquation, const QString & inAllowedChars,
     const QChar inDiffVar) const {

    int8_t countBraces = 0;
    uint16_t countFrom = 0;

    for (auto it = inEquation.constBegin();
              it != inEquation.constEnd(); ++it, ++countFrom) {

        // check for unallowed characters
        if (!it->isDigit() && !it->isLower() &&
            inAllowedChars.indexOf(*it) == -1)
            return (UNRECOGNIZED_CHAR);

        // check for VWADWRT on the right side of the equation
        if (*it == inDiffVar)
            return (DIFF_SIGN_MISPLACED);

        // check for matching curly braces
        if (*it == '{') {
            // opening brace must be preceded by a letter
            if (it == inEquation.constBegin() || !((it-1)->isLower()))
                return (NO_VARIABLE_TO_DIFF);
            ++countBraces;
            countFrom = 0;
        }
        if (*it == '}') {
            // just one character [in range 1-9] is allowed between braces
            if (countFrom != 2 || !isCharAllowed(*(it-1)))
                countBraces = 0;
            --countBraces;
        }
        if (countBraces < 0)
            return (NOT_MATCHING_BRACES);
    }
    if (countBraces != 0)
        return (NOT_MATCHING_BRACES);

    return (OK);
}

// slots

bool MainWindow::extractDiffVariable(QString inEquationLeftSide) const {

    // if no VWADWRT found => clear diffSign and return false
    int16_t position = inEquationLeftSide.indexOf(
                       QRegularExpression(findDiffVarRegex));

    if (position < 1) { // position == 0 would result in Out of Bounds access
        ui->lineEdit_DiffSign->clear();
        return (false);
    }

    // search for all occurrences of VWADWRT
    QChar order = 0;
    QString diffVariables = "";
    do {
        diffVariables += inEquationLeftSide[position-1];
        // detect highest order
        if (inEquationLeftSide[position+1] > order)
            order = inEquationLeftSide[position+1];

        inEquationLeftSide.remove(position-1,4);
    } while ((position = inEquationLeftSide.indexOf(
                         QRegularExpression(findDiffVarRegex))) > 0);

    QChar diffVar = multipleDiffVars;
    bool exit = true;

    if ((exit = determine::allCharsTheSame(diffVariables)))

        // if exactly one distinct VWADWRT found
        // => display diffSign and return true
        diffVar = diffVariables.at(0);\
    else
        // if more than one distinct VWADWRT found
        // => display warning and return false
        ui->lineEdit_EquationType->setText(QStringLiteral(
            "equation contains more than one variable "\
            "you are differentiating with respect to "));

    diffVariables = 'd'+diffVar+"[]";
    diffVariables.insert(diffVariables.length()-1,order);
    ui->lineEdit_DiffSign->setText(diffVariables);
    return (exit);
}

void MainWindow::preprocessing_l(const QString & inEquation) const {

    // check for syntax errors
    syntaxError typeOfError = UNSPECIFIED;
    if ((typeOfError = checkForSyntaxErrors(inEquation,
                       allowedSpecialCharacters+"{}")) != OK) {

        // changed background colour indicates check failure
        ui->lineEdit_EquationLeftSide->setStyleSheet(containsErrorBackground);
        // show explanation in equationType field
        ui->lineEdit_EquationType->setText(errorDescription(typeOfError));
    }
    else {
        ui->lineEdit_EquationLeftSide->
            setStyleSheet(QStringLiteral("background-color: white;"));
        ui->lineEdit_EquationType->clear();
        // call right-side-check in case it has already been containing
        // current input recognized as VWADWRT
        preprocessing_r(ui->lineEdit_EquationRightSide->text());
    }
    return;
}

void MainWindow::preprocessing_r(const QString & inEquation) const {

    // check for syntax errors
    syntaxError typeOfError = UNSPECIFIED;
    QChar diffVariable = '\0';
    if (!ui->lineEdit_DiffSign->text().isEmpty())
        diffVariable = ui->lineEdit_DiffSign->text().at(1);
    if ((typeOfError = checkForSyntaxErrors(inEquation,
                       allowedSpecialCharacters, diffVariable)) != OK) {

        // changed background colour indicates check failure
        ui->lineEdit_EquationRightSide->setStyleSheet(containsErrorBackground);
        // show explanation in equationType field
        ui->lineEdit_EquationType->setText(errorDescription(typeOfError, diffVariable));
    }
    else {
        ui->lineEdit_EquationRightSide->
            setStyleSheet(QStringLiteral("background-color: white;"));
        ui->lineEdit_EquationType->clear();
    }
    return;
}

void MainWindow::insertDiffAtCaret() const {

    const QString newText =
        ui->lineEdit_EquationLeftSide->text().left(cursorPosition) + "{1}" +
        ui->lineEdit_EquationLeftSide->text().mid(cursorPosition);
    ui->lineEdit_EquationLeftSide->setText(newText);
    ui->lineEdit_EquationLeftSide->setCursorPosition(cursorPosition+1);
    ui->lineEdit_EquationLeftSide->setSelection(cursorPosition+1,1);
    ui->lineEdit_EquationLeftSide->setFocus();
    return;
}

void MainWindow::enableSubmitButton() const {

    const bool disable = (ui->lineEdit_DiffSign->text().isEmpty()) ||
                         (ui->lineEdit_EquationLeftSide->styleSheet().
                          contains(containsErrorBackground)) ||
                         (ui->lineEdit_EquationRightSide->styleSheet().
                          contains(containsErrorBackground));

    ui->pushButton_SubmitEquation->setEnabled(!disable);
    return;
}

void MainWindow::enableAddConditionsButton() const {

    const bool enabled = !ui->lineEdit_DiffSign->text().isEmpty() &&
                       ui->lineEdit_DiffSign->text().at(1) != multipleDiffVars;
        ui->pushButton_AddConditions->setEnabled(enabled);

    return;
}

void MainWindow::checkIfOrderEqualsNoOfConditions() {

    // if number of initial conditions (rows in QVector)
    // does not equal the order of the equation ...
    if (!setOfInitConditions.isEmpty() && (ui->lineEdit_DiffSign->text().isEmpty()
        || setOfInitConditions.size() != ui->lineEdit_DiffSign->text().at(3))) {

        // ... erase initial conditions ...
        ui->pushButton_RemoveConditions->click();

        // ... and display a warning
        QMessageBox messageBox(QMessageBox::Warning, QStringLiteral("Warning"),
                 QStringLiteral("<b>Initial conditions have been erased.</b>"),
                               QMessageBox::Ok, this);
        messageBox.setInformativeText(QStringLiteral(
            "According to a recent change\n"\
            "in the order of your equation,\n"\
            "the original set of initial conditions\n"\
            "was not longer valid (or applicable)\n"\
            "and therefore has been erased."));
        messageBox.exec();
    }
    return;
}

int MainWindow::addInitialConditions() {

    // extract VWADWRT
    const QChar diffVariable = ui->lineEdit_DiffSign->text().at(1);
    // highest order = number of required conditions
    const uint8_t numberOfConditions =
        ui->lineEdit_DiffSign->text().at(3).digitValue();
    // if no conditions have been set yet, reserve storage space
    if (setOfInitConditions.isEmpty())
        setOfInitConditions.reserve(numberOfConditions);

    ConditionsWindow initialConditionsWindow(diffVariable, numberOfConditions,
                                             setOfInitConditions, this);
    return initialConditionsWindow.exec();
}

void MainWindow::setInitialConditions() {

    bool set = !setOfInitConditions.isEmpty();
    QString label = set ? QStringLiteral("initial conditions have been set") :
                          QStringLiteral("initial conditions not set");

    ui->pushButton_RemoveConditions->setEnabled(set);
    ui->lineEdit_InitialConditionsState->setText(label);
    format::setElementWidth(ui->lineEdit_InitialConditionsState);
    ui->pushButton_ApplyConditions->setEnabled(set);
    ui->pushButton_ApplyConditions->setChecked(set);

    return;
}

void MainWindow::applyConditionsIcon(const bool isChecked) const {

    QIcon * iconOK = new QIcon();
    if (isChecked)
        iconOK->addFile(QStringLiteral(":/icons/icons/ok.ico"),
                        QSize(), QIcon::Normal, QIcon::Off);
    ui->pushButton_ApplyConditions->setIcon(*iconOK);
    return;
}

void MainWindow::enableShowStepsButton() const {

    ui->pushButton_ShowSteps->
        setEnabled(!ui->textEdit_Result->toPlainText().isEmpty());
    return;
}

// process submitted equation
void MainWindow::submitEquation_clicked() {

    // if empty, assume that the right side of the equation is equal to zero
    if (ui->lineEdit_EquationRightSide->text().trimmed().isEmpty())
       ui->lineEdit_EquationRightSide->setText("0");

    // create Equation object with supplied parameters
    diffEquationToSolve = std::unique_ptr<Equation>(new Equation(
        ui->pushButton_ApplyConditions->isChecked(),
        setOfInitConditions, ui->lineEdit_DiffSign->text().at(1),
        ui->lineEdit_DiffSign->text().at(3).digitValue(),
        // remove whitespace from both sides of the equation
        ui->lineEdit_EquationLeftSide->text().remove(' '),
        ui->lineEdit_EquationRightSide->text().remove(' ')));
    diffEquationToSolve->saveStep(QStringLiteral("original input: ")+
                                 diffEquationToSolve->getEquationComplete());

    // VWADWRT could not be unambiguously determined
    if (diffEquationToSolve->getDiffVariable() == multipleDiffVars) {

        diffEquationToSolve->saveStep(QStringLiteral(
            "differentiate according to (more variables)"));
        diffEquationToSolve->saveStep(noSolution);
        ui->textEdit_Result->setPlainText(noSolution);
        return;
    }

    diffEquationToSolve->saveStep(QStringLiteral(
        "differentiate according to ")+diffEquationToSolve->getDiffVariable());

    // any parts of the equation containing other variables than VWADWRT
    // move to the right side of the equation
    while (diffEquationToSolve->moveToRightSide()) {

        ui->lineEdit_EquationLeftSide->
            setText(diffEquationToSolve->getEquationLeftSide());
        ui->lineEdit_EquationRightSide->
            setText(diffEquationToSolve->getEquationRightSide());

        // save (intermediate) state of regrouping stage
        diffEquationToSolve->saveStep(QStringLiteral(
            "regrouping: ")+diffEquationToSolve->getEquationComplete());

        // moving an element to the right side of the equation could have
        // introduced a new detectable error on the right side
        // note1: rules for both sides, e.g. allowed chars, are different
        // note2: function preprocessing_r is called automatically
        if (!ui->lineEdit_EquationType->text().isEmpty() ||
            ui->lineEdit_EquationLeftSide->text().isEmpty())
            return;
    }

    // reformat equation (only for output)
    int16_t position = -1;
    QString equationReformatted = diffEquationToSolve->getEquationLeftSide();
    // express derivatives as symbols [' signs or (n) superscripts]
    while ((position = equationReformatted.
           indexOf(QRegularExpression(findDiffVarRegex))) != -1) {

        equationReformatted.replace(position,3,
            format::expressDerivativesAsSymbols(
                equationReformatted.at(position+1).digitValue()));
    }
    diffEquationToSolve->saveStep(
        QStringLiteral("symbolically expressed as: ")+equationReformatted
                       +" = "+diffEquationToSolve->getEquationRightSide());

    // display initial conditions (if any)
    QString conditions = "initial conditions: ";
    if (!diffEquationToSolve->getInitialConditions().isEmpty()) {

        uint8_t rowNumber = 0;
        for (auto it: diffEquationToSolve->getInitialConditions())
            conditions += diffEquationToSolve->getDiffVariable()+it.display(rowNumber++)+", ";
        conditions.chop(2);
        if (!diffEquationToSolve->applyInitialConditions())
            conditions += " (not applied)";
    }
    else
        conditions += "not set";

    diffEquationToSolve->saveStep(conditions);

    // call parseEquation function
    if (!diffEquationToSolve->parseEquation(findDiffVarRegex)) {

        // if parsing fails, indicate error and return
        ui->lineEdit_EquationType->setText("parsing error");
        ui->lineEdit_EquationLeftSide->setStyleSheet(containsErrorBackground);
        ui->textEdit_Result->setPlainText("expression could not be parsed");
        return;
    }
    diffEquationToSolve->saveStep(QStringLiteral("equation successfully parsed"));

    // load tokens
    QString allTokens;
    for (auto it = diffEquationToSolve->getTokens().begin();
              it != diffEquationToSolve->getTokens().end(); ++it) {

        allTokens += it->getToken() + ", ";
    }
    allTokens.chop(2);
    diffEquationToSolve->saveStep(QStringLiteral("tokens: ")+allTokens);

    // determine type of equation
    bool canBeSolved = true;
    const QString equationType = diffEquationToSolve->determineType(canBeSolved);
    ui->lineEdit_EquationType->setText(equationType);
    diffEquationToSolve->saveStep(QStringLiteral("type: ")+equationType);

    // default configuration of fundamental matrix
    fundamentalMatrix.first = false;
    fundamentalMatrix.second.empty();

    // if equation could be (according to its type evaluation) solved
    if (!canBeSolved || // call this procedure to do the calculations
        !diffEquationToSolve->generalSolution(fundamentalMatrix)) {

        // description for equations that can't be solved (by this software)
        diffEquationToSolve->saveStep(noSolution);
        ui->textEdit_Result->setPlainText(noSolution);
        return;
    }

    // show result
    QString displaySolution =
        QStringLiteral("<span style=""font-size:14pt; font-weight:bold;"">") +
        diffEquationToSolve->getDiffVariable() + " = " +
        diffEquationToSolve->getFinalResult() + "</span>";
    ui->textEdit_Result->setHtml(displaySolution);

    return;
}

void MainWindow::clearEquation_clicked() {

    // remove initial conditions
    ui->pushButton_RemoveConditions->click();
    // reset cursor position
    cursorPosition = 0;
    // clear all elements
    ui->lineEdit_EquationLeftSide->clear();
    ui->lineEdit_EquationRightSide->clear();
    ui->lineEdit_EquationType->clear();
    ui->textEdit_Result->clear();

    return;
}

int MainWindow::showSteps_clicked() {

    // extract VWADWRT
    const QChar diffVariable = ui->lineEdit_DiffSign->text().at(1);

    QString displaySolution;
    for (auto it : diffEquationToSolve->getSolution())
        // replace "<" and ">" characters with escape sequences (HTML enforced)
        displaySolution +=
            it.replace(" < "," &lt; ").replace(" > "," &gt; ") + "<br>";

    StepsWindow solutionInStepsWindow(diffVariable, displaySolution,
                                      fundamentalMatrix, this);
    return solutionInStepsWindow.exec();
}

int MainWindow::action_About() {

    AboutWindow aboutWindow(this);
    return aboutWindow.exec();
}
