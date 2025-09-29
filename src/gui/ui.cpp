#include "ui.h"
#include "../symmetry.h"
#include "../util.h"
#include <QEvent>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QIcon>
#include <QLabel>
#include <QLineEdit>
#include <QMenuBar>
#include <QMessageBox>
#include <QMouseEvent>
#include <QPushButton>
#include <QSplitter>
#include <QStandardItem>
#include <QStandardItemModel>
#include <QTableView>
#include <QTableWidgetItem>
#include <QTextStream>
#include <QToolBar>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QVBoxLayout>

class ThermoModel : public QStandardItemModel
{
public:
    explicit ThermoModel(QObject* parent = nullptr) : QStandardItemModel(parent)
    {
        setColumnCount(11);
        setHorizontalHeaderLabels({"File Name",
                                   "T (K)",
                                   "P (atm)",
                                   "G (kJ/mol)",
                                   "H (kJ/mol)",
                                   "S (kJ/mol/K)",
                                   "ZPE (kJ/mol)",
                                   "Thermal U (kJ/mol)",
                                   "Thermal H (kJ/mol)",
                                   "Thermal G (kJ/mol)",
                                   "Electronic E (a.u.)"});
    }

    void addResult(const QString& file,
                   double         T,
                   double         P,
                   double         G,
                   double         H,
                   double         S,
                   double         ZPE,
                   double         corrU,
                   double         corrH,
                   double         corrG,
                   double         E)
    {
        QList<QStandardItem*> row;
        row << new QStandardItem(file) << new QStandardItem(QString::number(T)) << new QStandardItem(QString::number(P))
            << new QStandardItem(QString::number(G)) << new QStandardItem(QString::number(H))
            << new QStandardItem(QString::number(S)) << new QStandardItem(QString::number(ZPE))
            << new QStandardItem(QString::number(corrU)) << new QStandardItem(QString::number(corrH))
            << new QStandardItem(QString::number(corrG)) << new QStandardItem(QString::number(E));
        appendRow(row);
        files << file;
    }

    QStringList getFiles() const
    {
        return files;
    }

    void clear()
    {
        QStandardItemModel::clear();
        files.clear();
        setColumnCount(11);
        setHorizontalHeaderLabels({"File Name",
                                   "T (K)",
                                   "P (atm)",
                                   "G (kJ/mol)",
                                   "H (kJ/mol)",
                                   "S (kJ/mol/K)",
                                   "ZPE (kJ/mol)",
                                   "Thermal U (kJ/mol)",
                                   "Thermal H (kJ/mol)",
                                   "Thermal G (kJ/mol)",
                                   "Electronic E (a.u.)"});
    }

private:
    QStringList files;
};

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent)
{
    setWindowTitle("OpenThermo");
    setWindowIcon(QIcon(":/icon.png"));
    resize(1000, 700);

    // Create menu bar
    menuBar = new QMenuBar(this);
    setMenuBar(menuBar);

    // File menu
    QMenu* fileMenu = menuBar->addMenu("&File");
    openAction      = fileMenu->addAction("&Open Files...");
    openAction->setShortcut(QKeySequence::Open);
    openAction->setToolTip("Open quantum chemistry output files for analysis");
    openAction->setIcon(QIcon(":/icons/document-import.svg"));
    connect(openAction, &QAction::triggered, this, &MainWindow::openFiles);

    exportAction = fileMenu->addAction("&Export Results...");
    exportAction->setShortcut(QKeySequence("Ctrl+E"));
    exportAction->setToolTip("Export analysis results to CSV file");
    exportAction->setIcon(QIcon(":/icons/document-export.svg"));
    connect(exportAction, &QAction::triggered, this, &MainWindow::exportToCSV);

    fileMenu->addSeparator();
    exitAction = fileMenu->addAction("E&xit");
    exitAction->setShortcut(QKeySequence::Quit);
    exitAction->setToolTip("Exit the application");
    exitAction->setIcon(QIcon(":/icons/application-exit.svg"));
    connect(exitAction, &QAction::triggered, this, &QWidget::close);

    // Edit menu
    QMenu* editMenu = menuBar->addMenu("&Edit");
    settingsAction  = editMenu->addAction("&Preferences...");
    settingsAction->setShortcut(QKeySequence("Ctrl+,"));
    settingsAction->setToolTip("Configure analysis settings");

    connect(settingsAction, &QAction::triggered, this, &MainWindow::showSettings);

    // View menu
    QMenu* viewMenu    = menuBar->addMenu("&View");
    showFileListAction = viewMenu->addAction("Show &File List");
    showFileListAction->setCheckable(true);
    showFileListAction->setChecked(true);
    showFileListAction->setShortcut(QKeySequence("Ctrl+L"));
    showFileListAction->setToolTip("Show/hide the file list panel");
    connect(showFileListAction, &QAction::toggled, this, &MainWindow::toggleFileList);

    showParametersAction = viewMenu->addAction("Show &Parameters");
    showParametersAction->setCheckable(true);
    showParametersAction->setChecked(true);
    showParametersAction->setShortcut(QKeySequence("Ctrl+P"));
    showParametersAction->setToolTip("Show/hide the parameters panel");
    connect(showParametersAction, &QAction::toggled, this, &MainWindow::toggleParameters);

    // Analysis menu
    QMenu* analysisMenu = menuBar->addMenu("&Analysis");
    runAction           = analysisMenu->addAction("&Run Analysis");
    runAction->setShortcut(QKeySequence("Ctrl+R"));
    runAction->setToolTip("Run thermodynamic analysis on loaded files");
    runAction->setIcon(QIcon(":/icons/run.svg"));
    connect(runAction, &QAction::triggered, this, &MainWindow::runAnalysis);

    // Help menu
    QMenu* helpMenu = menuBar->addMenu("&Help");
    aboutAction     = helpMenu->addAction("&About");
    aboutAction->setToolTip("About OpenThermo");
    connect(aboutAction, &QAction::triggered, this, &MainWindow::showAbout);

    // Create toolbar
    toolBar = addToolBar("Main");
    toolBar->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    toolBar->addAction(openAction);
    toolBar->addAction(runAction);
    toolBar->addAction(exportAction);
    toolBar->addSeparator();
    toolBar->addAction(settingsAction);

    // Create central widget
    QWidget* centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);

    // Create splitter for main layout
    QSplitter* splitter = new QSplitter(Qt::Horizontal, centralWidget);
    splitter->setHandleWidth(1);  // Set uniform thin separator width
    splitter->setStyleSheet("QSplitter::handle { background-color: #d0d0d0; }"
                            "QGroupBox { border: 1px solid #cccccc; border-radius: 3px; margin-top: 20px; }"
                            "QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top center; padding: 0 "
                            "5px; font-weight: bold; }");

    // Left Panel (File List)
    QGroupBox*   leftPanel  = new QGroupBox("Loaded Files");
    QVBoxLayout* leftLayout = new QVBoxLayout(leftPanel);
    leftLayout->setContentsMargins(5, 5, 5, 5);

    fileTree = new QTreeWidget();
    fileTree->setHeaderHidden(true);  // Hide header since title is in group box
    fileTree->setToolTip("List of loaded quantum chemistry output files");
    leftLayout->addWidget(fileTree);

    splitter->addWidget(leftPanel);

    // Central Area
    QWidget*     mainArea   = new QWidget();
    QVBoxLayout* mainLayout = new QVBoxLayout(mainArea);
    searchBar               = new QLineEdit();
    searchBar->setPlaceholderText("Search results...");
    searchBar->setToolTip("Filter the results table");
    mainLayout->addWidget(searchBar);
    resultsTable = new QTableView();
    model        = new ThermoModel(this);
    proxyModel   = new QSortFilterProxyModel(this);
    proxyModel->setSourceModel(model);
    proxyModel->setFilterKeyColumn(0);  // Filter by file name column
    resultsTable->setModel(proxyModel);
    resultsTable->horizontalHeader()->setStretchLastSection(true);
    resultsTable->setToolTip("Thermodynamic analysis results");
    mainLayout->addWidget(resultsTable);
    splitter->addWidget(mainArea);

    // Right Panel (Parameters)
    QGroupBox*   rightPanel  = new QGroupBox("Current Parameters");
    QVBoxLayout* rightLayout = new QVBoxLayout(rightPanel);
    rightLayout->setContentsMargins(5, 5, 5, 5);

    parametersTable = new QTableWidget();
    parametersTable->setColumnCount(2);
    parametersTable->setHorizontalHeaderLabels({"Parameter", "Value"});
    parametersTable->horizontalHeader()->setStretchLastSection(true);
    parametersTable->setRowCount(15);
    parametersTable->setAlternatingRowColors(true);
    parametersTable->setSelectionBehavior(QAbstractItemView::SelectRows);
    parametersTable->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::SelectedClicked);

    // Temperature
    parametersTable->setItem(0, 0, new QTableWidgetItem("Temperature (K)"));
    QTableWidgetItem* tempItem = new QTableWidgetItem(QString::number(temperature));
    tempItem->setFlags(tempItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(0, 1, tempItem);

    // Pressure
    parametersTable->setItem(1, 0, new QTableWidgetItem("Pressure (atm)"));
    QTableWidgetItem* pressItem = new QTableWidgetItem(QString::number(pressure));
    pressItem->setFlags(pressItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(1, 1, pressItem);

    // Concentration
    parametersTable->setItem(2, 0, new QTableWidgetItem("Concentration (mol/L)"));
    QTableWidgetItem* concItem = new QTableWidgetItem(QString::number(concentration));
    concItem->setFlags(concItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(2, 1, concItem);

    // Low vibration method
    parametersTable->setItem(3, 0, new QTableWidgetItem("Low vib method"));
    QTableWidgetItem* vibItem = new QTableWidgetItem(lowVibMethod);
    vibItem->setFlags(vibItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(3, 1, vibItem);

    // Imaginary frequency threshold
    parametersTable->setItem(4, 0, new QTableWidgetItem("Imag freq threshold (cm⁻¹)"));
    QTableWidgetItem* imagItem = new QTableWidgetItem(QString::number(imagReal));
    imagItem->setFlags(imagItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(4, 1, imagItem);

    // ZPE scaling factor
    parametersTable->setItem(5, 0, new QTableWidgetItem("ZPE scaling factor"));
    QTableWidgetItem* zpeItem = new QTableWidgetItem(QString::number(sclZPE));
    zpeItem->setFlags(zpeItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(5, 1, zpeItem);

    // Heat scaling factor
    parametersTable->setItem(6, 0, new QTableWidgetItem("Heat scaling factor"));
    QTableWidgetItem* heatItem = new QTableWidgetItem(QString::number(sclHeat));
    heatItem->setFlags(heatItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(6, 1, heatItem);

    // Entropy scaling factor
    parametersTable->setItem(7, 0, new QTableWidgetItem("Entropy scaling factor"));
    QTableWidgetItem* entropyItem = new QTableWidgetItem(QString::number(sclEntropy));
    entropyItem->setFlags(entropyItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(7, 1, entropyItem);

    // CV scaling factor
    parametersTable->setItem(8, 0, new QTableWidgetItem("CV scaling factor"));
    QTableWidgetItem* cvItem = new QTableWidgetItem(QString::number(sclCV));
    cvItem->setFlags(cvItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(8, 1, cvItem);

    // Raising value for low frequencies
    parametersTable->setItem(9, 0, new QTableWidgetItem("Low freq raising value (cm⁻¹)"));
    QTableWidgetItem* ravibItem = new QTableWidgetItem(QString::number(ravib));
    ravibItem->setFlags(ravibItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(9, 1, ravibItem);

    // Interpolation parameter
    parametersTable->setItem(10, 0, new QTableWidgetItem("Interpolation parameter"));
    QTableWidgetItem* intpvibItem = new QTableWidgetItem(QString::number(intpvib));
    intpvibItem->setFlags(intpvibItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(10, 1, intpvibItem);

    // Calculation mode
    parametersTable->setItem(11, 0, new QTableWidgetItem("Calculation mode"));
    QTableWidgetItem* ipmodeItem = new QTableWidgetItem(QString::number(ipmode));
    ipmodeItem->setFlags(ipmodeItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(11, 1, ipmodeItem);

    // Mass assignment mode
    parametersTable->setItem(12, 0, new QTableWidgetItem("Mass assignment mode"));
    QTableWidgetItem* massmodItem = new QTableWidgetItem(QString::number(massmod));
    massmodItem->setFlags(massmodItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(12, 1, massmodItem);

    // Print vibration contributions
    parametersTable->setItem(13, 0, new QTableWidgetItem("Print vib contributions"));
    QTableWidgetItem* prtvibItem = new QTableWidgetItem(QString::number(prtvib));
    prtvibItem->setFlags(prtvibItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(13, 1, prtvibItem);

    // Output .otm file
    parametersTable->setItem(14, 0, new QTableWidgetItem("Output .otm file"));
    QTableWidgetItem* outotmItem = new QTableWidgetItem(QString::number(outotm));
    outotmItem->setFlags(outotmItem->flags() | Qt::ItemIsEditable);
    parametersTable->setItem(14, 1, outotmItem);

    rightLayout->addWidget(parametersTable);
    splitter->addWidget(rightPanel);

    // Connect parameter changes
    connect(parametersTable, &QTableWidget::cellChanged, this, &MainWindow::updateParameter);

    // Set size policies for consistent panel heights
    leftPanel->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
    rightPanel->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);

    // Set stretch factors for proper resizing
    splitter->setStretchFactor(0, 1);  // Left panel
    splitter->setStretchFactor(1, 3);  // Middle panel (main content)
    splitter->setStretchFactor(2, 1);  // Right panel

    splitter->setSizes({200, 500, 200});

    // Set central layout
    QVBoxLayout* centralLayout = new QVBoxLayout(centralWidget);
    centralLayout->addWidget(splitter);
    centralLayout->setContentsMargins(5, 5, 5, 5);

    // Connect signals
    connect(searchBar, &QLineEdit::textChanged, this, &MainWindow::filterResults);
}

void MainWindow::openFiles()
{
    QStringList files = QFileDialog::getOpenFileNames(
        this, "Select Quantum Output Files", "", "Output Files (*.log *.out *.otm *.txt);;All Files (*)");
    for (const QString& file : files)
    {
        try
        {
            SystemData sys;
            sys.inputfile = file.toStdString();
            // Determine program and load
            if (sys.inputfile.find(".otm") != std::string::npos)
            {
                LoadFile::loadotm(sys);
            }
            else
            {
                auto pcprog = util::deterprog(sys);
                if (pcprog == util::QuantumChemistryProgram::Gaussian)
                {
                    LoadFile::loadgau(sys);
                }
                else if (pcprog == util::QuantumChemistryProgram::Orca)
                {
                    LoadFile::loadorca(sys);
                }
                else if (pcprog == util::QuantumChemistryProgram::Gamess)
                {
                    LoadFile::loadgms(sys);
                }
                else if (pcprog == util::QuantumChemistryProgram::Nwchem)
                {
                    LoadFile::loadnw(sys);
                }
                else if (pcprog == util::QuantumChemistryProgram::Cp2k)
                {
                    LoadFile::loadCP2K(sys);
                }
                else if (pcprog == util::QuantumChemistryProgram::Xtb)
                {
                    LoadFile::loadxtb(sys);
                }
                else if (pcprog == util::QuantumChemistryProgram::Vasp)
                {
                    LoadFile::loadvasp(sys);
                }
                else
                {
                    QMessageBox::warning(this, "Error", QString("Unknown file format for %1").arg(file));
                    continue;
                }
            }
            // Set default settings
            sys.T       = temperature;
            sys.P       = pressure;
            sys.concstr = std::to_string(concentration);
            // Set low vib method
            if (lowVibMethod == "Harmonic")
            {
                sys.lowVibTreatment = LowVibTreatment::Harmonic;
            }
            else if (lowVibMethod == "Truhlar")
            {
                sys.lowVibTreatment = LowVibTreatment::Truhlar;
            }
            else if (lowVibMethod == "Grimme")
            {
                sys.lowVibTreatment = LowVibTreatment::Grimme;
            }
            else if (lowVibMethod == "Minenkov")
            {
                sys.lowVibTreatment = LowVibTreatment::Minenkov;
            }
            sys.imagreal = imagReal;
            sys.sclZPE   = sclZPE;
            sys.sclheat  = sclHeat;
            sys.sclS     = sclEntropy;
            sys.sclCV    = sclCV;
            sys.ravib    = ravib;
            sys.intpvib  = intpvib;
            sys.ipmode   = ipmode;
            sys.massmod  = massmod;
            sys.prtvib   = prtvib;
            sys.outotm   = outotm;

            // Set electron levels (similar to ensemble mode) if not already loaded
            if (sys.nelevel == 0)
            {
                sys.nelevel = 1;
                sys.elevel  = {0.0};
                int deg     = std::max(sys.spinmult, 1);
                sys.edegen  = {deg};
            }

            // Calculate total mass
            sys.totmass = 0.0;
            for (const auto& atom : sys.a)
            {
                sys.totmass += atom.mass;
            }

            systems.push_back(sys);

            QTreeWidgetItem* item     = new QTreeWidgetItem(fileTree);
            QString          fileName = QFileInfo(file).fileName();
            item->setText(0, QString("%1. %2").arg(fileTree->topLevelItemCount()).arg(fileName));
            item->setData(0, Qt::UserRole, file);
        }
        catch (const std::exception& e)
        {
            QMessageBox::warning(this, "Error", QString("Failed to load %1: %2").arg(file, e.what()));
        }
    }
}

void MainWindow::runAnalysis()
{
    if (systems.empty())
    {
        QMessageBox::warning(this, "Error", "No files loaded!");
        return;
    }

    model->clear();  // Clear previous results

    for (size_t i = 0; i < systems.size(); ++i)
    {
        SystemData& sys = systems[i];
        // Update settings
        sys.T       = temperature;
        sys.P       = pressure;
        sys.concstr = std::to_string(concentration);
        if (lowVibMethod == "Harmonic")
        {
            sys.lowVibTreatment = LowVibTreatment::Harmonic;
        }
        else if (lowVibMethod == "Truhlar")
        {
            sys.lowVibTreatment = LowVibTreatment::Truhlar;
        }
        else if (lowVibMethod == "Grimme")
        {
            sys.lowVibTreatment = LowVibTreatment::Grimme;
        }
        else if (lowVibMethod == "Minenkov")
        {
            sys.lowVibTreatment = LowVibTreatment::Minenkov;
        }
        sys.imagreal = imagReal;
        sys.sclZPE   = sclZPE;
        sys.sclheat  = sclHeat;
        sys.sclS     = sclEntropy;
        sys.sclCV    = sclCV;
        sys.ravib    = ravib;
        sys.intpvib  = intpvib;
        sys.ipmode   = ipmode;
        sys.massmod  = massmod;
        sys.prtvib   = prtvib;
        sys.outotm   = outotm;

        try
        {
            // Calculate inertia
            calc::calcinertia(sys);

            // Symmetry
            symmetry::SymmetryDetector symDetector;
            symDetector.ncenter = sys.a.size();
            symDetector.a       = sys.a;
            symDetector.a_index.resize(sys.a.size());
            for (size_t j = 0; j < sys.a.size(); ++j)
            {
                symDetector.a_index[j] = j;
            }
            symDetector.detectPG(sys.prtvib ? 1 : 0);
            sys.rotsym  = symDetector.rotsym;
            sys.PGlabel = symDetector.PGlabel;

            // Handle imaginary frequencies
            if (sys.imagreal != 0.0)
            {
                for (int j = 0; j < sys.nfreq; ++j)
                {
                    if (sys.wavenum[j] < 0 && std::abs(sys.wavenum[j]) < sys.imagreal)
                    {
                        sys.wavenum[j] = std::abs(sys.wavenum[j]);
                    }
                }
            }

            // Convert wavenum to freq
            std::vector<double> freq(sys.nfreq);
            for (int j = 0; j < sys.nfreq; ++j)
            {
                freq[j] = sys.wavenum[j] * wave2freq;
            }
            sys.freq = freq;

            // Calculate thermo
            double corrU, corrH, corrG, S, CV, CP, QV, Qbot, ZPE;
            calc::calcthermo(sys, temperature, pressure, corrU, corrH, corrG, S, CV, CP, QV, Qbot, ZPE);

            // Final values
            double G_final = sys.E + corrG / au2kJ_mol;
            double H_final = sys.E + corrH / au2kJ_mol;
            double S_final = S / 1000.0;  // Convert to kJ/mol/K

            QString fileName = QFileInfo(QString::fromStdString(sys.inputfile)).fileName();
            model->addResult(
                fileName, temperature, pressure, G_final, H_final, S_final, ZPE, corrU, corrH, corrG, sys.E);
        }
        catch (const std::exception& e)
        {
            QMessageBox::warning(
                this,
                "Error",
                QString("Failed to analyze %1: %2").arg(QString::fromStdString(sys.inputfile), e.what()));
        }
    }
}

void MainWindow::showSettings()
{
    SettingsDialog dialog(this);
    dialog.setTemperature(temperature);
    dialog.setPressure(pressure);
    dialog.setConcentration(concentration);
    dialog.setLowVibMethod(lowVibMethod);
    dialog.setImagReal(imagReal);

    if (dialog.exec() == QDialog::Accepted)
    {
        temperature   = dialog.getTemperature();
        pressure      = dialog.getPressure();
        concentration = dialog.getConcentration();
        lowVibMethod  = dialog.getLowVibMethod();
        imagReal      = dialog.getImagReal();
    }
}

void MainWindow::showAbout()
{
    QMessageBox::about(this,
                       "About OpenThermo",
                       "OpenThermo - Thermodynamic Analysis Tool\n\n"
                       "A tool for analyzing quantum chemistry output files\n"
                       "and calculating thermodynamic properties.\n\n"
                       "Version 0.001.1");
}

void MainWindow::exportToCSV()
{
    if (proxyModel->rowCount() == 0)
    {
        QMessageBox::information(this, "Export Results", "No results to export. Please run analysis first.");
        return;
    }

    QString fileName =
        QFileDialog::getSaveFileName(this, "Export Results to CSV", "", "CSV Files (*.csv);;All Files (*)");
    if (fileName.isEmpty())
        return;

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QMessageBox::critical(this, "Export Error", "Could not open file for writing.");
        return;
    }

    QTextStream out(&file);

    // Write header
    for (int col = 0; col < proxyModel->columnCount(); ++col)
    {
        if (col > 0)
            out << ",";
        out << "\"" << proxyModel->headerData(col, Qt::Horizontal).toString() << "\"";
    }
    out << "\n";

    // Write data
    for (int row = 0; row < model->rowCount(); ++row)
    {
        for (int col = 0; col < model->columnCount(); ++col)
        {
            if (col > 0)
                out << ",";
            QString data = model->data(model->index(row, col)).toString();
            // Escape quotes and wrap in quotes if contains comma or quote
            if (data.contains(",") || data.contains("\"") || data.contains("\n"))
            {
                data.replace("\"", "\"\"");
                out << "\"" << data << "\"";
            }
            else
            {
                out << data;
            }
        }
        out << "\n";
    }

    file.close();
    QMessageBox::information(this, "Export Complete", "Results exported successfully to " + fileName);
}

void MainWindow::updateParameter(int row, int column)
{
    if (column != 1)
        return;  // Only process changes to the value column

    QTableWidgetItem* item = parametersTable->item(row, column);
    if (!item)
        return;

    QString value = item->text();

    switch (row)
    {
        case 0: {  // Temperature
            bool   ok;
            double newTemp = value.toDouble(&ok);
            if (ok && newTemp > 0)
            {
                temperature = newTemp;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new temperature
                }
            }
            else
            {
                item->setText(QString::number(temperature));  // Reset to valid value
            }
            break;
        }
        case 1: {  // Pressure
            bool   ok;
            double newPress = value.toDouble(&ok);
            if (ok && newPress >= 0)
            {
                pressure = newPress;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new pressure
                }
            }
            else
            {
                item->setText(QString::number(pressure));  // Reset to valid value
            }
            break;
        }
        case 2: {  // Concentration
            bool   ok;
            double newConc = value.toDouble(&ok);
            if (ok && newConc >= 0)
            {
                concentration = newConc;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new concentration
                }
            }
            else
            {
                item->setText(QString::number(concentration));  // Reset to valid value
            }
            break;
        }
        case 3: {  // Low vibration method
            QStringList validMethods = {"Harmonic", "Truhlar", "Grimme", "Minenkov"};
            if (validMethods.contains(value))
            {
                lowVibMethod = value;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new method
                }
            }
            else
            {
                item->setText(lowVibMethod);  // Reset to valid value
            }
            break;
        }
        case 4: {  // Imaginary frequency threshold
            bool   ok;
            double newImag = value.toDouble(&ok);
            if (ok && newImag >= 0)
            {
                imagReal = newImag;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new threshold
                }
            }
            else
            {
                item->setText(QString::number(imagReal));  // Reset to valid value
            }
            break;
        }
        case 5: {  // ZPE scaling factor
            bool   ok;
            double newSclZPE = value.toDouble(&ok);
            if (ok && newSclZPE > 0)
            {
                sclZPE = newSclZPE;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new scaling
                }
            }
            else
            {
                item->setText(QString::number(sclZPE));  // Reset to valid value
            }
            break;
        }
        case 6: {  // Heat scaling factor
            bool   ok;
            double newSclHeat = value.toDouble(&ok);
            if (ok && newSclHeat > 0)
            {
                sclHeat = newSclHeat;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new scaling
                }
            }
            else
            {
                item->setText(QString::number(sclHeat));  // Reset to valid value
            }
            break;
        }
        case 7: {  // Entropy scaling factor
            bool   ok;
            double newSclEntropy = value.toDouble(&ok);
            if (ok && newSclEntropy > 0)
            {
                sclEntropy = newSclEntropy;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new scaling
                }
            }
            else
            {
                item->setText(QString::number(sclEntropy));  // Reset to valid value
            }
            break;
        }
        case 8: {  // CV scaling factor
            bool   ok;
            double newSclCV = value.toDouble(&ok);
            if (ok && newSclCV > 0)
            {
                sclCV = newSclCV;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new scaling
                }
            }
            else
            {
                item->setText(QString::number(sclCV));  // Reset to valid value
            }
            break;
        }
        case 9: {  // Low freq raising value
            bool   ok;
            double newRavib = value.toDouble(&ok);
            if (ok && newRavib >= 0)
            {
                ravib = newRavib;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new value
                }
            }
            else
            {
                item->setText(QString::number(ravib));  // Reset to valid value
            }
            break;
        }
        case 10: {  // Interpolation parameter
            bool   ok;
            double newIntpvib = value.toDouble(&ok);
            if (ok && newIntpvib >= 0)
            {
                intpvib = newIntpvib;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new value
                }
            }
            else
            {
                item->setText(QString::number(intpvib));  // Reset to valid value
            }
            break;
        }
        case 11: {  // Calculation mode
            bool ok;
            int  newIpmode = value.toInt(&ok);
            if (ok && (newIpmode == 0 || newIpmode == 1))
            {
                ipmode = newIpmode;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new mode
                }
            }
            else
            {
                item->setText(QString::number(ipmode));  // Reset to valid value
            }
            break;
        }
        case 12: {  // Mass assignment mode
            bool ok;
            int  newMassmod = value.toInt(&ok);
            if (ok && newMassmod >= 0 && newMassmod <= 3)
            {
                massmod = newMassmod;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new mode
                }
            }
            else
            {
                item->setText(QString::number(massmod));  // Reset to valid value
            }
            break;
        }
        case 13: {  // Print vib contributions
            bool ok;
            int  newPrtvib = value.toInt(&ok);
            if (ok && (newPrtvib == 0 || newPrtvib == 1))
            {
                prtvib = newPrtvib;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new setting
                }
            }
            else
            {
                item->setText(QString::number(prtvib));  // Reset to valid value
            }
            break;
        }
        case 14: {  // Output .otm file
            bool ok;
            int  newOutotm = value.toInt(&ok);
            if (ok && (newOutotm == 0 || newOutotm == 1))
            {
                outotm = newOutotm;
                if (!systems.empty())
                {
                    runAnalysis();  // Re-run analysis with new setting
                }
            }
            else
            {
                item->setText(QString::number(outotm));  // Reset to valid value
            }
            break;
        }
    }
}

void MainWindow::toggleFileList()
{
    fileTree->setVisible(showFileListAction->isChecked());
}

void MainWindow::toggleParameters()
{
    // Find the parameters pane in the splitter
    QSplitter* splitter = centralWidget()->findChild<QSplitter*>();
    if (splitter)
    {
        for (int i = 0; i < splitter->count(); ++i)
        {
            QWidget* widget = splitter->widget(i);
            if (widget != fileTree && widget != resultsTable->parentWidget())
            {
                // This should be the parameters pane
                widget->setVisible(showParametersAction->isChecked());
                break;
            }
        }
    }
}

void MainWindow::filterResults(const QString& text)
{
    proxyModel->setFilterFixedString(text);
    proxyModel->setFilterCaseSensitivity(Qt::CaseInsensitive);
}
