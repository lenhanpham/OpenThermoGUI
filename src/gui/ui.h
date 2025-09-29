#ifndef UI_H
#define UI_H

#include "../calc.h"
#include "../loadfile.h"
#include "settingsdialog.h"
#include <QAction>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QHBoxLayout>
#include <QLineEdit>
#include <QMainWindow>
#include <QMenuBar>

#include <QSortFilterProxyModel>
#include <QTableView>
#include <QTableWidget>
#include <QToolBar>
#include <QTreeWidget>
#include <QVBoxLayout>

class ThermoModel;

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow(QWidget* parent = nullptr);

private slots:
    void openFiles();
    void runAnalysis();
    void showSettings();
    void showAbout();
    void exportToCSV();
    void updateParameter(int row, int column);
    void toggleFileList();
    void toggleParameters();
    void filterResults(const QString& text);

private:
    // Menu and toolbar
    QMenuBar* menuBar;
    QToolBar* toolBar;
    QAction*  openAction;
    QAction*  runAction;
    QAction*  settingsAction;
    QAction*  exportAction;
    QAction*  exitAction;
    QAction*  aboutAction;
    QAction*  showFileListAction;
    QAction*  showParametersAction;

    // Main widgets
    QTreeWidget*           fileTree;
    QTableView*            resultsTable;
    QLineEdit*             searchBar;
    ThermoModel*           model;
    QSortFilterProxyModel* proxyModel;

    // Parameters panel
    QTableWidget* parametersTable;

    // Settings
    double  temperature   = 298.15;
    double  pressure      = 1.0;
    double  concentration = 0.0;
    QString lowVibMethod  = "Truhlar";
    double  imagReal      = 0.0;
    double  sclZPE        = 1.0;
    double  sclHeat       = 1.0;
    double  sclEntropy    = 1.0;
    double  sclCV         = 1.0;
    double  ravib         = 100.0;
    double  intpvib       = 100.0;
    int     ipmode        = 0;
    int     massmod       = 3;
    int     prtvib        = 0;
    int     outotm        = 0;

    // Loaded systems
    std::vector<SystemData> systems;
};

#endif