#include "settingsdialog.h"
#include <QHBoxLayout>
#include <QLabel>
#include <QVBoxLayout>

SettingsDialog::SettingsDialog(QWidget* parent) : QDialog(parent)
{
    setWindowTitle("Settings");
    setModal(true);

    QFormLayout* formLayout = new QFormLayout();

    // Temperature
    temperatureSpin = new QDoubleSpinBox();
    temperatureSpin->setRange(0.0, 10000.0);
    temperatureSpin->setValue(298.15);
    temperatureSpin->setSuffix(" K");
    formLayout->addRow("Temperature (T):", temperatureSpin);

    // Pressure
    pressureSpin = new QDoubleSpinBox();
    pressureSpin->setRange(0.0, 1000.0);
    pressureSpin->setValue(1.0);
    pressureSpin->setSuffix(" atm");
    formLayout->addRow("Pressure (P):", pressureSpin);

    // Concentration
    concentrationSpin = new QDoubleSpinBox();
    concentrationSpin->setRange(0.0, 100.0);
    concentrationSpin->setValue(0.0);
    concentrationSpin->setSuffix(" mol/L");
    formLayout->addRow("Concentration:", concentrationSpin);

    // Low vibration method
    lowVibCombo = new QComboBox();
    lowVibCombo->addItem("Harmonic");
    lowVibCombo->addItem("Truhlar");
    lowVibCombo->addItem("Grimme");
    lowVibCombo->addItem("Minenkov");
    lowVibCombo->setCurrentText("Truhlar");
    formLayout->addRow("Low vib method:", lowVibCombo);

    // Imaginary frequency threshold
    imagRealSpin = new QDoubleSpinBox();
    imagRealSpin->setRange(0.0, 100.0);
    imagRealSpin->setValue(0.0);
    imagRealSpin->setSuffix(" cm⁻¹");
    formLayout->addRow("Imag freq threshold:", imagRealSpin);

    // Buttons
    buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    QVBoxLayout* mainLayout = new QVBoxLayout();
    mainLayout->addLayout(formLayout);
    mainLayout->addWidget(buttonBox);

    setLayout(mainLayout);
}