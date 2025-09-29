#ifndef SETTINGSDIALOG_H
#define SETTINGSDIALOG_H

#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFormLayout>

class SettingsDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SettingsDialog(QWidget* parent = nullptr);

    double getTemperature() const
    {
        return temperatureSpin->value();
    }
    double getPressure() const
    {
        return pressureSpin->value();
    }
    double getConcentration() const
    {
        return concentrationSpin->value();
    }
    QString getLowVibMethod() const
    {
        return lowVibCombo->currentText();
    }
    double getImagReal() const
    {
        return imagRealSpin->value();
    }

    void setTemperature(double t)
    {
        temperatureSpin->setValue(t);
    }
    void setPressure(double p)
    {
        pressureSpin->setValue(p);
    }
    void setConcentration(double c)
    {
        concentrationSpin->setValue(c);
    }
    void setLowVibMethod(const QString& method)
    {
        lowVibCombo->setCurrentText(method);
    }
    void setImagReal(double ir)
    {
        imagRealSpin->setValue(ir);
    }

private:
    QDoubleSpinBox*   temperatureSpin;
    QDoubleSpinBox*   pressureSpin;
    QDoubleSpinBox*   concentrationSpin;
    QComboBox*        lowVibCombo;
    QDoubleSpinBox*   imagRealSpin;
    QDialogButtonBox* buttonBox;
};

#endif  // SETTINGSDIALOG_H