/*
 * 文件名: paramselectdialog.cpp
 * 文件作用: 参数选择配置对话框的具体实现
 * 功能描述:
 * 1. [核心修改] 调整列顺序: 显示->当前数值->单位->参数名称->拟合变量->下限->上限->滚轮步长。
 * 2. 保持滚轮步长6位小数精度。
 * 3. 保持滚轮拦截功能。
 */

#include "paramselectdialog.h"
#include "ui_paramselectdialog.h"
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QDebug>

ParamSelectDialog::ParamSelectDialog(const QList<FitParameter> &params, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ParamSelectDialog),
    m_params(params)
{
    ui->setupUi(this);
    this->setWindowTitle("拟合参数配置");
    connect(ui->btnOk, &QPushButton::clicked, this, &ParamSelectDialog::onConfirm);
    connect(ui->btnCancel, &QPushButton::clicked, this, &ParamSelectDialog::onCancel);
    ui->btnCancel->setAutoDefault(false);
    initTable();
}

ParamSelectDialog::~ParamSelectDialog()
{
    delete ui;
}

bool ParamSelectDialog::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::Wheel) {
        if (qobject_cast<QAbstractSpinBox*>(obj)) {
            return true;
        }
    }
    return QDialog::eventFilter(obj, event);
}

void ParamSelectDialog::initTable()
{
    QStringList headers;
    // [修改] 新的列顺序
    // 0:显示, 1:数值, 2:单位, 3:名称, 4:拟合变量, 5:下限, 6:上限, 7:步长
    headers << "显示" << "当前数值" << "单位" << "参数名称" << "拟合变量" << "下限" << "上限" << "滚轮步长";
    ui->tableWidget->setColumnCount(headers.size());
    ui->tableWidget->setHorizontalHeaderLabels(headers);
    ui->tableWidget->setRowCount(m_params.size());

    QString checkBoxStyle =
        "QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #cccccc; border-radius: 3px; background-color: white; }"
        "QCheckBox::indicator:checked { background-color: #0078d7; border-color: #0078d7; }"
        "QCheckBox::indicator:hover { border-color: #0078d7; }";

    for(int i = 0; i < m_params.size(); ++i) {
        const FitParameter& p = m_params[i];

        // Col 0: 显示 (Visible)
        QWidget* pWidgetVis = new QWidget();
        QHBoxLayout* pLayoutVis = new QHBoxLayout(pWidgetVis);
        QCheckBox* chkVis = new QCheckBox();
        chkVis->setChecked(p.isVisible);
        chkVis->setStyleSheet(checkBoxStyle);
        pLayoutVis->addWidget(chkVis);
        pLayoutVis->setAlignment(Qt::AlignCenter);
        pLayoutVis->setContentsMargins(0,0,0,0);
        ui->tableWidget->setCellWidget(i, 0, pWidgetVis);

        // Col 1: 当前数值 (Value)
        QDoubleSpinBox* spinVal = new QDoubleSpinBox();
        spinVal->setRange(-9e9, 9e9);
        spinVal->setDecimals(6);
        spinVal->setValue(p.value);
        spinVal->setFrame(false);
        spinVal->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 1, spinVal);

        // Col 2: 单位 (Unit)
        QString dummy, dummy2, dummy3, unitStr;
        FittingParameterChart::getParamDisplayInfo(p.name, dummy, dummy2, dummy3, unitStr);
        if(unitStr == "无因次" || unitStr == "小数") unitStr = "-";
        QTableWidgetItem* unitItem = new QTableWidgetItem(unitStr);
        unitItem->setFlags(unitItem->flags() & ~Qt::ItemIsEditable);
        ui->tableWidget->setItem(i, 2, unitItem);

        // Col 3: 参数名称 (Name)
        QString displayNameFull = QString("%1 (%2)").arg(p.displayName).arg(p.name);
        QTableWidgetItem* nameItem = new QTableWidgetItem(displayNameFull);
        nameItem->setFlags(nameItem->flags() & ~Qt::ItemIsEditable);
        nameItem->setData(Qt::UserRole, p.name);
        ui->tableWidget->setItem(i, 3, nameItem);

        // Col 4: 拟合变量 (Fit) - 允许选择哪些参数参与拟合
        QWidget* pWidgetFit = new QWidget();
        QHBoxLayout* pLayoutFit = new QHBoxLayout(pWidgetFit);
        QCheckBox* chkFit = new QCheckBox();
        chkFit->setChecked(p.isFit);
        chkFit->setStyleSheet(checkBoxStyle);
        pLayoutFit->addWidget(chkFit);
        pLayoutFit->setAlignment(Qt::AlignCenter);
        pLayoutFit->setContentsMargins(0,0,0,0);

        // 基础参数允许勾选，但默认未勾选；模型参数默认勾选
        // LfD 特殊处理：只读，不允许拟合
        if (p.name == "LfD") {
            chkFit->setEnabled(false);
            chkFit->setChecked(false);
        }
        ui->tableWidget->setCellWidget(i, 4, pWidgetFit);

        // 联动逻辑：如果勾选“拟合”，则强制“显示”
        connect(chkFit, &QCheckBox::checkStateChanged, [chkVis](Qt::CheckState state){
            if (state == Qt::Checked) {
                chkVis->setChecked(true);
                chkVis->setEnabled(false);
                chkVis->setStyleSheet("QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #ccc; border-radius: 3px; background-color: #e0e0e0; } "
                                      "QCheckBox::indicator:checked { background-color: #80bbeb; border-color: #80bbeb; }");
            } else {
                chkVis->setEnabled(true);
                chkVis->setStyleSheet(
                    "QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #cccccc; border-radius: 3px; background-color: white; }"
                    "QCheckBox::indicator:checked { background-color: #0078d7; border-color: #0078d7; }"
                    "QCheckBox::indicator:hover { border-color: #0078d7; }"
                    );
            }
        });

        // 初始化状态样式
        if (p.isFit) {
            chkVis->setChecked(true);
            chkVis->setEnabled(false);
            chkVis->setStyleSheet("QCheckBox::indicator { width: 20px; height: 20px; border: 1px solid #ccc; border-radius: 3px; background-color: #e0e0e0; } "
                                  "QCheckBox::indicator:checked { background-color: #80bbeb; border-color: #80bbeb; }");
        }

        // Col 5: 下限 (Min)
        QDoubleSpinBox* spinMin = new QDoubleSpinBox();
        spinMin->setRange(-9e9, 9e9);
        spinMin->setDecimals(6);
        spinMin->setValue(p.min);
        spinMin->setFrame(false);
        spinMin->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 5, spinMin);

        // Col 6: 上限 (Max)
        QDoubleSpinBox* spinMax = new QDoubleSpinBox();
        spinMax->setRange(-9e9, 9e9);
        spinMax->setDecimals(6);
        spinMax->setValue(p.max);
        spinMax->setFrame(false);
        spinMax->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 6, spinMax);

        // Col 7: 滚轮步长 (Step)
        QDoubleSpinBox* spinStep = new QDoubleSpinBox();
        spinStep->setRange(0.0, 10000.0);
        spinStep->setDecimals(6);
        spinStep->setValue(p.step);
        spinStep->setFrame(false);
        spinStep->installEventFilter(this);
        ui->tableWidget->setCellWidget(i, 7, spinStep);
    }

    ui->tableWidget->resizeColumnsToContents();
    // 调整列宽分配，让名称列自适应
    ui->tableWidget->horizontalHeader()->setSectionResizeMode(3, QHeaderView::Stretch); // Col 3 is Name
}

void ParamSelectDialog::collectData()
{
    for(int i = 0; i < ui->tableWidget->rowCount(); ++i) {
        if(i >= m_params.size()) break;

        // 0: Visible
        QWidget* wVis = ui->tableWidget->cellWidget(i, 0);
        if (wVis) {
            QCheckBox* chkVis = wVis->findChild<QCheckBox*>();
            if(chkVis) m_params[i].isVisible = chkVis->isChecked();
        }

        // 1: Value
        QDoubleSpinBox* spinVal = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 1));
        if(spinVal) m_params[i].value = spinVal->value();

        // 4: Fit
        QWidget* wFit = ui->tableWidget->cellWidget(i, 4);
        if (wFit) {
            QCheckBox* chkFit = wFit->findChild<QCheckBox*>();
            if(chkFit) m_params[i].isFit = chkFit->isChecked();
        }

        // 5: Min
        QDoubleSpinBox* spinMin = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 5));
        if(spinMin) m_params[i].min = spinMin->value();

        // 6: Max
        QDoubleSpinBox* spinMax = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 6));
        if(spinMax) m_params[i].max = spinMax->value();

        // 7: Step
        QDoubleSpinBox* spinStep = qobject_cast<QDoubleSpinBox*>(ui->tableWidget->cellWidget(i, 7));
        if(spinStep) m_params[i].step = spinStep->value();
    }
}

QList<FitParameter> ParamSelectDialog::getUpdatedParams() const
{
    return m_params;
}

void ParamSelectDialog::onConfirm()
{
    collectData();
    accept();
}

void ParamSelectDialog::onCancel()
{
    reject();
}
