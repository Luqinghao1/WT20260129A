/*
 * wt_modelwidget.h
 * 文件作用: 压裂水平井复合页岩油模型界面类头文件
 * 功能描述:
 * 1. 定义界面类，负责用户交互和结果显示。
 * 2. 声明界面初始化、信号连接和计算流程控制函数。
 * 3. 包含绘图相关逻辑。
 */

#ifndef WT_MODELWIDGET_H
#define WT_MODELWIDGET_H

#include <QWidget>
#include <QVector>
#include <QMap>
#include <tuple>
#include "modelsolver01-06.h"

namespace Ui {
class WT_ModelWidget;
}

class WT_ModelWidget : public QWidget
{
    Q_OBJECT

public:
    using ModelType = ModelSolver01_06::ModelType;
    // ModelCurveData: <时间, 压力, 导数>
    using ModelCurveData = std::tuple<QVector<double>, QVector<double>, QVector<double>>;

    explicit WT_ModelWidget(ModelType type, QWidget *parent = nullptr);
    ~WT_ModelWidget();

    // 外部设置/获取
    QString getModelName() const;
    void setHighPrecision(bool high);

signals:
    // 计算完成后发送信号，传递参数供外部保存
    void calculationCompleted(const QString& modelName, const QMap<QString, double>& params);
    // 请求切换模型
    void requestModelSelection();

private slots:
    void onCalculateClicked();
    void onResetParameters();
    void onExportData();
    void onShowPointsToggled(bool checked);
    void onDependentParamsChanged(); // 处理参数联动

private:
    void initUi();
    void initChart();
    void setupConnections();
    void runCalculation();

    // 辅助函数
    QVector<double> parseInput(const QString& text);
    void setInputText(class QLineEdit* edit, double value);
    void plotCurve(const ModelCurveData& data, const QString& name, QColor color, bool isSensitivity);
    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params, const QVector<double>& providedTime);

private:
    Ui::WT_ModelWidget *ui;
    ModelSolver01_06* m_solver;
    ModelType m_type;
    bool m_highPrecision;

    // 结果缓存
    QVector<double> res_tD;
    QVector<double> res_pD;
    QVector<double> res_dpD;

    QList<QColor> m_colorList;
};

#endif // WT_MODELWIDGET_H
