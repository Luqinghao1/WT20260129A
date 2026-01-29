/*
 * 文件名: fittingchart.h
 * 文件作用: 拟合绘图管理类头文件
 * 功能描述:
 * 1. 统一管理 LogLog(双对数), SemiLog(半对数), Cartesian(笛卡尔) 三个图表的绘制。
 * 2. 实现不同试井类型（降落/恢复）的差异化绘图逻辑。
 * 3. 针对恢复试井，在半对数图上绘制 Horner Plot 并计算初始地层压力。
 */

#ifndef FITTINGCHART_H
#define FITTINGCHART_H

#include <QObject>
#include <QVector>
#include "mousezoom.h"
#include "fittingdatadialog.h" // 引入 FittingDataSettings

class FittingChart : public QObject
{
    Q_OBJECT
public:
    explicit FittingChart(QObject *parent = nullptr);

    // 初始化图表指针和默认设置
    void initializeCharts(MouseZoom* logLog, MouseZoom* semiLog, MouseZoom* cartesian);

    // 更新观测数据
    // rawP: 原始实测压力（非压差）
    void setObservedData(const QVector<double>& t, const QVector<double>& deltaP,
                         const QVector<double>& deriv, const QVector<double>& rawP);

    // 更新设置（主要用于获取试井类型和tp）
    void setSettings(const FittingDataSettings& settings);

    // 绘制所有曲线（包含实测和理论模型）
    // t_model, p_model, d_model: 理论曲线数据
    void plotAll(const QVector<double>& t_model, const QVector<double>& p_model, const QVector<double>& d_model, bool isModelValid);

    // 绘制抽样点（调试用）
    void plotSampledPoints(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d);

    // 获取计算出的初始压力（仅恢复试井有效）
    double getCalculatedInitialPressure() const;

private:
    MouseZoom* m_plotLogLog;
    MouseZoom* m_plotSemiLog;
    MouseZoom* m_plotCartesian;

    // 观测数据
    QVector<double> m_obsT;
    QVector<double> m_obsDeltaP;
    QVector<double> m_obsDeriv;
    QVector<double> m_obsRawP; // 原始压力

    FittingDataSettings m_settings;
    double m_calculatedPi; // 计算出的初始地层压力

    // 内部绘图函数
    void plotLogLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel);
    void plotSemiLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel);
    void plotCartesian(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel);

    // 辅助：Horner Plot 计算与拟合
    // 返回计算出的 Pi
    double calculateHornerPressure();

    // 辅助：显示结果到双对数图
    void showResultOnLogPlot();
};

#endif // FITTINGCHART_H
