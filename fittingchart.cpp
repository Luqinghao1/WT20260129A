/*
 * 文件名: fittingchart.cpp
 * 文件作用: 拟合绘图管理类实现文件
 * 功能描述:
 * 1. 实现 LogLog, SemiLog, Cartesian 三个图表的绘制逻辑。
 * 2. 针对压力恢复试井，实现 Horner Plot (P vs lg((tp+dt)/dt)) 绘制与线性拟合。
 * 3. 自动计算初始地层压力并在图表上显示。
 */

#include "fittingchart.h"
#include <cmath>
#include <algorithm>
#include <QDebug>

FittingChart::FittingChart(QObject *parent)
    : QObject(parent), m_plotLogLog(nullptr), m_plotSemiLog(nullptr), m_plotCartesian(nullptr), m_calculatedPi(0.0)
{
}

void FittingChart::initializeCharts(MouseZoom* logLog, MouseZoom* semiLog, MouseZoom* cartesian)
{
    m_plotLogLog = logLog;
    m_plotSemiLog = semiLog;
    m_plotCartesian = cartesian;
}

void FittingChart::setObservedData(const QVector<double>& t, const QVector<double>& deltaP,
                                   const QVector<double>& deriv, const QVector<double>& rawP)
{
    m_obsT = t;
    m_obsDeltaP = deltaP;
    m_obsDeriv = deriv;
    m_obsRawP = rawP;
}

void FittingChart::setSettings(const FittingDataSettings& settings)
{
    m_settings = settings;
}

void FittingChart::plotAll(const QVector<double>& t_model, const QVector<double>& p_model, const QVector<double>& d_model, bool isModelValid)
{
    if (!m_plotLogLog || !m_plotSemiLog || !m_plotCartesian) return;

    plotLogLog(t_model, p_model, d_model, isModelValid);
    plotSemiLog(t_model, p_model, d_model, isModelValid);
    plotCartesian(t_model, p_model, d_model, isModelValid);

    m_plotLogLog->replot();
    m_plotSemiLog->replot();
    m_plotCartesian->replot();
}

void FittingChart::plotLogLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel)
{
    MouseZoom* plot = m_plotLogLog;
    plot->clearGraphs();
    plot->clearItems(); // 清除旧的文本框

    // 1. 实测数据
    QVector<double> vt, vp, vd;
    for(int i=0; i<m_obsT.size(); ++i) {
        if (m_obsT[i] > 1e-10 && m_obsDeltaP[i] > 1e-10) {
            vt << m_obsT[i];
            vp << m_obsDeltaP[i];
            if (i < m_obsDeriv.size() && m_obsDeriv[i] > 1e-10) vd << m_obsDeriv[i];
            else vd << 1e-10; // log scale placeholder
        }
    }

    plot->addGraph(); // 0: 实测压差
    plot->graph(0)->setData(vt, vp);
    plot->graph(0)->setPen(Qt::NoPen);
    plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QColor(0, 100, 0), 6));
    plot->graph(0)->setName("实测压差");

    plot->addGraph(); // 1: 实测导数
    plot->graph(1)->setData(vt, vd);
    plot->graph(1)->setPen(Qt::NoPen);
    plot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssTriangle, Qt::magenta, 6));
    plot->graph(1)->setName("实测导数");

    // 2. 理论曲线
    if (hasModel) {
        QVector<double> vtm, vpm, vdm;
        for(int i=0; i<tm.size(); ++i) {
            if(tm[i] > 1e-10) {
                vtm << tm[i];
                vpm << (pm[i] > 1e-10 ? pm[i] : 1e-10);
                vdm << (dm[i] > 1e-10 ? dm[i] : 1e-10);
            }
        }
        plot->addGraph(); // 2: 理论压差
        plot->graph(2)->setData(vtm, vpm);
        plot->graph(2)->setPen(QPen(Qt::red, 2));
        plot->graph(2)->setName("理论压差");

        plot->addGraph(); // 3: 理论导数
        plot->graph(3)->setData(vtm, vdm);
        plot->graph(3)->setPen(QPen(Qt::blue, 2));
        plot->graph(3)->setName("理论导数");
    }

    // 3. 设置坐标轴范围
    plot->xAxis->setLabel("时间 Time (h)");
    plot->yAxis->setLabel("压差 & 导数 (MPa)");
    plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    plot->rescaleAxes();
    // 稍微扩展范围
    plot->xAxis->scaleRange(1.1, plot->xAxis->range().center());
    plot->yAxis->scaleRange(1.1, plot->yAxis->range().center());

    // 4. 显示计算结果 (如果是恢复试井)
    if (m_settings.testType == Test_Buildup && m_calculatedPi > 1e-6) {
        showResultOnLogPlot();
    }
}

void FittingChart::plotSemiLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel)
{
    MouseZoom* plot = m_plotSemiLog;
    plot->clearGraphs();
    Q_UNUSED(dm); // 半对数通常不画导数

    // 判断模式：压力降落(Drawdown) vs 压力恢复(Buildup)
    if (m_settings.testType == Test_Buildup && m_settings.producingTime > 0) {
        // === 压力恢复 Horner Plot 模式 ===

        QVector<double> hornerX, hornerY;
        double tp = m_settings.producingTime;

        for(int i=0; i<m_obsT.size(); ++i) {
            double dt = m_obsT[i];
            if (dt > 1e-6 && i < m_obsRawP.size()) {
                double val = (tp + dt) / dt;
                if (val > 0) {
                    hornerX << log10(val); // 绘制的是 log 值
                    hornerY << m_obsRawP[i];
                }
            }
        }

        plot->addGraph();
        plot->graph(0)->setData(hornerX, hornerY);
        plot->graph(0)->setPen(Qt::NoPen);
        plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QColor(0, 0, 180), 5));
        plot->graph(0)->setName("实测压力");

        // 计算并绘制拟合直线
        m_calculatedPi = calculateHornerPressure();
        if (m_calculatedPi > 0) {
            // [修复] 移除了未使用的 minX, maxX, lastX 变量
            if (!hornerX.isEmpty()) {
                double startX = hornerX.first(); // 早期时间，X大

                // 为了画线，取两端
                QVector<double> lineX, lineY;
                lineX << startX << 0.0;
                // 斜率 m = (Y - Pi) / X
                // Y = mX + Pi
                // 简单拟合用于绘图
                int nPoints = hornerX.size();
                if (nPoints > 2) {
                    double x1 = hornerX[nPoints-1]; double y1 = hornerY[nPoints-1]; // 最晚时间点
                    double x2 = hornerX[nPoints/2]; double y2 = hornerY[nPoints/2]; // 中间点
                    if (std::abs(x1-x2) > 1e-6) {
                        double m = (y1 - y2) / (x1 - x2);
                        // Y(0) = Pi
                        // Y(startX) = m * startX + Pi
                        lineY << (m * startX + m_calculatedPi) << m_calculatedPi;
                    } else {
                        lineY << hornerY.first() << m_calculatedPi;
                    }
                } else {
                    lineY << hornerY.first() << m_calculatedPi;
                }

                plot->addGraph();
                plot->graph(1)->setData(lineX, lineY);
                plot->graph(1)->setPen(QPen(Qt::red, 2, Qt::DashLine));
                plot->graph(1)->setName("Horner 拟合线");
            }
        }

        plot->xAxis->setLabel("Horner 时间比 lg((tp+dt)/dt)");
        plot->yAxis->setLabel("地层压力 P (MPa)");
        plot->xAxis->setScaleType(QCPAxis::stLinear); // 因为我们手动取了 log10
        plot->yAxis->setScaleType(QCPAxis::stLinear);
        plot->xAxis->setRangeReversed(true); // Horner 图习惯 X 轴从大到小
        plot->rescaleAxes();
        // 确保能看到 X=0
        double upperX = plot->xAxis->range().upper;
        plot->xAxis->setRange(upperX, 0.0);

    } else {
        // === 压力降落 (原逻辑) ===
        QVector<double> vt, vp;
        for(int i=0; i<m_obsT.size(); ++i) {
            if(m_obsT[i] > 1e-10) {
                vt << m_obsT[i];
                vp << m_obsDeltaP[i];
            }
        }

        plot->addGraph();
        plot->graph(0)->setData(vt, vp);
        plot->graph(0)->setPen(Qt::NoPen);
        plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QColor(0, 100, 0), 6));
        plot->graph(0)->setName("实测压差");

        if (hasModel) {
            QVector<double> vtm, vpm;
            for(int i=0; i<tm.size(); ++i) {
                if(tm[i] > 1e-10) {
                    vtm << tm[i];
                    vpm << (pm[i] > 1e-10 ? pm[i] : 1e-10);
                }
            }
            plot->addGraph();
            plot->graph(1)->setData(vtm, vpm);
            plot->graph(1)->setPen(QPen(Qt::red, 2));
            plot->graph(1)->setName("理论压差");
        }

        plot->xAxis->setLabel("时间 Time (h)");
        plot->yAxis->setLabel("压差 Delta P (MPa)");
        plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
        plot->yAxis->setScaleType(QCPAxis::stLinear);
        plot->xAxis->setRangeReversed(false);
        plot->rescaleAxes();
    }
}

void FittingChart::plotCartesian(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel)
{
    MouseZoom* plot = m_plotCartesian;
    plot->clearGraphs();
    Q_UNUSED(dm);

    QVector<double> vt, vp;
    // 笛卡尔图通常画压差
    for(int i=0; i<m_obsT.size(); ++i) {
        vt << m_obsT[i];
        vp << m_obsDeltaP[i];
    }

    plot->addGraph();
    plot->graph(0)->setData(vt, vp);
    plot->graph(0)->setPen(Qt::NoPen);
    plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QColor(0, 100, 0), 6));
    plot->graph(0)->setName("实测压差");

    if (hasModel) {
        QVector<double> vtm, vpm;
        for(int i=0; i<tm.size(); ++i) {
            vtm << tm[i];
            vpm << pm[i];
        }
        plot->addGraph();
        plot->graph(1)->setData(vtm, vpm);
        plot->graph(1)->setPen(QPen(Qt::red, 2));
        plot->graph(1)->setName("理论压差");
    }

    plot->xAxis->setLabel("时间 Time (h)");
    plot->yAxis->setLabel("压差 Delta P (MPa)");
    plot->xAxis->setScaleType(QCPAxis::stLinear);
    plot->yAxis->setScaleType(QCPAxis::stLinear);
    plot->rescaleAxes();
}

void FittingChart::plotSampledPoints(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d)
{
    if (!m_plotLogLog) return;

    QCPGraph* gP = m_plotLogLog->addGraph();
    gP->setData(t, p);
    gP->setPen(Qt::NoPen);
    gP->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(QColor(0, 100, 0)), QBrush(QColor(0, 100, 0)), 6));
    gP->setName("抽样压差");

    QCPGraph* gD = m_plotLogLog->addGraph();
    gD->setData(t, d);
    gD->setPen(Qt::NoPen);
    gD->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssTriangle, QPen(Qt::magenta), QBrush(Qt::magenta), 6));
    gD->setName("抽样导数");
}

double FittingChart::calculateHornerPressure()
{
    // 简单的线性回归，基于最后一部分数据 (径向流阶段)
    if (m_obsT.isEmpty() || m_obsRawP.isEmpty() || m_settings.producingTime <= 0) return 0.0;

    QVector<double> X, Y;
    for(int i=0; i<m_obsT.size(); ++i) {
        double dt = m_obsT[i];
        if(dt > 1e-5 && i < m_obsRawP.size()) {
            double val = (m_settings.producingTime + dt) / dt;
            if(val > 0) {
                X << log10(val);
                Y << m_obsRawP[i];
            }
        }
    }

    if (X.size() < 5) return 0.0;

    // 选取最后 30% 的点进行拟合
    int nPoints = X.size();
    int fitCount = nPoints * 0.3;
    if (fitCount < 3) fitCount = nPoints;

    double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
    int count = 0;

    for(int i = nPoints - fitCount; i < nPoints; ++i) {
        sumX += X[i];
        sumY += Y[i];
        sumXY += X[i] * Y[i];
        sumXX += X[i] * X[i];
        count++;
    }

    if (std::abs(count * sumXX - sumX * sumX) < 1e-9) return 0.0;

    double slope = (count * sumXY - sumX * sumY) / (count * sumXX - sumX * sumX);
    double intercept = (sumY - slope * sumX) / count;

    // Horner Extrapolation: log = 0, P = intercept
    return intercept;
}

double FittingChart::getCalculatedInitialPressure() const
{
    return m_calculatedPi;
}

void FittingChart::showResultOnLogPlot()
{
    if(!m_plotLogLog) return;

    QCPItemText *textLabel = new QCPItemText(m_plotLogLog);
    textLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignRight);
    textLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
    textLabel->position->setCoords(0.95, 0.05); // 右上角
    textLabel->setText(QString("Horner推算Pi: %1 MPa").arg(m_calculatedPi, 0, 'f', 2));
    textLabel->setFont(QFont("Microsoft YaHei", 10, QFont::Bold));
    textLabel->setColor(Qt::red);
    textLabel->setBrush(QBrush(QColor(255, 255, 255, 200)));
    textLabel->setPadding(QMargins(5, 5, 5, 5));
    textLabel->setPen(QPen(Qt::black));
}
