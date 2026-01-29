/*
 * 文件名: wt_fittingwidget.h
 * 文件作用: 试井拟合分析主界面头文件
 * 功能描述:
 * 1. 管理 FittingCore 实例，处理核心拟合任务。
 * 2. 负责界面元素的初始化、布局与交互响应。
 * 3. 负责数据的导入导出与报告生成。
 * 4. 集成 FittingChart 进行绘图管理。
 */

#ifndef WT_FITTINGWIDGET_H
#define WT_FITTINGWIDGET_H

#include <QWidget>
#include <QMap>
#include <QVector>
#include <QJsonObject>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <QResizeEvent>

#include "modelmanager.h"
#include "fittingparameterchart.h"
#include "chartwidget.h"
#include "mousezoom.h"
#include "fittingcore.h"
#include "fittingsamplingdialog.h"
#include "fittingreport.h"
#include "fittingchart.h"

namespace Ui {
class FittingWidget;
}

class FittingWidget : public QWidget
{
    Q_OBJECT

public:
    explicit FittingWidget(QWidget *parent = nullptr);
    ~FittingWidget();

    // 设置模型管理器
    void setModelManager(ModelManager* m);

    // 设置项目数据模型集合
    void setProjectDataModels(const QMap<QString, QStandardItemModel*>& models);

    // [新增] 兼容性接口：仅设置压差数据 (rawP 默认为空)
    // 解决 fittingpage.cpp 中的编译错误
    void setObservedData(const QVector<double>& t, const QVector<double>& deltaP, const QVector<double>& d);

    // 完整接口：设置所有观测数据 (包含原始压力)
    void setObservedData(const QVector<double>& t, const QVector<double>& deltaP,
                         const QVector<double>& d, const QVector<double>& rawP);

    void updateBasicParameters();
    QJsonObject getJsonState() const;
    void loadFittingState(const QJsonObject& root);
    QString getPlotImageBase64(MouseZoom* plot);

protected:
    void resizeEvent(QResizeEvent* event) override;

signals:
    void sigRequestSave();

private slots:
    void on_btnLoadData_clicked();
    void on_btnSelectParams_clicked();
    void on_btnResetParams_clicked();
    void on_btn_modelSelect_clicked();
    void onSliderWeightChanged(int value);
    void on_btnRunFit_clicked();
    void on_btnStop_clicked();
    void onFitFinished();
    void onOpenSamplingSettings();
    void on_btnExportData_clicked();
    void on_btnExportReport_clicked();
    void onExportCurveData();
    void on_btnImportModel_clicked();
    void on_btnSaveFit_clicked();
    void updateModelCurve(const QMap<QString, double>* explicitParams = nullptr);
    void onIterationUpdate(double err, const QMap<QString,double>& p, const QVector<double>& t, const QVector<double>& p_curve, const QVector<double>& d_curve);

private:
    Ui::FittingWidget *ui;
    ModelManager* m_modelManager;
    FittingCore* m_core;
    FittingChart* m_chartManager;

    QMdiArea* m_mdiArea;
    ChartWidget* m_chartLogLog;
    ChartWidget* m_chartSemiLog;
    ChartWidget* m_chartCartesian;
    QMdiSubWindow* m_subWinLogLog;
    QMdiSubWindow* m_subWinSemiLog;
    QMdiSubWindow* m_subWinCartesian;
    MouseZoom* m_plotLogLog;
    MouseZoom* m_plotSemiLog;
    MouseZoom* m_plotCartesian;

    FittingParameterChart* m_paramChart;
    QMap<QString, QStandardItemModel*> m_dataMap;
    ModelManager::ModelType m_currentModelType;

    // 观测数据
    QVector<double> m_obsTime;
    QVector<double> m_obsDeltaP;
    QVector<double> m_obsDerivative;
    QVector<double> m_obsRawP;

    bool m_isFitting;
    bool m_isCustomSamplingEnabled;
    QList<SamplingInterval> m_customIntervals;

    // 内部初始化
    void setupPlot();
    void initializeDefaultModel();
    void layoutCharts();

    QVector<double> parseSensitivityValues(const QString& text);
};

#endif // WT_FITTINGWIDGET_H
