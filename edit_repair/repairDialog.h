#ifndef REPAIRDIALOG_H
#define REPAIRDIALOG_H


#include "ui_repairDialog.h"
#include "editrepair.h"
#include <QtGui>
#include <QWidget>
#include <QtGui/QDockWidget>

class EditRepairPlugin;

class RepairDialog : public QDockWidget{
	Q_OBJECT
public:
	//GetPointsDialog(QWidget *parent);
	RepairDialog(QWidget *parent,EditRepairPlugin *_edit);
	~RepairDialog(){};

	//Ui::PlanTable ui;
	Ui::Form ui;
	virtual void closeEvent ( QCloseEvent * event )	;

// protected:
// 	void paintEvent(QPaintEvent *event);

signals:
	void closing();
	void updateMeshSetVisibilities();

};

#endif

