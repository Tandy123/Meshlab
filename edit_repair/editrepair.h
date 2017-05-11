/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef EDIT_REPAIR_PLUGIN_H
#define EDIT_REPAIR_PLUGIN_H

#include <QObject>
#include <common/interfaces.h>
#include "repairDialog.h"
#include "Teethtree.h"
#include "circle.h"

// #include <vcg/Eigen/Dense>//自带的Eigen有bug
// #include <vcg/Eigen/Sparse>
// #include <vcg/Eigen/Core>
#include <Dense>
#include <Sparse>
#include <Core>

class RepairDialog;

class ProLine{//自定义数据结果，用于存储投影线
public:
	vcg::Point3f lp;
	vcg::Point3f ldir;
};


class EditRepairPlugin : public QObject, public MeshEditInterface
{
	Q_OBJECT
	Q_INTERFACES(MeshEditInterface)
		
public:
	enum TranMode{
		NONE,
		IMPORTCOORD,
		APPLY,
		ERRORMAP,
		CUT,
		FINDNESTP,
		EXTENDB,
		LINKP
	};
    EditRepairPlugin();
    virtual ~EditRepairPlugin() {
		delete stagesDialog;
	}

    static const QString Info();

    bool StartEdit(MeshDocument &_md, GLArea *gla );
    void EndEdit(MeshModel &/*m*/, GLArea * /*parent*/){};
    void Decorate(MeshModel &/*m*/, GLArea * /*parent*/, QPainter *p);
    void mousePressEvent(QMouseEvent *, MeshModel &, GLArea * ) {};
    void mouseMoveEvent(QMouseEvent *, MeshModel &, GLArea * ) {};
    void mouseReleaseEvent(QMouseEvent *event, MeshModel &/*m*/, GLArea * );

	virtual void LayerChanged(MeshDocument &md, MeshModel &oldMeshModel, GLArea *parent){};

	

    QPoint cur;
	QFont qFont;
    bool haveToPick;
	bool importset;

	GLArea *gla;

	
private:
	TranMode current_mode;
	RepairDialog *stagesDialog;
	MeshModel *mesh;
	MeshModel *mesh2;
	MeshModel *mesh_tem;//template mesh
	MeshModel *mesh_tem_d;//template mesh deformation
	MeshModel *repaired_mesh;
	MeshModel *mesh_copy;
	MeshDocument *md;

	float DIF;
	QString currentfilepath;
	std::vector <IdMesh> MeshSet;
	vcg::Point3f barycentre;
	vcg::Point3f xcoord;
	vcg::Point3f ycoord;
	vcg::Point3f zcoord;
	vcg::Point3f barycentreTem;
	vcg::Point3f xcoordTem;
	vcg::Point3f ycoordTem;
	vcg::Point3f zcoordTem;
	double xTrans;
	double yTrans;
	double zTrans;
	double scaleTrans;
	std::vector <int> borders;
	std::vector <int> borders_root;
	std::vector <int> teethBordersInOrder;
	std::vector <int> teethBordersInOrder2;
	std::vector <int> rootBordersInOrder;

	std::vector <vcg::Point3f> vc1;
	std::vector <vcg::Point3f> vc2;
	std::vector <vcg::Point3f> vc3;

	std::vector <int> crownBoundaryT;
	std::vector <int> crownInnerT;
	std::vector <int> rootT;
	std::vector <int> rootBoundaryT;

	int teethBordersInOrderStart;
	int teethBordersInOrderEnd;
	int rootBordersInOrderEnd;

	bool getSE;

	std::vector <Circle*> circles;
	vcg::Point3f angleV1;
	vcg::Point3f angleV2;
	int flag1,flag2;
	int flagH1,flagH2;

	TranMode CurrentMode() {return current_mode;}
	void SetMode(TranMode m) {current_mode = m;}
	bool intersectTri(ProLine pl,vcg::Point3f a, vcg::Point3f b, vcg::Point3f c, vcg::Point3f fn, float tmax, float tmin,float &t);
	bool intersectTri2(ProLine pl,vcg::Point3f a, vcg::Point3f b, vcg::Point3f c, vcg::Point3f fn, float tmax, float tmin,float &t);
	void sortTeeth();

	void DrawLocalAxis();
	void DrawPyramid(vcg::Point3f vec,int loc,vcg::Point3f localaxis[3]);

	bool IntersectionPlaneMesh(CMeshO & m,vcg::Plane3f  pl,CMeshO & em);
	void FindNeastestP();
	void DrawColorBar(QPainter *p);
	//float SignedDistancePlanePoint(const Plane3f & plane, const Point3f & point);

public slots:
	void slotProject();
	void slotDeformationPro();
	void slotDeformationPro2();
	void slotDeformationPro3();
	void slotDeformationNonPro();
	void slotDeformationAuto();
	void slotDeformationNonPro2();
	void slotDeformationNonPro3();
	void slotDeformationNonPro4();
	void slotDeformationAuto2();
	void slotDeletePro();
	void slotImportSet();
	void slotFillGap();
	void slotFillGap1();
	void slotFillGap2();
	void slotFindBoundary();//Qian
	void slotFindBoundary1();
	void slotSortBorder();
	void slotSortBorder2();
	void slotCopyMesh();
	void slotImportSetTem();
	void slotAlignTem();
	void slotAdjustX(double);
	void slotAdjustY(double);
	void slotAdjustZ(double);
	void slotAdjustScale(double);
	void slotSaveAdjust();
	void slotProjectBarycentre();
	void slotProjectBarycentre2();
	void slotProjectBarycentre3();
	void slotProjectBarycentre4();
	void slotExportMesh();
	void slotDeformationLaplacian();
	void slotGetBoundaryTem();
	void slotSmoothRoot();
	void slotSmoothGap();
	void slotSmoothGap2();

	void slotDefomationTem();
	void slotRepairByTem();


	void slotCut();
	void slotCut1();

	void slotMerge();
	void slotNormalize();
	void slotExportPoints();
	void slotExportPoints2();
	void slotErrorMap();
	void slotErrorMap1();
	void slotErrorMap2();
	void slotInvertNormal();
	void slotInvertNormalFace();

	void slotFindNeastestPoints();
	void slotExtendBorder();
	void slotLinkPoints();
	void slotImportInsectionP();
};

#endif
