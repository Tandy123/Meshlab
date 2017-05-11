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
/****************************************************************************
History
$Log: meshedit.cpp,v $
****************************************************************************/
#include <QtGui>

#include <math.h>
#include <stdlib.h>
#include <meshlab/glarea.h>
#include "editrepair.h"
#include <wrap/gl/pick.h>
#include<vcg/complex/trimesh/append.h>
#include<wrap/qt/gl_label.h>
#include <vcg/simplex/face/jumping_pos.h>
#include <vcg/complex/trimesh/clean.h>

#include <vcg/space/index/aabb_binary_tree/aabb_binary_tree.h>
#include <vcg/complex/trimesh/smooth.h>

#include <fstream> 
#include<sstream> 

#include <time.h>

using namespace std;
using namespace vcg;
using namespace Eigen;
vector<vector<int>> faceids;

typedef vcg::AABBBinaryTreeIndex<CFaceO, float, vcg::EmptyClass> AIndex;

inline float ComMax(float a, float b){
	float c;
	if (a < b)
	{
		c = b;
	}else{
		c = a;
	}
	return c;
}

inline float ComMin(float a, float b){
	float c;
	if (a > b)
	{
		c = b;
	}else{
		c = a;
	}
	return c;
}

MatrixXf RotationVector(Point3f &srcv,Point3f &dstv){
	vcg::Point3f axis = srcv^dstv; //旋转轴
	float lenaxis = sqrt(axis.X()*axis.X()+axis.Y()*axis.Y()+axis.Z()*axis.Z());
	axis /= lenaxis;

	float lensrcv =sqrt(srcv.X()*srcv.X() + srcv.Y()*srcv.Y() + srcv.Z()*srcv.Z());
	srcv /= lensrcv;
	float lendstv =sqrt(dstv.X()*dstv.X() + dstv.Y()*dstv.Y() + dstv.Z()*dstv.Z());
	dstv /= lendstv;
	float angle = acos(srcv*dstv);

	Eigen::MatrixXf matrix(3,3);
	matrix(0,0) = (1.0-cos(angle))*axis.X()*axis.X()+cos(angle);
	matrix(0,1) = (1.0-cos(angle))*axis.X()*axis.Y()+sin(angle)*axis.Z();
	matrix(0,2) = (1.0-cos(angle))*axis.X()*axis.Z()-sin(angle)*axis.Y();
	matrix(1,0) = (1.0-cos(angle))*axis.X()*axis.Y()-sin(angle)*axis.Z();
	matrix(1,1) = (1.0-cos(angle))*axis.Y()*axis.Y()+cos(angle);
	matrix(1,2) = (1.0-cos(angle))*axis.Y()*axis.Z()+sin(angle)*axis.X();
	matrix(2,0) = (1.0-cos(angle))*axis.X()*axis.Z()+sin(angle)*axis.Y();
	matrix(2,1) = (1.0-cos(angle))*axis.Y()*axis.Z()-sin(angle)*axis.X();
	matrix(2,2) = (1.0-cos(angle))*axis.Z()*axis.Z()+cos(angle);

	Vector3f vtemp(srcv.X(),srcv.Y(),srcv.Z());
	MatrixXf mf(3,1);
	mf =matrix.transpose()*vtemp;//注意不要漏掉转置
	Point3f dstv1 = Point3f(mf(0,0),mf(1,0),mf(2,0));	

	float dis = Distance(dstv,dstv1);//前提是dstv必须是单位向量

	if (dis>0.01)
	{
		angle = -1.0 * angle;
		matrix(0,0) = (1.0-cos(angle))*axis.X()*axis.X()+cos(angle);
		matrix(0,1) = (1.0-cos(angle))*axis.X()*axis.Y()+sin(angle)*axis.Z();
		matrix(0,2) = (1.0-cos(angle))*axis.X()*axis.Z()-sin(angle)*axis.Y();
		matrix(1,0) = (1.0-cos(angle))*axis.X()*axis.Y()-sin(angle)*axis.Z();
		matrix(1,1) = (1.0-cos(angle))*axis.Y()*axis.Y()+cos(angle);
		matrix(1,2) = (1.0-cos(angle))*axis.Y()*axis.Z()+sin(angle)*axis.X();
		matrix(2,0) = (1.0-cos(angle))*axis.X()*axis.Z()+sin(angle)*axis.Y();
		matrix(2,1) = (1.0-cos(angle))*axis.Y()*axis.Z()-sin(angle)*axis.X();
		matrix(2,2) = (1.0-cos(angle))*axis.Z()*axis.Z()+cos(angle);
	}
	return matrix.transpose();
}

EditRepairPlugin::EditRepairPlugin() {
	qFont.setFamily("Helvetica");
	qFont.setPixelSize(12);
	stagesDialog = 0;
}

const QString EditRepairPlugin::Info() 
{
	return tr("Cut the Teeth.");
}

void EditRepairPlugin::mouseReleaseEvent(QMouseEvent * event, MeshModel &/*m*/, GLArea * gla)
{
	gla->update();
	cur=event->pos();
	haveToPick = true;
}

void EditRepairPlugin::DrawPyramid(vcg::Point3f vec,int loc,vcg::Point3f localaxis[3]){
	unsigned int a,b;
	for(unsigned int i=0;i<3;i++){
		if(loc != i){
			a = i;
			break;
		}
	}
	for(unsigned int i=0;i<3;i++){
		if(loc != i && a!= i){
			b = i;
			break;
		}
	}
	std::vector<vcg::Point3f> m_pyramid,m_quad;
	vcg::Point3f vertex;
	m_pyramid.resize(12);
	m_quad.resize(4);  
	vcg::Point3f offsetv;
	offsetv.X() = localaxis[loc].X() *0.1;offsetv.Y() = localaxis[loc].Y() *0.1;offsetv.Z() = localaxis[loc].Z() *0.1;
	vec += offsetv;

	m_pyramid.push_back(vec);
	vertex.X() = vec.X() + 0.3*localaxis[a].X() - 0.5*localaxis[loc].X();
	vertex.Y() = vec.Y() + 0.3*localaxis[a].Y() - 0.5*localaxis[loc].Y();
	vertex.Z() = vec.Z() + 0.3*localaxis[a].Z() - 0.5*localaxis[loc].Z();
	m_quad.push_back(vertex);
	m_pyramid.push_back(vertex); 
	vertex.X() = vec.X() + 0.3*localaxis[b].X() - 0.5*localaxis[loc].X();
	vertex.Y() = vec.Y() + 0.3*localaxis[b].Y() - 0.5*localaxis[loc].Y();
	vertex.Z() = vec.Z() + 0.3*localaxis[b].Z() - 0.5*localaxis[loc].Z();
	m_pyramid.push_back(vertex); 
	m_quad.push_back(vertex);

	m_pyramid.push_back(vec);
	vertex.X() = vec.X() + 0.3*localaxis[b].X() - 0.5*localaxis[loc].X();
	vertex.Y() = vec.Y() + 0.3*localaxis[b].Y() - 0.5*localaxis[loc].Y();
	vertex.Z() = vec.Z() + 0.3*localaxis[b].Z() - 0.5*localaxis[loc].Z();
	m_pyramid.push_back(vertex); 
	vertex.X() = vec.X() - 0.3*localaxis[a].X() - 0.5*localaxis[loc].X();
	vertex.Y() = vec.Y() - 0.3*localaxis[a].Y() - 0.5*localaxis[loc].Y();
	vertex.Z() = vec.Z() - 0.3*localaxis[a].Z() - 0.5*localaxis[loc].Z();
	m_pyramid.push_back(vertex); 
	m_quad.push_back(vertex);

	m_pyramid.push_back(vec);
	vertex.X() = vec.X() - 0.3*localaxis[a].X() - 0.5*localaxis[loc].X();
	vertex.Y() = vec.Y() - 0.3*localaxis[a].Y() - 0.5*localaxis[loc].Y();
	vertex.Z() = vec.Z() - 0.3*localaxis[a].Z() - 0.5*localaxis[loc].Z();
	m_pyramid.push_back(vertex); 
	vertex.X() = vec.X() - 0.3*localaxis[b].X() - 0.5*localaxis[loc].X();
	vertex.Y() = vec.Y() - 0.3*localaxis[b].Y() - 0.5*localaxis[loc].Y();
	vertex.Z() = vec.Z() - 0.3*localaxis[b].Z() - 0.5*localaxis[loc].Z();
	m_pyramid.push_back(vertex);
	m_quad.push_back(vertex);

	m_pyramid.push_back(vec);
	vertex.X() = vec.X() - 0.3*localaxis[b].X() - 0.5*localaxis[loc].X();
	vertex.Y() = vec.Y() - 0.3*localaxis[b].Y() - 0.5*localaxis[loc].Y();
	vertex.Z() = vec.Z() - 0.3*localaxis[b].Z() - 0.5*localaxis[loc].Z();
	m_pyramid.push_back(vertex); 
	vertex.X() = vec.X() + 0.3*localaxis[a].X() - 0.5*localaxis[loc].X();
	vertex.Y() = vec.Y() + 0.3*localaxis[a].Y() - 0.5*localaxis[loc].Y();
	vertex.Z() = vec.Z() + 0.3*localaxis[a].Z() - 0.5*localaxis[loc].Z();
	m_pyramid.push_back(vertex);
	glColor3f(0.0f, 0.0f, 0.0f);
	glBegin(GL_TRIANGLES);
	for(unsigned int i=0 ;i<m_pyramid.size();i++){
		glVertex3f(m_pyramid[i].X(),m_pyramid[i].Y(),m_pyramid[i].Z());  
	}
	glEnd();
	glBegin(GL_QUADS);
	for(unsigned int i=0 ;i<m_quad.size();i++){
		glVertex3f(m_quad[i].X(),m_quad[i].Y(),m_quad[i].Z());  
	}
	glEnd();
}

void EditRepairPlugin::DrawLocalAxis(){
	vcg::Point3f _CMPCA[3];
	_CMPCA[0] = zcoord;
	_CMPCA[1] = xcoord;
	_CMPCA[2] = ycoord;

	glLineWidth(3.0);
	glColor3f(1.0f, 0.0f, 0.0f);
	//z轴
	glBegin(GL_LINES);
	glVertex3f(barycentre.X() - 10*zcoord.X(),barycentre.Y() - 10*zcoord.Y(),barycentre.Z() - 10*zcoord.Z());
	glVertex3f(barycentre.X() + 10*zcoord.X(),barycentre.Y() + 10*zcoord.Y(),barycentre.Z() + 10*zcoord.Z());
	glEnd();
	DrawPyramid(vcg::Point3f(barycentre.X() + 10*zcoord.X(),barycentre.Y() + 10*zcoord.Y(),barycentre.Z() + 10*zcoord.Z()),0,_CMPCA);

	glColor3f(0.0f, 1.0f, 0.0f);
	//x轴
	glBegin(GL_LINES);
	glVertex3f(barycentre.X() - 10*xcoord.X(),barycentre.Y() - 10*xcoord.Y(),barycentre.Z() - 10*xcoord.Z());
	glVertex3f(barycentre.X() + 10*xcoord.X(),barycentre.Y() + 10*xcoord.Y(),barycentre.Z() + 10*xcoord.Z());
	glEnd();
	DrawPyramid(vcg::Point3f(barycentre.X() + 10*xcoord.X(),barycentre.Y() + 10*xcoord.Y(),barycentre.Z() + 10*xcoord.Z()),1,_CMPCA);


	glColor3f(0.0f, 0.0f, 1.0f);
	//y轴
	glBegin(GL_LINES);
	glVertex3f(barycentre.X() - 10*ycoord.X(),barycentre.Y() - 10*ycoord.Y(),barycentre.Z() - 10*ycoord.Z());
	glVertex3f(barycentre.X() + 10*ycoord.X(),barycentre.Y() + 10*ycoord.Y(),barycentre.Z() + 10*ycoord.Z());
	glEnd();
	DrawPyramid(vcg::Point3f(barycentre.X() + 10*ycoord.X(),barycentre.Y() + 10*ycoord.Y(),barycentre.Z() + 10*ycoord.Z()),2,_CMPCA);
}

void EditRepairPlugin::Decorate(MeshModel &m, GLArea * gla, QPainter *p)
{
	if(CurrentMode() == TranMode::FINDNESTP){
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Red);
		glPointSize(7.0f);
		glBegin(GL_POINTS);
		glVertex(mesh->cm.vert[flag1].P());
		glVertex(mesh->cm.vert[flag2].P());
		glVertex(mesh2->cm.vert[flagH1].P());
		glVertex(mesh2->cm.vert[flagH2].P());
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}
	if ( importset )
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		/*		if(isrepairing_axis)*/
		DrawLocalAxis();
		// 		else
		// 			DrawRepairAxis(rotationaxis);
		glPopAttrib();
		glPopMatrix();
	}
	if(borders.size()>0 && teethBordersInOrder.size() == 0)
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Red);
		glPointSize(7.0f);
		glBegin(GL_POINTS);
		for(int i =0;i<borders.size();i++)
			glVertex(mesh->cm.vert[borders[i]].P());
		glEnd();
		glPopMatrix();
	}
	if(borders_root.size()>0 && rootBordersInOrder.size() == 0)
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Red);
		glPointSize(7.0f);
		glBegin(GL_POINTS);
		for(int i =0;i<borders_root.size();i++)
			glVertex(mesh_tem_d->cm.vert[borders_root[i]].P());
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}

	if(teethBordersInOrder.size()>0 && !getSE)
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Green);
		glPointSize(5.0f);
		glBegin(GL_POINTS);
		for(int i =0;i<teethBordersInOrder.size();i++)
			glVertex(mesh->cm.vert[teethBordersInOrder[i]].P());
		glEnd();
		glPopAttrib();
		glPopMatrix();

		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Blue);
		glPointSize(10.0f);
		glBegin(GL_POINTS);
		glVertex(mesh->cm.vert[teethBordersInOrder[0]].P());
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}
	if(teethBordersInOrder2.size()>0 && !getSE)
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Green);
		glPointSize(5.0f);
		glBegin(GL_POINTS);
		for(int i =0;i<teethBordersInOrder2.size();i++)
			glVertex(mesh2->cm.vert[teethBordersInOrder2[i]].P());
		glEnd();
		glPopAttrib();
		glPopMatrix();

		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Blue);
		glPointSize(10.0f);
		glBegin(GL_POINTS);
		glVertex(mesh2->cm.vert[teethBordersInOrder2[0]].P());
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}

	if(rootBordersInOrder.size()>0 && !getSE)
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Green);
		glPointSize(5.0f);
		glBegin(GL_POINTS);
		for(int i = 0;i<rootBordersInOrder.size();i++)
			glVertex(mesh_tem_d->cm.vert[rootBordersInOrder[i]].P());
		glEnd();
		glPopAttrib();
		glPopMatrix();

		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Blue);
		glPointSize(10.0f);
		glBegin(GL_POINTS);
		glVertex(mesh_tem_d->cm.vert[rootBordersInOrder[0]].P());
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}
	if (getSE){
		if(teethBordersInOrder.size()>0)
		{
			glPushMatrix();
			glMultMatrix(m.cm.Tr);
			glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
			glColor(Color4b::Green);
			glPointSize(5.0f);
			glBegin(GL_POINTS);
			for(int i =teethBordersInOrderStart;i<teethBordersInOrderEnd+1;i++)
				glVertex(mesh->cm.vert[teethBordersInOrder[i]].P());
			glEnd();
			glPopAttrib();
			glPopMatrix();
		}

		if(rootBordersInOrder.size()>0)
		{
			glPushMatrix();
			glMultMatrix(m.cm.Tr);
			glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
			glColor(Color4b::Green);
			glPointSize(5.0f);
			glBegin(GL_POINTS);
			for(int i = 0;i<rootBordersInOrderEnd+1;i++)
				glVertex(mesh_tem_d->cm.vert[rootBordersInOrder[i]].P());
			glEnd();
			glPopAttrib();
			glPopMatrix();
		}
	}
	if(vc1.size()>0)
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Yellow);
		glPointSize(7.0f);
		glBegin(GL_POINTS);
		for(int i =0;i<vc1.size();i++)
			glVertex(vc1[i]);
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}
	if(vc2.size()>0)
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Green);
		glPointSize(7.0f);
		glBegin(GL_POINTS);
		for(int i =0;i<vc2.size();i++)
			glVertex(vc2[i]);
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}
	if(vc3.size()>0)
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Blue);
		glPointSize(7.0f);
		glBegin(GL_POINTS);
		for(int i =0;i<vc3.size();i++)
			glVertex(vc3[i]);
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}
	if(crownBoundaryT.size()>0 && stagesDialog->ui.checkBoxShowPoints->checkState()== Qt::Checked)
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Blue);
		glPointSize(7.0f);
		glBegin(GL_POINTS);
		for(int i =0;i<crownBoundaryT.size();i++)
			glVertex(mesh_tem->cm.vert[crownBoundaryT[i]].P());
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}
	if(crownInnerT.size()>0 && stagesDialog->ui.checkBoxShowPoints->checkState()== Qt::Checked)
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Yellow);
		glPointSize(7.0f);
		glBegin(GL_POINTS);
		for(int i =0;i<crownInnerT.size();i++)
			glVertex(mesh_tem->cm.vert[crownInnerT[i]].P());
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}
	if(rootBoundaryT.size()>0 && stagesDialog->ui.checkBoxShowPoints->checkState()== Qt::Checked)
	{
		glPushMatrix();
		glMultMatrix(m.cm.Tr);
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT | GL_DEPTH_BUFFER_BIT);
		glColor(Color4b::Green);
		glPointSize(7.0f);
		glBegin(GL_POINTS);
		for(int i =0;i<rootBoundaryT.size();i++)
			glVertex(mesh_tem->cm.vert[rootBoundaryT[i]].P());
		glEnd();
		glPopAttrib();
		glPopMatrix();
	}
	if ( CurrentMode() == TranMode::ERRORMAP)
	{
		DrawColorBar(p);
	}
}

void EditRepairPlugin::slotProject()
{
	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	float t;
	for(vit = mesh_tem->cm.vert.begin(); vit != mesh_tem->cm.vert.end(); ++vit){
		for(fi = mesh->cm.face.begin(); fi != mesh->cm.face.end(); ++fi){
			ProLine pl;
			pl.lp = vit->P();
			pl.ldir = vit->N();
			if (intersectTri(pl,fi->V(0)->P(),fi->V(1)->P(),fi->V(2)->P(), fi->N(), 1.0, -1.0, t))
			{
				vit->C() = Color4b::Red;
				break;
			}
		}
	}
	mesh_tem->updateDataMask(MeshModel::MM_VERTCOLOR);
}

bool EditRepairPlugin::intersectTri(ProLine pl,vcg::Point3f pa, vcg::Point3f pb, vcg::Point3f pc, vcg::Point3f fn, float tmax, float tmin, float &t){
	float a = pa.X()-pb.X();
	float b = pa.Y()-pb.Y();
	float c = pa.Z()-pb.Z();
	float d = pa.X()-pc.X();
	float e = pa.Y()-pc.Y();
	float f = pa.Z()-pc.Z();
	float g = pl.ldir.X();
	float h = pl.ldir.Y();
	float i = pl.ldir.Z();
	float j = pa.X()-pl.lp.X();
	float k = pa.Y()-pl.lp.Y();
	float l = pa.Z()-pl.lp.Z();

	float M = a*(e*i-h*f)+b*(g*f-d*i)+c*(d*h-e*g);

	t = -(f*(a*k-j*b)+e*(j*c-a*l)+d*(b*l-k*c))/M;
	if(t<tmin || t>tmax)
		return false;
	float gama = (i*(a*k-j*b)+h*(j*c-a*l)+g*(b*l-k*c))/M;
	if (gama<0 || gama>1)
	{
		return false;
	}
	float beta = (j*(e*i-h*f)+k*(g*f-d*i)+l*(d*h-e*g))/M;
	if (beta<0 || beta>1-gama)
	{
		return false;
	}
	if (fn*pl.ldir<0)
	{
		return false;
	}
	//vcg::Point3f faceNormal = (pa+pb+pc)
	return true;
}

bool EditRepairPlugin::intersectTri2(ProLine pl,vcg::Point3f pa, vcg::Point3f pb, vcg::Point3f pc, vcg::Point3f fn, float tmax, float tmin, float &t){
	float a = pa.X()-pb.X();
	float b = pa.Y()-pb.Y();
	float c = pa.Z()-pb.Z();
	float d = pa.X()-pc.X();
	float e = pa.Y()-pc.Y();
	float f = pa.Z()-pc.Z();
	float g = pl.ldir.X();
	float h = pl.ldir.Y();
	float i = pl.ldir.Z();
	float j = pa.X()-pl.lp.X();
	float k = pa.Y()-pl.lp.Y();
	float l = pa.Z()-pl.lp.Z();

	float M = a*(e*i-h*f)+b*(g*f-d*i)+c*(d*h-e*g);

	t = -(f*(a*k-j*b)+e*(j*c-a*l)+d*(b*l-k*c))/M;
	if(t<tmin || t>tmax)
		return false;
	float gama = (i*(a*k-j*b)+h*(j*c-a*l)+g*(b*l-k*c))/M;
	if (gama<0 || gama>1)
	{
		return false;
	}
	float beta = (j*(e*i-h*f)+k*(g*f-d*i)+l*(d*h-e*g))/M;
	if (beta<0 || beta>1-gama)
	{
		return false;
	}
	if (fn*pl.ldir<0)
	{
		return false;
	}
	//vcg::Point3f faceNormal = (pa+pb+pc)
	return true;
}


void EditRepairPlugin::slotDeformationPro()
{
	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	float t;
	for(vit = mesh_tem_d->cm.vert.begin(); vit != mesh_tem_d->cm.vert.end(); ++vit){
		float tmin = 1000;
		ProLine pl;
		pl.lp = vit->P();
		pl.ldir = vit->N();
		for(fi = mesh->cm.face.begin(); fi != mesh->cm.face.end(); ++fi){
			if (intersectTri(pl,fi->V(0)->P(),fi->V(1)->P(),fi->V(2)->P(), fi->N(), 1.0, -1.0, t))
			{
				if (fabs(t)<fabs(tmin))
				{
					tmin = t;
				}
			}
		}
		if (tmin<1000)
		{
			vit->P() = pl.lp + pl.ldir*tmin;
			vit->C() = Color4b::Red;
		}
	}
	mesh_tem_d->updateDataMask(MeshModel::MM_VERTCOLOR);
	// 	vcg::tri::UpdateBounding<CMeshO>::Box(mesh_tem->cm);
	// 	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_tem->cm);
}

void EditRepairPlugin::slotDeformationPro2()
{
	clock_t start,finish;
	double totaltime;
	start=clock();

	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	float t;

	AIndex gIndex;
	gIndex.Set(mesh->cm.face.begin(), mesh->cm.face.end());

	for(vit = mesh_tem_d->cm.vert.begin(); vit != mesh_tem_d->cm.vert.end(); ++vit){
		if (!vit->IsD())
		{
			ProLine pl;
			pl.lp = vit->P();
			pl.ldir = (vit->P()-barycentre)/Distance(vit->P(),barycentre);

			const bool TEST_BACK_FACES = true;
			vcg::RayTriangleIntersectionFunctor<TEST_BACK_FACES> rayIntersector;
			const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();
			vcg::Ray3<AIndex::ScalarType, false> ray(barycentre, pl.ldir);

			AIndex::ObjPtr isectFace;
			AIndex::ScalarType rayT;
			vcg::EmptyClass a;
			isectFace = gIndex.DoRay(rayIntersector, a , ray, maxDist, rayT);

			if (isectFace != 0) {
				vit->P() = barycentre + pl.ldir*rayT;
				vit->C() = Color4b::Red;
				continue;
			}
		}
	}
	mesh_tem_d->updateDataMask(MeshModel::MM_VERTCOLOR);

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n此程序的运行时间为"<<totaltime<<endl;

	gla->update();
	// 	vcg::tri::UpdateBounding<CMeshO>::Box(mesh_tem->cm);
	// 	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_tem->cm);
}

void EditRepairPlugin::slotDeformationPro3()
{
	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	float t;
	for(vit = mesh_tem_d->cm.vert.begin(); vit != mesh_tem_d->cm.vert.end(); ++vit){
		float tmin = 1000;
		ProLine pl;
		pl.lp = vit->P();
		pl.ldir = (vit->P()-barycentre)/Distance(vit->P(),barycentre);
		for(fi = mesh->cm.face.begin(); fi != mesh->cm.face.end(); ++fi){
			vcg::Point3f fdir = (fi->V(0)->P()-barycentre)/Distance(fi->V(0)->P(),barycentre);
			if (intersectTri2(pl,fi->V(0)->P(),fi->V(1)->P(),fi->V(2)->P(), fdir, 10.0, -10.0, t))
			{
				if (fabs(t)<fabs(tmin))
				{
					tmin = t;
				}
			}
		}
		if (tmin<1000)
		{
			vit->P() = pl.lp + pl.ldir*tmin;
			vit->C() = Color4b::Red;
		}
	}
	mesh_tem_d->updateDataMask(MeshModel::MM_VERTCOLOR);
	// 	vcg::tri::UpdateBounding<CMeshO>::Box(mesh_tem->cm);
	// 	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_tem->cm);
}

void EditRepairPlugin::slotDeformationNonPro()
{
	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vitd;
	mesh_tem_d->updateDataMask(MeshModel::MM_VERTFACETOPO);
	mesh_tem_d->updateDataMask(MeshModel::MM_FACEFACETOPO);
	vit = mesh_tem->cm.vert.begin();
	for(vitd = mesh_tem_d->cm.vert.begin(); vitd != mesh_tem_d->cm.vert.end(); ++vitd){
		vcg::Color4f c1 = Color4f::Construct(vitd->C());
		if ( c1 != Color4<float>(Color4<float>::Red))//if white color?
		{
			vcg::Point3f tran = vcg::Point3f(0.0,0.0,0.0);
			int n = 0;
			if(!(*vitd).IsD()){//判断该点是否存在
				CVertexO *vii = &(*vitd);
				vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
				CVertexO* firstV = pos.VFlip();//取和vitd相邻的第一个点
				CVertexO* tempV= NULL;
				do
				{
					pos.NextE();
					tempV = pos.VFlip();
					int num = tempV - &*(mesh_tem_d->cm.vert.begin());
					tran += (tempV->P() - (vit + num)->P());
					n++;
				}while (firstV != tempV);//从第一个点开始按顺序遍历，知道回到第一个点为止
				int num2 = &*vitd - &*(mesh_tem_d->cm.vert.begin());
				vcg::Point3f temp = (vit + num2)->P() + tran/n;
				DIF += Distance(temp,vitd->P());
				vitd->P() = temp;
			}
		}
	}
}

void EditRepairPlugin::slotDeformationAuto()
{
	clock_t start,finish;
	double totaltime;
	start=clock();
	int i = 0;
	DIF = 1000;
	while (DIF > 0.2)
	{
		DIF = 0.0;
		slotDeformationNonPro();
		cout<<"slotDeformationNonPro "<<i<<" :"<<DIF<<endl;
		i++;
	}
	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n此程序的运行时间为"<<totaltime<<endl;
}

void EditRepairPlugin::slotDeformationNonPro2()
{
	CMeshO::VertexIterator vit;
	CVertexO* vitd;
	std::vector <vcg::Point3f> crownBoundary;
	vit = mesh_tem->cm.vert.begin();

	for (int i = 0; i < rootBoundaryT.size(); i++)
	{
		vcg::Point3f tran = vcg::Point3f(0.0,0.0,0.0);
		int n = 0;
		vitd = &mesh_tem_d->cm.vert[rootBoundaryT[i]];
		if(!(*vitd).IsD()){//判断该点是否存在
			CVertexO *vii = &(*vitd);
			vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
			CVertexO* firstV = pos.VFlip();//取和vitd相邻的第一个点
			CVertexO* tempV= NULL;
			crownBoundary.clear();
			do
			{
				pos.NextE();
				tempV = pos.VFlip();
				if (Color4f::Construct(tempV->C())== Color4<float>(Color4<float>::Red))
				{
					crownBoundary.push_back(tempV->P());
				}
				int num = tempV - &*(mesh_tem_d->cm.vert.begin());
				tran += (tempV->P() - (vit + num)->P());
				n++;
			}while (firstV != tempV);//从第一个点开始按顺序遍历，知道回到第一个点为止
			int num2 = &*vitd - &*(mesh_tem_d->cm.vert.begin());
			vcg::Point3f temp = (vit + num2)->P() + tran/n;
			if (crownBoundary.size()>0)
			{
				for ( int i = 0; i < crownBoundary.size();i++)
				{
					temp = temp + crownBoundary[i];
				}
				temp = temp / (crownBoundary.size()+1);
			}

			DIF += Distance(temp,vitd->P());
			vitd->P() = temp;
		}
	}
}

void EditRepairPlugin::slotDeformationNonPro3()
{
	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vitd;
	mesh_tem_d->updateDataMask(MeshModel::MM_VERTFACETOPO);
	mesh_tem_d->updateDataMask(MeshModel::MM_FACEFACETOPO);
	vit = mesh_tem->cm.vert.begin();
	for(vitd = mesh_tem_d->cm.vert.begin(); vitd != mesh_tem_d->cm.vert.end(); ++vitd){
		vcg::Color4f c1 = Color4f::Construct(vitd->C());
		if ( c1 != Color4<float>(Color4<float>::Red))//if red color?
		{
			vcg::Point3f btran = vcg::Point3f(0.0,0.0,0.0);
			vcg::Point3f tran = vcg::Point3f(0.0,0.0,0.0);
			float dis = 0.0;
			int n = 0;
			if(!(*vitd).IsD()){//判断该点是否存在
				CVertexO *vii = &(*vitd);
				vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
				CVertexO* firstV = pos.VFlip();//取和vitd相邻的第一个点
				CVertexO* tempV= NULL;
				vcg::Point3f dirv;
				do
				{
					pos.NextE();
					tempV = pos.VFlip();
					vcg::Color4f c2 = Color4f::Construct(tempV->C());
					if (c2 == Color4<float>(Color4<float>::Red))
					{
					}
					int num = tempV - &*(mesh_tem_d->cm.vert.begin());
					tran += (tempV->P() - (vit + num)->P());
					dirv = (vitd->P()-barycentre)/Distance(vitd->P(),barycentre);
					if ((tempV->P() - (vit + num)->P())*dirv>0)
					{
						dis += Distance(tempV->P(),(vit + num)->P());
					}else{
						dis -= Distance(tempV->P(),(vit + num)->P());
					}
					n++;
				}while (firstV != tempV);//从第一个点开始按顺序遍历，知道回到第一个点为止
				int num2 = &*vitd - &*(mesh_tem_d->cm.vert.begin());
				//vcg::Point3f temp = (vit + num2)->P() + tran/n;
				vcg::Point3f temp = (vit + num2)->P() + dirv*(dis/n);
				DIF += Distance(temp,vitd->P());
				vitd->P() = temp;
			}
		}
	}
}

void EditRepairPlugin::slotDeformationNonPro4()
{
	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vitd;
	std::vector <vcg::Point3f> crownBoundary;
	vit = mesh_tem->cm.vert.begin();

	for(vitd = mesh_tem_d->cm.vert.begin(); vitd != mesh_tem_d->cm.vert.end(); ++vitd){
		vcg::Color4f c1 = Color4f::Construct(vitd->C());
		if ( c1 != Color4<float>(Color4<float>::Red))//if white color?
		{
			vcg::Point3f tran = vcg::Point3f(0.0,0.0,0.0);
			int n = 0;
			if(!(*vitd).IsD()){//判断该点是否存在
				CVertexO *vii = &(*vitd);
				vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
				CVertexO* firstV = pos.VFlip();//取和vitd相邻的第一个点
				CVertexO* tempV= NULL;
				crownBoundary.clear();
				do
				{
					pos.NextE();
					tempV = pos.VFlip();
					if (Color4f::Construct(tempV->C())== Color4<float>(Color4<float>::Red))
					{
						crownBoundary.push_back(tempV->P());
					}
					int num = tempV - &*(mesh_tem_d->cm.vert.begin());
					tran += (tempV->P() - (vit + num)->P());
					n++;
				}while (firstV != tempV);//从第一个点开始按顺序遍历，知道回到第一个点为止
				int num2 = &*vitd - &*(mesh_tem_d->cm.vert.begin());
				vcg::Point3f temp = (vit + num2)->P() + tran/n;
				if (crownBoundary.size()>0)
				{
					for ( int i = 0; i < crownBoundary.size();i++)
					{
						temp = temp + crownBoundary[i];
					}
					temp = temp / (crownBoundary.size()+1);
				}

				DIF += Distance(temp,vitd->P());
				vitd->P() = temp;
			}
		}
	}
}


void EditRepairPlugin::slotDeformationAuto2()
{
	clock_t start,finish;
	double totaltime;
	start=clock();
	mesh_tem_d->updateDataMask(MeshModel::MM_VERTFACETOPO);
	mesh_tem_d->updateDataMask(MeshModel::MM_FACEFACETOPO);
	int i = 0;
	DIF = 1000;
	while (DIF > 0.2)
	{
		DIF = 0.0;
		slotDeformationNonPro2();
		cout<<"slotDeformationNonPro "<<i<<" :"<<DIF<<endl;
		i++;
	}
	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n此程序的运行时间为"<<totaltime<<endl;
}

void EditRepairPlugin::slotDeletePro()
{
	//mesh_tem_d = mesh_tem;
	CMeshO::VertexIterator vitd;
	CMeshO::FaceIterator fitd;
	vcg::Color4f c1;
	for(fitd = mesh_tem_d->cm.face.begin(); fitd != mesh_tem_d->cm.face.end(); ++fitd){
		if(!(*fitd).IsD()){//判断该点是否存在
			for (int i = 0;i<3;i++)
			{
				c1 = Color4f::Construct((fitd->V(i))->C());

				if(c1==Color4<float>(Color4<float>::Red)){
					fitd->SetS();
				} 
			}
		}
	}

	CMeshO::FaceIterator fi;
	CMeshO::VertexIterator vi;
	tri::UpdateSelection<CMeshO>::ClearVertex(mesh_tem_d->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceStrict(mesh_tem_d->cm);
	for(fi=mesh_tem_d->cm.face.begin();fi!=mesh_tem_d->cm.face.end();++fi)
		if(!(*fi).IsD() && (*fi).IsS() )
			tri::Allocator<CMeshO>::DeleteFace(mesh_tem_d->cm,*fi);
	for(vi=mesh_tem_d->cm.vert.begin();vi!=mesh_tem_d->cm.vert.end();++vi)
		if(!(*vi).IsD() && (*vi).IsS() )
			tri::Allocator<CMeshO>::DeleteVertex(mesh_tem_d->cm,*vi);
	mesh_tem_d->clearDataMask(MeshModel::MM_FACEFACETOPO);
	mesh_tem_d->clearDataMask(MeshModel::MM_FACEFLAGBORDER);
}

void EditRepairPlugin::slotImportSet()
{
	for(int i = 0; i < md->meshList.size(); i++){
		QString qstr = md->meshList[i]->shortName();
		char*  ch;
		QByteArray ba = qstr.toLatin1();
		ch=ba.data();
		if (ch[3] == '.')
		{
			mesh = md->meshList[i];
		}
	}
	QString dir = QFileDialog::getOpenFileName(this->gla, tr("Open Set File..."),
		currentfilepath,
		tr("Setting (*.set)"));

	int index = dir.lastIndexOf("/");
	currentfilepath = dir;
	currentfilepath.remove(index,dir.length()-index);

	std::string str = dir.toStdString();
	ifstream inputs;
	std::string oneline;
	inputs.open(str.c_str());
	if (inputs.is_open())
	{
		int i = -1;
		while(getline(inputs, oneline)){
			istringstream strStream(oneline);
			std::string kind,name;
			strStream >> kind;
			double ax,ay,az;
			if (kind == "#")
			{
				i = -1;
				strStream >>name;
				if (name.compare(mesh->shortName().toStdString())==0)
				{
					i = 0;
				}
			}
			if(kind == "m" && i != -1){
				strStream >> ax >> ay >> az;
				barycentre.X() = ax;
				barycentre.Y() = ay;
				barycentre.Z() = az;
			}
			if(kind == "x" && i != -1){
				strStream >> ax >> ay >> az;
				xcoord.X() = ax;
				xcoord.Y() = ay;
				xcoord.Z() = az;
			}
			if(kind == "y" && i != -1){
				strStream >> ax >> ay >> az;
				ycoord.X() = ax;
				ycoord.Y() = ay;
				ycoord.Z() = az;
			}
			if(kind == "z" && i != -1){
				strStream >> ax >> ay >> az;
				zcoord.X() = ax;
				zcoord.Y() = ay;
				zcoord.Z() = az;
			}
		}
		inputs.close();

		float max = -100000;
		float min = 100000;
		for (int i = 0 ; i < mesh->cm.vert.size(); i++)
		{
			float zpro = (mesh->cm.vert[i].P() - barycentre)* zcoord;
			if (zpro > max)
			{
				max = zpro;
			}
			if (zpro < min)
			{
				min = zpro;
			}
		}
		stagesDialog->ui.toothinfo->setText(QString("Max: ")+QString("%1").arg(max)+QString("; Min: ") + QString("%1").arg(min));
		stagesDialog->ui.toothinfo->adjustSize();


		if(CurrentMode() != TranMode::IMPORTCOORD)
			SetMode(TranMode::IMPORTCOORD);
		importset = true;
	}
}

void EditRepairPlugin::slotFindBoundary()
{
	mesh->updateDataMask(MeshModel::MM_FACEFLAGBORDER);
	borders.clear();
	CMeshO::VertexIterator vi;
	int i = 0;
	for(vi = mesh->cm.vert.begin(); vi != mesh->cm.vert.end(); ++vi){
		if( !(*vi).IsD() )
		{
			if((*vi).IsB() )
			{
				borders.push_back(i);
			}
		}
		i++;
	}
	cout<<"slotFindBoundary1"<<endl;
	mesh_tem_d = mesh2;//test
	mesh_tem_d->updateDataMask(MeshModel::MM_FACEFLAGBORDER);
	borders_root.clear();
	i = 0;
	for(vi = mesh_tem_d->cm.vert.begin(); vi != mesh_tem_d->cm.vert.end(); ++vi){
		if( !(*vi).IsD() )
		{
			if((*vi).IsB() )
			{
				borders_root.push_back(i);
			}
		}
		i++;
	}
	cout<<"slotFindBoundary2"<<endl;
}

void EditRepairPlugin::slotFindBoundary1()
{
	mesh->updateDataMask(MeshModel::MM_FACEFLAGBORDER);
	borders.clear();
	CMeshO::VertexIterator vi;
	int i = 0;
	for(vi = mesh->cm.vert.begin(); vi != mesh->cm.vert.end(); ++vi){
		if( !(*vi).IsD() )
		{
			if((*vi).IsB() )
			{
				borders.push_back(i);
			}
		}
		i++;
	}
	//mesh_tem_d = mesh_tem;//test
	mesh_tem_d->updateDataMask(MeshModel::MM_FACEFLAGBORDER);
	borders_root.clear();
	i = 0;
	for(vi = mesh_tem_d->cm.vert.begin(); vi != mesh_tem_d->cm.vert.end(); ++vi){
		if( !(*vi).IsD() )
		{
			if((*vi).IsB() )
			{
				borders_root.push_back(i);
			}
		}
		i++;
	}
}


void EditRepairPlugin::slotCopyMesh()
{
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	std::vector<CVertexO> vv_mm;
	std::vector<CFaceO> vf_mm;
	std::vector<int> VertexId(mesh_tem_d->cm.vert.size());
	int numvert = 0;

	for (vi=mesh_tem_d->cm.vert.begin();vi!=mesh_tem_d->cm.vert.end();vi++){
		if (!(*vi).IsD())
		{
			VertexId[vi-mesh_tem_d->cm.vert.begin()]=numvert;
			vv_mm.push_back(*vi);
			numvert++;
		}
	}
	for(fi=mesh_tem_d->cm.face.begin(); fi!=mesh_tem_d->cm.face.end(); ++fi)
	{
		if( !(*fi).IsD() ){
			vf_mm.push_back(*fi);
		}
	}
	//for (int i = 0; i<mesh_tem_d->cm.face.size();i++){
	//	if (!mesh_tem_d->cm.face[i].IsD())
	//	{
	//		vf_mm.push_back(mesh_tem_d->cm.face[i]);
	//	}
	//}

	mesh_copy = md->addNewMesh("mesh_copy","",false);

	mesh_copy->cm.Clear();
	int VN,FN;
	VN=vv_mm.size();
	FN=vf_mm.size();
	vcg::tri::Allocator<CMeshO>::AddVertices(mesh_copy->cm,VN);
	vcg::tri::Allocator<CMeshO>::AddFaces(mesh_copy->cm,FN);
	CMeshO::VertexPointer *ivp = new CMeshO::VertexPointer[VN];
	int cnt=0;
	vi= mesh_copy->cm.vert.begin();
	std::vector<CVertexO>::const_iterator itv;
	for(itv=vv_mm.begin();itv!=vv_mm.end();itv++){
		ivp[cnt]=&*vi;
		(*vi).P()= (*itv).P();
		++vi;
		++cnt;
	}
	fi=mesh_copy->cm.face.begin();
	std::vector<CFaceO>::const_iterator itf;
	for (itf = vf_mm.begin(); itf!= vf_mm.end(); itf++)
	{
		for(int k=0;k<(*itf).VN();k++)
		{
			int vInd = VertexId[(*itf).V(k)-&*(mesh_tem_d->cm.vert.begin())];
			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}

	vcg::tri::UpdateBounding<CMeshO>::Box(mesh_copy->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_copy->cm);
	mesh_tem_d = mesh_copy;
}


void EditRepairPlugin::slotSortBorder()
{
	float nearestProFront = 1000.0;
	int frontP = 0;
	for (int j = 0; j<borders.size(); j++)
	{
		Point3f p1 = mesh->cm.vert[borders[j]].P()-barycentre;
		float xPro = p1*xcoord;
		float yPro = p1*ycoord;
		float zPro = p1*zcoord;
		int id =borders[j];
		if (xPro>=0 && fabs(yPro)<nearestProFront)
		{
			nearestProFront = fabs(yPro);
			frontP = id;
		}
	}
	cout<<"slotSortBorder1"<<endl;
	mesh->updateDataMask(MeshModel::MM_VERTFACETOPO);
	mesh->updateDataMask(MeshModel::MM_FACEFACETOPO);

	teethBordersInOrder.clear();
	teethBordersInOrder.push_back(frontP);

	int findtheFirst = 0;
	CVertexO* V1= NULL;

	CVertexO *vi = &(mesh->cm.vert[frontP]);
	vcg::face::JumpingPos<CFaceO> pos(vi->VFp(),vi);
	CVertexO* firstV = pos.VFlip();
	CVertexO* tempV= NULL;
	CVertexO* tempVLast= NULL;
	CVertexO* tempVLastLast= NULL;
	Point3f p1,p2;
	do 
	{
		pos.NextE();
		tempV = pos.VFlip();
		// 			cout<<"vi:"<< vi->P().X()<< vi->P().Y()<< vi->P().Z()<<endl;
		// 			cout<<"firstV:"<< firstV->P().X()<< firstV->P().Y()<< firstV->P().Z()<<endl;
		// 			cout<<"tempV:"<< tempV->P().X()<< tempV->P().Y()<< tempV->P().Z()<<endl;
		if (tempVLast != NULL && tempVLastLast != NULL && tempVLast->IsB()){
			Point3f ptemp = tempVLast->P()-vi->P();
			if (ptemp*(ycoord)<0)//沿Y轴负方向
			{
				if (tempV == tempVLastLast)
				{
					break;
					//cout<<"slotSortBorder2"<<endl;
				}
			}
		}
		tempVLastLast = tempVLast;
		tempVLast = tempV;
	} while (1);

	cout<<"slotSortBorder3"<<endl;

	CVertexO* Vp= vi;

	int num = 1;
	while (tempVLast != vi)
	{
		teethBordersInOrder.push_back(tempVLast-&*(mesh->cm.vert.begin()));
		CVertexO *vi0 = tempVLast;
		vcg::face::JumpingPos<CFaceO> pos0(vi0->VFp(),vi0);
		CVertexO* firstV0 = pos0.VFlip();
		CVertexO* tempV0= NULL;//访问V的一环邻域
		CVertexO* tempV0Last= NULL;
		CVertexO* tempV0LastLast= NULL;
		do 
		{

			pos0.NextE();
			tempV0 = pos0.VFlip();
			if (tempV0Last != NULL && tempV0LastLast != NULL && tempV0Last->IsB() && tempV0Last != Vp)
			{
				if (tempV0 == tempV0LastLast)
				{
					break;
				}
			}
			tempV0LastLast = tempV0Last;
			tempV0Last = tempV0;
		} while (1);
		Vp = tempVLast;
		tempVLast = tempV0Last;
		num++;
		if (num>1000)
		{
			break;
		}
	}

	cout<<"slotSortBorder4"<<endl;

	//root
	vcg::tri::UpdateBounding<CMeshO>::Box(mesh_tem_d->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_tem_d->cm);
	mesh_tem_d->clearDataMask(MeshModel::MM_FACEFACETOPO);
	mesh_tem_d->clearDataMask(MeshModel::MM_FACEFLAGBORDER);
	mesh_tem_d->updateDataMask(MeshModel::MM_VERTFACETOPO);
	mesh_tem_d->updateDataMask(MeshModel::MM_FACEFACETOPO);
	mesh_tem_d->updateDataMask(MeshModel::MM_FACEFLAGBORDER);

	nearestProFront = -1000.0;
	frontP = 0;
	for (int j = 0; j<borders_root.size(); j++)
	{
		Point3f p1 = mesh_tem_d->cm.vert[borders_root[j]].P()-barycentre;
		float xPro = p1*xcoord;
		float yPro = p1*ycoord;
		float zPro = p1*zcoord;
		int id =borders_root[j];
		if (xPro>nearestProFront)
		{
			nearestProFront = xPro;
			frontP = id;
		}
	}
	//mesh_tem_d->updateDataMask(MeshModel::MM_VERTFACETOPO);
	//mesh_tem_d->updateDataMask(MeshModel::MM_FACEFACETOPO);

	rootBordersInOrder.clear();
	rootBordersInOrder.push_back(frontP);

	findtheFirst = 0;
	V1= NULL;

	vi = &(mesh_tem_d->cm.vert[frontP]);
	vcg::face::JumpingPos<CFaceO> pos1(vi->VFp(),vi);
	firstV = pos1.VFlip();
	tempV= NULL;
	tempVLast= NULL;
	tempVLastLast= NULL;
	p1,p2;
	do 
	{
		pos1.NextE();
		tempV = pos1.VFlip();
		if (tempVLast != NULL && tempVLastLast != NULL && tempVLast->IsB()){
			Point3f ptemp = tempVLast->P()-vi->P();
			if (ptemp*(zcoord)>0)
			{
				if (tempV == tempVLastLast)
				{
					break;
				}
			}
		}
		tempVLastLast = tempVLast;
		tempVLast = tempV;
	} while (1);

	cout<<"slotSortBorder5"<<endl;//从这往下有问题

	Vp= vi;

	num = 1;
	while (tempVLast != vi)
	{
		//cout<<tempVLast-&*(mesh_tem_d->cm.vert.begin())<<endl;
		rootBordersInOrder.push_back(tempVLast-&*(mesh_tem_d->cm.vert.begin()));
		CVertexO *vi0 = tempVLast;
		vcg::face::JumpingPos<CFaceO> pos0(vi0->VFp(),vi0);
		CVertexO* firstV0 = pos0.VFlip();
		CVertexO* tempV0= NULL;//访问V的一环邻域
		CVertexO* tempV0Last= NULL;
		CVertexO* tempV0LastLast= NULL;
		do 
		{

			pos0.NextE();
			tempV0 = pos0.VFlip();
			if (tempV0Last != NULL && tempV0LastLast != NULL && tempV0Last->IsB() && tempV0Last != Vp)
			{
				if (tempV0 == tempV0LastLast)
				{
					//cout<<"if (tempV0 == tempV0LastLast)"<<endl;
					break;
				}
			}
			tempV0LastLast = tempV0Last;
			tempV0Last = tempV0;
		} while (1);
		Vp = tempVLast;
		tempVLast = tempV0Last;
		num++;
		if (num>1000)
		{
			break;
		}
		//cout<<num<<endl;
		//break;
	}
	cout<<"slotSortBorder6"<<endl;
}

void EditRepairPlugin::slotSortBorder2()
{
	float nearestProFront = 1000.0;
	int frontP = 0;
	for (int j = 0; j<borders.size(); j++)
	{
		Point3f p1 = mesh->cm.vert[borders[j]].P()-barycentre;
		float xPro = p1*xcoord;
		float yPro = p1*ycoord;
		float zPro = p1*zcoord;
		int id =borders[j];
		if (xPro>=0 && fabs(yPro)<nearestProFront)
		{
			nearestProFront = fabs(yPro);
			frontP = id;
		}
	}
	//cout<<"slotSortBorder1"<<endl;
	mesh->updateDataMask(MeshModel::MM_VERTFACETOPO);
	mesh->updateDataMask(MeshModel::MM_FACEFACETOPO);

	teethBordersInOrder.clear();
	teethBordersInOrder.push_back(frontP);

	int findtheFirst = 0;
	CVertexO* V1= NULL;

	CVertexO *vi = &(mesh->cm.vert[frontP]);
	vcg::face::JumpingPos<CFaceO> pos(vi->VFp(),vi);
	CVertexO* firstV = pos.VFlip();
	CVertexO* tempV= NULL;
	CVertexO* tempVLast= NULL;
	CVertexO* tempVLastLast= NULL;
	Point3f p1,p2;
	do 
	{
		pos.NextE();
		tempV = pos.VFlip();
		// 			cout<<"vi:"<< vi->P().X()<< vi->P().Y()<< vi->P().Z()<<endl;
		// 			cout<<"firstV:"<< firstV->P().X()<< firstV->P().Y()<< firstV->P().Z()<<endl;
		// 			cout<<"tempV:"<< tempV->P().X()<< tempV->P().Y()<< tempV->P().Z()<<endl;
		if (tempVLast != NULL && tempVLastLast != NULL && tempVLast->IsB()){
			Point3f ptemp = tempVLast->P()-vi->P();
			if (ptemp*(ycoord)<0)//沿Y轴负方向
			{
				if (tempV == tempVLastLast)
				{
					break;
					//cout<<"slotSortBorder2"<<endl;
				}
			}
		}
		tempVLastLast = tempVLast;
		tempVLast = tempV;
	} while (1);

	//cout<<"slotSortBorder3"<<endl;

	CVertexO* Vp= vi;

	int num = 1;
	while (tempVLast != vi)
	{
		teethBordersInOrder.push_back(tempVLast-&*(mesh->cm.vert.begin()));
		CVertexO *vi0 = tempVLast;
		vcg::face::JumpingPos<CFaceO> pos0(vi0->VFp(),vi0);
		CVertexO* firstV0 = pos0.VFlip();
		CVertexO* tempV0= NULL;//访问V的一环邻域
		CVertexO* tempV0Last= NULL;
		CVertexO* tempV0LastLast= NULL;
		do 
		{

			pos0.NextE();
			tempV0 = pos0.VFlip();
			if (tempV0Last != NULL && tempV0LastLast != NULL && tempV0Last->IsB() && tempV0Last != Vp)
			{
				if (tempV0 == tempV0LastLast)
				{
					break;
				}
			}
			tempV0LastLast = tempV0Last;
			tempV0Last = tempV0;
		} while (1);
		Vp = tempVLast;
		tempVLast = tempV0Last;
		num++;
		if (num>1000)
		{
			break;
		}
	}

	//cout<<"slotSortBorder4"<<endl;

	//root
	vcg::tri::UpdateBounding<CMeshO>::Box(mesh_tem->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_tem->cm);
	mesh_tem->clearDataMask(MeshModel::MM_FACEFACETOPO);
	mesh_tem->clearDataMask(MeshModel::MM_FACEFLAGBORDER);
	mesh_tem->updateDataMask(MeshModel::MM_VERTFACETOPO);
	mesh_tem->updateDataMask(MeshModel::MM_FACEFACETOPO);
	mesh_tem->updateDataMask(MeshModel::MM_FACEFLAGBORDER);

	nearestProFront = 1000.0;
	frontP = 0;
	for (int j = 0; j<borders_root.size(); j++)
	{
		Point3f p1 = mesh_tem->cm.vert[borders_root[j]].P()-barycentre;
		float xPro = p1*xcoord;
		float yPro = p1*ycoord;
		float zPro = p1*zcoord;
		int id =borders_root[j];
		if (xPro>=0 && fabs(yPro)<nearestProFront)
		{
			nearestProFront = fabs(yPro);
			frontP = id;
		}
	}
	//mesh_tem_d->updateDataMask(MeshModel::MM_VERTFACETOPO);
	//mesh_tem_d->updateDataMask(MeshModel::MM_FACEFACETOPO);

	rootBordersInOrder.clear();
	rootBordersInOrder.push_back(frontP);

	findtheFirst = 0;
	V1= NULL;

	vi = &(mesh_tem->cm.vert[frontP]);
	vcg::face::JumpingPos<CFaceO> pos1(vi->VFp(),vi);
	firstV = pos1.VFlip();
	tempV= NULL;
	tempVLast= NULL;
	tempVLastLast= NULL;
	p1,p2;
	do 
	{
		pos1.NextE();
		tempV = pos1.VFlip();
		if (tempVLast != NULL && tempVLastLast != NULL && tempVLast->IsB()){
			Point3f ptemp = tempVLast->P()-vi->P();
			if (ptemp*(ycoord)<0)
			{
				if (tempV == tempVLastLast)
				{
					break;
				}
			}
		}
		tempVLastLast = tempVLast;
		tempVLast = tempV;
	} while (1);

	//cout<<"slotSortBorder5"<<endl;//从这往下有问题

	Vp= vi;

	num = 1;
	while (tempVLast != vi)
	{
		//cout<<tempVLast-&*(mesh_tem_d->cm.vert.begin())<<endl;
		rootBordersInOrder.push_back(tempVLast-&*(mesh_tem->cm.vert.begin()));
		CVertexO *vi0 = tempVLast;
		vcg::face::JumpingPos<CFaceO> pos0(vi0->VFp(),vi0);
		CVertexO* firstV0 = pos0.VFlip();
		CVertexO* tempV0= NULL;//访问V的一环邻域
		CVertexO* tempV0Last= NULL;
		CVertexO* tempV0LastLast= NULL;
		do 
		{

			pos0.NextE();
			tempV0 = pos0.VFlip();
			if (tempV0Last != NULL && tempV0LastLast != NULL && tempV0Last->IsB() && tempV0Last != Vp)
			{
				if (tempV0 == tempV0LastLast)
				{
					//cout<<"if (tempV0 == tempV0LastLast)"<<endl;
					break;
				}
			}
			tempV0LastLast = tempV0Last;
			tempV0Last = tempV0;
		} while (1);
		Vp = tempVLast;
		tempVLast = tempV0Last;
		num++;
		if (num>1000)
		{
			break;
		}
		//cout<<num<<endl;
		//break;
	}
	//cout<<"slotSortBorder6"<<endl;
}

void EditRepairPlugin::FindNeastestP()
{
	float dis = 10000;
	int flag = 0;
	int i;
	for (i = 0; i < teethBordersInOrder.size(); i++)
	{
		if (Distance(mesh->cm.vert[teethBordersInOrder[i]].P(), mesh_tem_d->cm.vert[rootBordersInOrder[0]].P()) < dis)
		{
			dis = Distance(mesh->cm.vert[teethBordersInOrder[i]].P(), mesh_tem_d->cm.vert[rootBordersInOrder[0]].P());
			flag = i;
		}
	}
	int maxflag = 0;
	for (i = flag ;i < teethBordersInOrder.size(); i++)
	{
		float mindis = 10000;
		int minflag = 0;
		for (int j = 0; j < rootBordersInOrder.size(); j++)
		{
			if (Distance(mesh->cm.vert[teethBordersInOrder[i]].P(), mesh_tem_d->cm.vert[rootBordersInOrder[j]].P()) < mindis)
			{
				mindis = Distance(mesh->cm.vert[teethBordersInOrder[i]].P(), mesh_tem_d->cm.vert[rootBordersInOrder[j]].P());
				minflag = j;
			}
		}
		if (mindis < 0.5)//can change
		{
			if (minflag > maxflag)
			{
				maxflag = minflag;
			}
		}else{
			break;
		}
	}
	teethBordersInOrderStart = flag; 
	teethBordersInOrderEnd = i - 1; 
	rootBordersInOrderEnd = maxflag; 

}

void EditRepairPlugin::slotFillGap()
{
	FindNeastestP();
	getSE = true;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	std::vector<CVertexO> vv_mm;
	std::vector<CFaceO> vf_mm;

	for (vi=mesh->cm.vert.begin();vi!=mesh->cm.vert.end();vi++){
		vv_mm.push_back(*vi);
	}
	for (vi=mesh_tem_d->cm.vert.begin();vi!=mesh_tem_d->cm.vert.end();vi++){
		vv_mm.push_back(*vi);
	}
	for (int i = 0; i<mesh->cm.face.size();i++){
		vf_mm.push_back(mesh->cm.face[i]);
	}
	for (int i = 0; i<mesh_tem_d->cm.face.size();i++){
		vf_mm.push_back(mesh_tem_d->cm.face[i]);
	}
	repaired_mesh = md->addNewMesh("repaired_mesh","",false);

	repaired_mesh->cm.Clear();
	int VN,FN;
	VN=vv_mm.size();
	FN=vf_mm.size()+teethBordersInOrderEnd - teethBordersInOrderStart + rootBordersInOrderEnd-1;

	//FN=vf_mm.size()+teethBordersInOrder.size()+rootBordersInOrder.size();
	vcg::tri::Allocator<CMeshO>::AddVertices(repaired_mesh->cm,VN);
	vcg::tri::Allocator<CMeshO>::AddFaces(repaired_mesh->cm,FN);
	CMeshO::VertexPointer *ivp = new CMeshO::VertexPointer[VN];
	int cnt=0;
	vi= repaired_mesh->cm.vert.begin();
	std::vector<CVertexO>::const_iterator itv;
	for(itv=vv_mm.begin();itv!=vv_mm.end();itv++){
		ivp[cnt]=&*vi;
		(*vi).P()= (*itv).P();
		++vi;
		++cnt;
	}
	fi=repaired_mesh->cm.face.begin();
	std::vector<CFaceO>::const_iterator itf;
	for(int i = 0; i<mesh->cm.face.size(); i++){
		for (int k = 0 ;k<3;k++)
		{
			int vInd = vf_mm[i].V(k)-&*(mesh->cm.vert.begin());

			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}
	int f1 = mesh->cm.face.size();
	int v1 = mesh->cm.vert.size();
	for(int i = f1; i<vf_mm.size(); i++){
		for (int k = 0 ;k<3;k++)
		{
			int vInd = vf_mm[i].V(k)-&*(mesh_tem_d->cm.vert.begin())+v1;
			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}



	int time = vf_mm.size();
	// 	int vInd;
	// 	vInd = teethBordersInOrder[0];
	// 	(*fi).V(0)= ivp[vInd];
	// 	vInd = teethBordersInOrder[1];
	// 	(*fi).V(1)= ivp[vInd];
	// 	vInd = rootBordersInOrder[0]+v1;
	// 	(*fi).V(2)= ivp[vInd];
	// 	++fi;
	// 	time++;

	for (int i = teethBordersInOrderStart; i < teethBordersInOrderEnd ; i++)
	{
		int id1 = i;
		int id2 = (i+1)%teethBordersInOrder.size();
		float dis1 = Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[0]].P());
		float dis2 = Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[0]].P());
		int r1 = 0;
		int r2 = 0;
		for (int j = 1; j< rootBordersInOrderEnd+1; j++)
		{
			if (Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[j]].P())<dis1)
			{
				dis1 = Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[j]].P());
				r1 = j;
			}
			if (Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[j]].P())<dis2)
			{
				dis2 = Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[j]].P());
				r2 = j;
			}
		}

		int vInd;
		vInd = teethBordersInOrder[id1];
		(*fi).V(0)= ivp[vInd];
		vInd = teethBordersInOrder[id2];
		(*fi).V(1)= ivp[vInd];
		vInd = rootBordersInOrder[r1]+v1;
		(*fi).V(2)= ivp[vInd];
		++fi;
		time++;

		if (r1 < r2)
		{
			for (int j = r1; j<r2; j++)
			{
				vInd = rootBordersInOrder[j%(rootBordersInOrder.size())]+v1;
				(*fi).V(0)= ivp[vInd];
				vInd = teethBordersInOrder[id2];
				(*fi).V(1)= ivp[vInd];
				vInd = rootBordersInOrder[(j+1)%(rootBordersInOrder.size())]+v1;
				(*fi).V(2)= ivp[vInd];
				++fi;
				time++;
			}
		}else if (r1 > r2)
		{
			for (int j = r1; j<r2+rootBordersInOrder.size(); j++)
			{
				vInd = rootBordersInOrder[j%(rootBordersInOrder.size())]+v1;
				(*fi).V(0)= ivp[vInd];
				vInd = teethBordersInOrder[id2];
				(*fi).V(1)= ivp[vInd];
				vInd = rootBordersInOrder[(j+1)%(rootBordersInOrder.size())]+v1;
				(*fi).V(2)= ivp[vInd];
				++fi;
				time++;
			}
		}
	}

	cout<<"Success!"<<time<<" "<<FN<<endl;

	vcg::tri::UpdateBounding<CMeshO>::Box(repaired_mesh->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(repaired_mesh->cm);

}

void EditRepairPlugin::slotFillGap1()
{
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	std::vector<CVertexO> vv_mm;
	std::vector<CFaceO> vf_mm;

	for (vi=mesh->cm.vert.begin();vi!=mesh->cm.vert.end();vi++){
		vv_mm.push_back(*vi);
	}
	for (vi=mesh_tem_d->cm.vert.begin();vi!=mesh_tem_d->cm.vert.end();vi++){
		vv_mm.push_back(*vi);
	}
	for (int i = 0; i<mesh->cm.face.size();i++){
		vf_mm.push_back(mesh->cm.face[i]);
	}
	for (int i = 0; i<mesh_tem_d->cm.face.size();i++){
		vf_mm.push_back(mesh_tem_d->cm.face[i]);
	}
	repaired_mesh = md->addNewMesh("repaired_mesh","",false);

	repaired_mesh->cm.Clear();
	int VN,FN;
	VN=vv_mm.size();
	FN=vf_mm.size()+teethBordersInOrder.size()+rootBordersInOrder.size();
	vcg::tri::Allocator<CMeshO>::AddVertices(repaired_mesh->cm,VN);
	vcg::tri::Allocator<CMeshO>::AddFaces(repaired_mesh->cm,FN);
	CMeshO::VertexPointer *ivp = new CMeshO::VertexPointer[VN];
	int cnt=0;
	vi= repaired_mesh->cm.vert.begin();
	std::vector<CVertexO>::const_iterator itv;
	for(itv=vv_mm.begin();itv!=vv_mm.end();itv++){
		ivp[cnt]=&*vi;
		(*vi).P()= (*itv).P();
		++vi;
		++cnt;
	}
	fi=repaired_mesh->cm.face.begin();
	std::vector<CFaceO>::const_iterator itf;
	for(int i = 0; i<mesh->cm.face.size(); i++){
		for (int k = 0 ;k<3;k++)
		{
			int vInd = vf_mm[i].V(k)-&*(mesh->cm.vert.begin());

			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}
	int f1 = mesh->cm.face.size();
	int v1 = mesh->cm.vert.size();
	for(int i = f1; i<vf_mm.size(); i++){
		for (int k = 0 ;k<3;k++)
		{
			int vInd = vf_mm[i].V(k)-&*(mesh_tem_d->cm.vert.begin())+v1;
			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}



	int time = vf_mm.size();
	// 	int vInd;
	// 	vInd = teethBordersInOrder[0];
	// 	(*fi).V(0)= ivp[vInd];
	// 	vInd = teethBordersInOrder[1];
	// 	(*fi).V(1)= ivp[vInd];
	// 	vInd = rootBordersInOrder[0]+v1;
	// 	(*fi).V(2)= ivp[vInd];
	// 	++fi;
	// 	time++;

	for (int i = 0; i < teethBordersInOrder.size(); i++)
	{
		int id1 = i;
		int id2 = (i+1)%teethBordersInOrder.size();
		float dis1 = Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[0]].P());
		float dis2 = Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[0]].P());
		int r1 = 0;
		int r2 = 0;
		for (int j = 1; j< rootBordersInOrder.size(); j++)
		{
			if (Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[j]].P())<dis1)
			{
				dis1 = Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[j]].P());
				r1 = j;
			}
			if (Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[j]].P())<dis2)
			{
				dis2 = Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh_tem_d->cm.vert[rootBordersInOrder[j]].P());
				r2 = j;
			}
		}

		int vInd;
		vInd = teethBordersInOrder[id1];
		(*fi).V(0)= ivp[vInd];
		vInd = teethBordersInOrder[id2];
		(*fi).V(1)= ivp[vInd];
		vInd = rootBordersInOrder[r1]+v1;
		(*fi).V(2)= ivp[vInd];
		++fi;
		time++;

		if (r1 < r2)
		{
			for (int j = r1; j<r2; j++)
			{
				vInd = rootBordersInOrder[j%(rootBordersInOrder.size())]+v1;
				(*fi).V(0)= ivp[vInd];
				vInd = teethBordersInOrder[id2];
				(*fi).V(1)= ivp[vInd];
				vInd = rootBordersInOrder[(j+1)%(rootBordersInOrder.size())]+v1;
				(*fi).V(2)= ivp[vInd];
				++fi;
				time++;
			}
		}else if (r1 > r2)
		{
			for (int j = r1; j<r2+rootBordersInOrder.size(); j++)
			{
				vInd = rootBordersInOrder[j%(rootBordersInOrder.size())]+v1;
				(*fi).V(0)= ivp[vInd];
				vInd = teethBordersInOrder[id2];
				(*fi).V(1)= ivp[vInd];
				vInd = rootBordersInOrder[(j+1)%(rootBordersInOrder.size())]+v1;
				(*fi).V(2)= ivp[vInd];
				++fi;
				time++;
			}
		}
	}

	cout<<"Success!"<<time<<" "<<FN<<endl;

	vcg::tri::UpdateBounding<CMeshO>::Box(repaired_mesh->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(repaired_mesh->cm);
}

void EditRepairPlugin::slotFillGap2()
{
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	std::vector<CVertexO> vv_mm;
	std::vector<CFaceO> vf_mm;

	for (vi=mesh->cm.vert.begin();vi!=mesh->cm.vert.end();vi++){
		vv_mm.push_back(*vi);
	}
	for (vi=mesh_tem->cm.vert.begin();vi!=mesh_tem->cm.vert.end();vi++){
		vv_mm.push_back(*vi);
	}
	for (int i = 0; i<mesh->cm.face.size();i++){
		vf_mm.push_back(mesh->cm.face[i]);
	}
	for (int i = 0; i<mesh_tem->cm.face.size();i++){
		vf_mm.push_back(mesh_tem->cm.face[i]);
	}

	MeshModel *newm;
	newm =  md->addNewMesh(mesh_tem->fullName(),"",false);

	newm->cm.Clear();
	int VN,FN;
	VN=vv_mm.size();
	FN=vf_mm.size()+teethBordersInOrder.size()+rootBordersInOrder.size();
	vcg::tri::Allocator<CMeshO>::AddVertices(newm->cm,VN);
	vcg::tri::Allocator<CMeshO>::AddFaces(newm->cm,FN);
	CMeshO::VertexPointer *ivp = new CMeshO::VertexPointer[VN];
	int cnt=0;
	vi= newm->cm.vert.begin();
	std::vector<CVertexO>::const_iterator itv;
	for(itv=vv_mm.begin();itv!=vv_mm.end();itv++){
		ivp[cnt]=&*vi;
		(*vi).P()= (*itv).P();
		++vi;
		++cnt;
	}
	fi=newm->cm.face.begin();
	std::vector<CFaceO>::const_iterator itf;
	for(int i = 0; i<mesh->cm.face.size(); i++){
		for (int k = 0 ;k<3;k++)
		{
			int vInd = vf_mm[i].V(k)-&*(mesh->cm.vert.begin());

			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}
	int f1 = mesh->cm.face.size();
	int v1 = mesh->cm.vert.size();
	for(int i = f1; i<vf_mm.size(); i++){
		for (int k = 0 ;k<3;k++)
		{
			int vInd = vf_mm[i].V(k)-&*(mesh_tem->cm.vert.begin())+v1;
			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}

	for (int i = 0; i < teethBordersInOrder.size(); i++)
	{
		int id1 = i;
		int id2 = (i+1)%teethBordersInOrder.size();
		float dis1 = Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh_tem->cm.vert[rootBordersInOrder[0]].P());
		float dis2 = Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh_tem->cm.vert[rootBordersInOrder[0]].P());
		int r1 = 0;
		int r2 = 0;
		for (int j = 1; j< rootBordersInOrder.size(); j++)
		{
			if (Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh_tem->cm.vert[rootBordersInOrder[j]].P())<dis1)
			{
				dis1 = Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh_tem->cm.vert[rootBordersInOrder[j]].P());
				r1 = j;
			}
			if (Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh_tem->cm.vert[rootBordersInOrder[j]].P())<dis2)
			{
				dis2 = Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh_tem->cm.vert[rootBordersInOrder[j]].P());
				r2 = j;
			}
		}

		int vInd;
		vInd = teethBordersInOrder[id1];
		(*fi).V(0)= ivp[vInd];
		vInd = teethBordersInOrder[id2];
		(*fi).V(1)= ivp[vInd];
		vInd = rootBordersInOrder[r1]+v1;
		(*fi).V(2)= ivp[vInd];
		++fi;

		if (r1 < r2)
		{
			for (int j = r1; j<r2; j++)
			{
				vInd = rootBordersInOrder[j%(rootBordersInOrder.size())]+v1;
				(*fi).V(0)= ivp[vInd];
				vInd = teethBordersInOrder[id2];
				(*fi).V(1)= ivp[vInd];
				vInd = rootBordersInOrder[(j+1)%(rootBordersInOrder.size())]+v1;
				(*fi).V(2)= ivp[vInd];
				++fi;
			}
		}else if (r1 > r2)
		{
			for (int j = r1; j<r2+rootBordersInOrder.size(); j++)
			{
				vInd = rootBordersInOrder[j%(rootBordersInOrder.size())]+v1;
				(*fi).V(0)= ivp[vInd];
				vInd = teethBordersInOrder[id2];
				(*fi).V(1)= ivp[vInd];
				vInd = rootBordersInOrder[(j+1)%(rootBordersInOrder.size())]+v1;
				(*fi).V(2)= ivp[vInd];
				++fi;
			}
		}
	}

	vcg::tri::UpdateBounding<CMeshO>::Box(newm->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(newm->cm);
	md->delMesh(mesh_tem);

	newm->clearDataMask(MeshModel::MM_FACEMARK);
	newm->updateDataMask(MeshModel::MM_FACEMARK);
	newm->updateDataMask(MeshModel::MM_FACEFACETOPO);

	mesh_tem = newm;
}

void EditRepairPlugin::sortTeeth()//Sort TeethL from L to R and Sort TeethU from R to L（每次只能单独处理上牙或者下牙）
{
	vector<IdMesh>::iterator it;
	MeshSet.clear();
	foreach(MeshModel *mm, md->meshList)
	{
		tri::UpdateBounding<CMeshO>::Box(mm->cm);
		mm->clearDataMask(MeshModel::MM_FACEMARK);
		mm->updateDataMask(MeshModel::MM_VERTFACETOPO | MeshModel::MM_FACEMARK|MeshModel::MM_VERTMARK|MeshModel::MM_FACEFACETOPO);

		if (!mm->hasDataMask(MeshModel::MM_VERTCOLOR))
		{
			mm->updateDataMask(MeshModel::MM_VERTCOLOR);
			//tri::UpdateColor<CMeshO>::VertexConstant(mm->cm,Color4b(150, 150, 150, 255));
		}
		tri::InitFaceIMark(mm->cm);
		tri::InitVertexIMark(mm->cm);

		QString qstr = mm->shortName();
		char*  ch;
		QByteArray ba = qstr.toLatin1();
		ch=ba.data();
		/*string name = string (ch);*/
		if (ch[0] == 'L')
		{
			if (ch[1]=='L') //左侧
			{
				char numc = ch[2];
				int num = 10 - (numc-'0');
				IdMesh im;
				im.id = num;
				im.mm = mm;
				im.haveBottom = false;
				im.deleteFn = 0;
				int insertflag = 0;
				for (it = MeshSet.begin();it != MeshSet.end();it++)
				{
					if (it->id>num)
					{	
						MeshSet.insert(it,im);
						insertflag++;
						break;
					}
				}
				if(insertflag == 0){
					MeshSet.push_back(im);
				}
			}else if (ch[1]=='R') //右侧
			{
				char numc = ch[2];
				int num = 10 + (numc-'0');
				IdMesh im;
				im.id = num;
				im.mm = mm;	
				im.haveBottom = false;
				im.deleteFn = 0;
				int insertflag = 0;
				for (it = MeshSet.begin();it != MeshSet.end();it++)
				{
					if (it->id>num)
					{	
						MeshSet.insert(it,im);
						insertflag++;
						break;
					}
				}
				if(insertflag == 0){
					MeshSet.push_back(im);
				}
			}
		}else if (ch[0] == 'U')
		{
			if (ch[1]=='R') //右侧
			{
				char numc = ch[2];
				int num = 10 - (numc-'0');
				IdMesh im;
				im.id = num;
				im.mm = mm;
				im.haveBottom = false;
				im.deleteFn = 0;
				int insertflag = 0;
				for (it = MeshSet.begin();it != MeshSet.end();it++)
				{
					if (it->id>num)
					{	
						MeshSet.insert(it,im);
						insertflag++;
						break;
					}
				}
				if(insertflag == 0){
					MeshSet.push_back(im);
				}
			}else if (ch[1]=='L') //左侧
			{
				char numc = ch[2];
				int num = 10 + (numc-'0');
				IdMesh im;
				im.id = num;
				im.mm = mm;	
				im.haveBottom = false;
				im.deleteFn = 0;
				int insertflag = 0;
				for (it = MeshSet.begin();it != MeshSet.end();it++)
				{
					if (it->id>num)
					{	
						MeshSet.insert(it,im);
						insertflag++;
						break;
					}
				}
				if(insertflag == 0){
					MeshSet.push_back(im);
				}
			}
		}
	}
	for (int i = 0 ;i< MeshSet.size();i++)
	{
		MeshSet[i].id = i;
	}
}

void EditRepairPlugin::slotImportSetTem()
{
	QString dir = QFileDialog::getOpenFileName(this->gla, tr("Open Set File..."),
		currentfilepath,
		tr("Setting (*.set)"));

	int index = dir.lastIndexOf("/");
	currentfilepath = dir;
	currentfilepath.remove(index,dir.length()-index);

	std::string str = dir.toStdString();
	ifstream inputs;
	std::string oneline;
	inputs.open(str.c_str());
	int i = -1;
	while(getline(inputs, oneline)){
		istringstream strStream(oneline);
		std::string kind,name;
		strStream >> kind;
		double ax,ay,az;
		if (kind == "#")
		{
			i = -1;
			strStream >>name;
			if (name.compare(mesh_tem->shortName().toStdString())==0)
			{
				i = 0;
			}
		}
		if(kind == "m" && i != -1){
			strStream >> ax >> ay >> az;
			barycentreTem.X() = ax;
			barycentreTem.Y() = ay;
			barycentreTem.Z() = az;
		}
		if(kind == "x" && i != -1){
			strStream >> ax >> ay >> az;
			xcoordTem.X() = ax;
			xcoordTem.Y() = ay;
			xcoordTem.Z() = az;
		}
		if(kind == "y" && i != -1){
			strStream >> ax >> ay >> az;
			ycoordTem.X() = ax;
			ycoordTem.Y() = ay;
			ycoordTem.Z() = az;
		}
		if(kind == "z" && i != -1){
			strStream >> ax >> ay >> az;
			zcoordTem.X() = ax;
			zcoordTem.Y() = ay;
			zcoordTem.Z() = az;
		}
	}
	inputs.close();
}

void EditRepairPlugin::slotAlignTem()
{
	CMeshO::VertexIterator vi;
	MatrixXf matrixZ = RotationVector(zcoordTem,zcoord);
	Eigen::Vector3f axisX(xcoordTem.X(),xcoordTem.Y(),xcoordTem.Z());
	Eigen::MatrixXf axisX1(3,1);
	axisX1 = matrixZ * axisX;
	vcg::Point3f axisX2= vcg::Point3f(axisX1(0,0),axisX1(1,0),axisX1(2,0));
	MatrixXf matrixX = RotationVector(axisX2,xcoord);
	for (vi = mesh_tem->cm.vert.begin(); vi != mesh_tem->cm.vert.end(); vi++)
	{
		vcg::Point3f ptemp = vi->P()-barycentreTem;
		Eigen::Vector3f p(ptemp.X(),ptemp.Y(),ptemp.Z());
		Eigen::MatrixXf rp(3,1);
		rp = matrixX * (matrixZ * p);
		ptemp= vcg::Point3f(rp(0,0),rp(1,0),rp(2,0));
		vi->P() = ptemp + barycentre;
	}
	for (vi = mesh_tem_d->cm.vert.begin(); vi != mesh_tem_d->cm.vert.end(); vi++)
	{
		vcg::Point3f ptemp = vi->P()-barycentreTem;
		Eigen::Vector3f p(ptemp.X(),ptemp.Y(),ptemp.Z());
		Eigen::MatrixXf rp(3,1);
		rp = matrixX * (matrixZ * p);
		ptemp= vcg::Point3f(rp(0,0),rp(1,0),rp(2,0));
		vi->P() = ptemp ;
	}

	vcg::tri::UpdateBounding<CMeshO>::Box(mesh_tem->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_tem->cm);
	vcg::tri::UpdateBounding<CMeshO>::Box(mesh_tem_d->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_tem_d->cm);
}

void EditRepairPlugin::slotAdjustX(double delta)
{
	CMeshO::VertexIterator vi;
	CMeshO::VertexIterator vii;

	for (vi = mesh_tem->cm.vert.begin(),vii = mesh_tem_d->cm.vert.begin(); vi != mesh_tem->cm.vert.end(); vi++,vii++)
	{
		xTrans = delta;
		vi->P() = vii->P() * scaleTrans + barycentre + xcoord * xTrans + ycoord * yTrans + zcoord* zTrans;
	}
	//vcg::tri::UpdateBounding<CMeshO>::Box(mesh_tem->cm);
	//vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_tem->cm);
	gla->update();
}

void EditRepairPlugin::slotAdjustY(double delta)
{
	CMeshO::VertexIterator vi;
	CMeshO::VertexIterator vii;

	for (vi = mesh_tem->cm.vert.begin(),vii = mesh_tem_d->cm.vert.begin(); vi != mesh_tem->cm.vert.end(); vi++,vii++)
	{
		yTrans = delta;
		vi->P() = vii->P() * scaleTrans + barycentre + xcoord * xTrans + ycoord * yTrans + zcoord* zTrans;
	}
	//vcg::tri::UpdateBounding<CMeshO>::Box(mesh_tem->cm);
	//vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_tem->cm);
	gla->update();
}

void EditRepairPlugin::slotAdjustZ(double delta)
{
	CMeshO::VertexIterator vi;
	CMeshO::VertexIterator vii;

	for (vi = mesh_tem->cm.vert.begin(),vii = mesh_tem_d->cm.vert.begin(); vi != mesh_tem->cm.vert.end(); vi++,vii++)
	{
		zTrans = delta;
		vi->P() = vii->P() * scaleTrans + barycentre + xcoord * xTrans + ycoord * yTrans + zcoord* zTrans;
	}
	//vcg::tri::UpdateBounding<CMeshO>::Box(mesh_tem->cm);
	//vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_tem->cm);
	gla->update();
}

void EditRepairPlugin::slotAdjustScale(double delta)
{
	CMeshO::VertexIterator vi;
	CMeshO::VertexIterator vii;

	for (vi = mesh_tem->cm.vert.begin(),vii = mesh_tem_d->cm.vert.begin(); vi != mesh_tem->cm.vert.end(); vi++,vii++)
	{
		scaleTrans = delta;
		vi->P() = vii->P() * scaleTrans + barycentre + xcoord * xTrans + ycoord * yTrans + zcoord* zTrans;
	}
	//vcg::tri::UpdateBounding<CMeshO>::Box(mesh_tem->cm);
	//vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_tem->cm);
	gla->update();
}

void EditRepairPlugin::slotSaveAdjust()
{
}

void EditRepairPlugin::slotProjectBarycentre2()
{
	clock_t start,finish;
	double totaltime;
	start=clock();

	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	float t;
	for(vit = mesh_tem->cm.vert.begin(); vit != mesh_tem->cm.vert.end(); ++vit){
		for(fi = mesh->cm.face.begin(); fi != mesh->cm.face.end(); ++fi){
			ProLine pl;
			pl.lp = vit->P();
			pl.ldir = (vit->P()-barycentre)/Distance(vit->P(),barycentre);
			vcg::Point3f fdir = (fi->V(0)->P()-barycentre)/Distance(fi->V(0)->P(),barycentre);
			if (intersectTri2(pl,fi->V(0)->P(),fi->V(1)->P(),fi->V(2)->P(), fdir, 6.0, -6.0, t))
			{
				vit->C() = Color4b::Red;
				break;
			}
		}
	}
	mesh_tem->updateDataMask(MeshModel::MM_VERTCOLOR);

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n此程序的运行时间为"<<totaltime<<endl;
}

void EditRepairPlugin::slotProjectBarycentre3()
{
	clock_t start,finish;
	double totaltime;
	start=clock();

	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	float t;
	float t_xmin,t_xmax,t_ymin,t_ymax,t_zmin,t_zmax;
	float x_min = mesh->cm.bbox.min.X();
	float x_max = mesh->cm.bbox.max.X();
	float y_min = mesh->cm.bbox.min.Y();
	float y_max = mesh->cm.bbox.max.Y();
	float z_min = mesh->cm.bbox.min.Z();
	float z_max = mesh->cm.bbox.max.Z();
	cout<<"x_min:"<<x_min<<endl;
	cout<<"x_max:"<<x_max<<endl;
	cout<<"y_min:"<<y_min<<endl;
	cout<<"y_max:"<<y_max<<endl;
	cout<<"z_min:"<<z_min<<endl;
	cout<<"z_max:"<<z_max<<endl;

	for(vit = mesh_tem->cm.vert.begin(); vit != mesh_tem->cm.vert.end(); ++vit){
		ProLine pl;
		pl.lp = vit->P();
		pl.ldir = (vit->P()-barycentre)/Distance(vit->P(),barycentre);

		float a = 1/pl.ldir.X();
		if (a >= 0)
		{
			t_xmin = a*(x_min - pl.lp.X());
			t_xmax = a*(x_max - pl.lp.X());
		}else{
			t_xmin = a*(x_max - pl.lp.X());
			t_xmax = a*(x_min - pl.lp.X());
		}

		a = 1/pl.ldir.Y();
		if (a >= 0)
		{
			t_ymin = a*(y_min - pl.lp.Y());
			t_ymax = a*(y_max - pl.lp.Y());
		}else{
			t_ymin = a*(y_max - pl.lp.Y());
			t_ymax = a*(y_min - pl.lp.Y());
		}

		a = 1/pl.ldir.Z();
		if (a >= 0)
		{
			t_zmin = a*(z_min - pl.lp.Z());
			t_zmax = a*(z_max - pl.lp.Z());
		}else{
			t_zmin = a*(z_max - pl.lp.Z());
			t_zmax = a*(z_min - pl.lp.Z());
		}

		float d0;
		float d1;

		float d2;
		float d3;
		if (t_xmin > t_ymin)
		{
			d0 = t_xmin;
			if (t_xmax < t_ymax)
			{
				d1 = t_xmax;
			}else{
				d1 = t_ymax;
			}
		}else{
			d0 = t_ymin;
			if (t_xmax < t_ymax)
			{
				d1 = t_xmax;
			}else{
				d1 = t_ymax;
			}
		}
		if (d0 >= d1)
		{
			cout<<"continue1"<<endl;
			continue;
		}else{
			if (t_zmin > d0)
			{
				d2 = t_zmin;
				if (t_zmax < d1)
				{
					d3 = t_zmax;
				}else{
					d3 = d1;
				}
			}else{
				d2 = d0;
				if (t_zmax < d1)
				{
					d3 = t_zmax;
				}else{
					d3 = d1;
				}
			}
			if (d2 >= d3)
			{
				cout<<"continue2"<<endl;
				continue;
			}
		}
		if ((t_xmin>t_ymax)||(t_ymin>t_xmax)||(t_xmin>t_zmax)||(t_zmin>t_xmax)||(t_zmin>t_ymax)||(t_ymin>t_zmax))
		{
			cout<<"continue2"<<endl;
			continue;
		}else{
			for(fi = mesh->cm.face.begin(); fi != mesh->cm.face.end(); ++fi){
				vcg::Point3f fdir = (fi->V(0)->P()-barycentre)/Distance(fi->V(0)->P(),barycentre);
				if (intersectTri2(pl,fi->V(0)->P(),fi->V(1)->P(),fi->V(2)->P(), fdir, 6.0, -6.0, t))
				{
					vit->C() = Color4b::Red;
					break;
				}
			}
		}

	}
	mesh_tem->updateDataMask(MeshModel::MM_VERTCOLOR);

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n此程序的运行时间为"<<totaltime<<endl;
}

void EditRepairPlugin::slotProjectBarycentre4()
{
	clock_t start,finish;
	double totaltime;
	start=clock();

	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	float t;
	float t_xmin,t_xmax,t_ymin,t_ymax,t_zmin,t_zmax;
	float x_min = mesh->cm.bbox.min.X();
	float x_max = mesh->cm.bbox.max.X();
	float y_min = mesh->cm.bbox.min.Y();
	float y_max = mesh->cm.bbox.max.Y();
	float z_min = mesh->cm.bbox.min.Z();
	float z_max = mesh->cm.bbox.max.Z();

	vc1.clear();
	vc2.clear();
	vc3.clear();

	int i = 0;
	int j = 0;
	int k = 0;
	for(vit = mesh_tem->cm.vert.begin(); vit != mesh_tem->cm.vert.end(); ++vit){
		ProLine pl;
		pl.lp = vit->P();
		pl.ldir = (vit->P()-barycentre)/Distance(vit->P(),barycentre);

		vcg::Point3f dirfrac = vcg::Point3f(1.0/pl.ldir.X(),1.0/pl.ldir.Y(),1.0/pl.ldir.Z());

		t_xmin = dirfrac.X()*(x_min - pl.lp.X());
		t_xmax = dirfrac.X()*(x_max - pl.lp.X());
		t_ymin = dirfrac.Y()*(y_min - pl.lp.Y());
		t_ymax = dirfrac.Y()*(y_max - pl.lp.Y());
		t_zmin = dirfrac.Z()*(z_min - pl.lp.Z());
		t_zmax = dirfrac.Z()*(z_max - pl.lp.Z());

		float tmin = ComMax(ComMax(ComMin(t_xmin, t_xmax), ComMin(t_ymin, t_ymax)),ComMin(t_zmin, t_zmax));
		float tmax = ComMin(ComMin(ComMax(t_xmin, t_xmax), ComMax(t_ymin, t_ymax)),ComMax(t_zmin, t_zmax));

		// if tmin > tmax, ray doesn't intersect AABB
		if (tmin > tmax)
		{
			i++;
			vc1.push_back(vit->P());
			cout<<"continue1_"<<i<<endl;
			continue;
		}

		for(fi = mesh->cm.face.begin(); fi != mesh->cm.face.end(); ++fi){
			vcg::Point3f fdir = (fi->V(0)->P()-barycentre)/Distance(fi->V(0)->P(),barycentre);
			if (intersectTri2(pl,fi->V(0)->P(),fi->V(1)->P(),fi->V(2)->P(), fdir, 6.0, -6.0, t))
			{
				j++;
				vit->C() = Color4b::Red;
				break;
			}
		}
		k++;
	}
	cout<<i<<" "<<j <<" "<<k-j<<endl;
	cout<<"x_min:"<<x_min<<endl;
	cout<<"x_max:"<<x_max<<endl;
	cout<<"y_min:"<<y_min<<endl;
	cout<<"y_max:"<<y_max<<endl;
	cout<<"z_min:"<<z_min<<endl;
	cout<<"z_max:"<<z_max<<endl;
	mesh_tem->updateDataMask(MeshModel::MM_VERTCOLOR);

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n此程序的运行时间为"<<totaltime<<endl;
}

void EditRepairPlugin::slotProjectBarycentre()
{
	clock_t start,finish;
	double totaltime;
	start=clock();

	AIndex gIndex;
	gIndex.Set(mesh->cm.face.begin(), mesh->cm.face.end());

	CMeshO::VertexIterator vit;

	vc1.clear();
	vc2.clear();
	vc3.clear();

	for(vit = mesh_tem->cm.vert.begin(); vit != mesh_tem->cm.vert.end(); ++vit){
		ProLine pl;
		pl.lp = vit->P();
		pl.ldir = (vit->P()-barycentre)/Distance(vit->P(),barycentre);

		const bool TEST_BACK_FACES = true;
		vcg::RayTriangleIntersectionFunctor<TEST_BACK_FACES> rayIntersector;
		const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();
		vcg::Ray3<AIndex::ScalarType, false> ray(barycentre, pl.ldir);

		AIndex::ObjPtr isectFace;
		AIndex::ScalarType rayT;
		vcg::EmptyClass a;
		isectFace = gIndex.DoRay(rayIntersector, a , ray, maxDist, rayT);

		if (isectFace != 0) {
			vit->C() = Color4b::Red;
			continue;
		}
	}
	mesh_tem->updateDataMask(MeshModel::MM_VERTCOLOR);

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n此程序的运行时间为"<<totaltime<<endl;
}

void EditRepairPlugin::slotExportMesh()
{
	CMeshO::VertexIterator vit;
	QString filetype = ".txt";
	QString filename = QFileDialog::getSaveFileName(this->gla,tr("Save file"),"*"+filetype);
	if (filename != "")
	{
		QFile f(filename);
		if(!f.open(QIODevice::WriteOnly | QIODevice::Text))
		{
			cout << "Open failed." << endl;
			return;
		}
		QTextStream txtOutput(&f);
		txtOutput <<"#Teeth Points"<<endl;
		//txtOutput <<"#Teeth Gum Line Points outside"<<endl;
		for (vit = mesh_tem->cm.vert.begin(); vit != mesh_tem->cm.vert.end(); ++vit)
		{
			txtOutput <<"{"<<vit->P().X()<<","<<vit->P().Y()<<","<<vit->P().Z()<<"},";
		}
		f.close();
	}
}

void EditRepairPlugin::slotDeformationLaplacian(){
	int i,j,k;
	CMeshO::VertexIterator vi;
	int degree;//顶点度数
	int v_weight = 1000;//权重，越大点越固定
	int VN;//处理顶点的数量
	int flag;
	CMeshO::PerVertexAttributeHandle<int> col = tri::Allocator<CMeshO>::AddPerVertexAttribute<int>(mesh_tem_d->cm ,std::string("v_col"));//给每个顶点增加col属性，以便快速找到该顶点在矩阵中的列号
	for(vi=mesh_tem_d->cm.vert.begin(),i=0; vi!=mesh_tem_d->cm.vert.end(); ++vi,++i)
	{
		if (!(*vi).IsD())
		{				
			col[vi] = -1;//初始化列号属性
		}
	}
	VN=mesh_tem->cm.vert.size();
	j = 0;
	for (i = 0;i < crownBoundaryT.size();i++,j++)
	{
		flag = crownBoundaryT[i];
		vi = mesh_tem_d->cm.vert.begin();
		vi+=flag;
		col[vi] = j;
	}
	for (i = 0;i < rootT.size();i++,j++)
	{
		flag = rootT[i];
		vi = mesh_tem_d->cm.vert.begin();
		vi+=flag;
		col[vi] = j;
	}
	for (i = 0;i < crownInnerT.size();i++,j++)
	{
		flag = crownInnerT[i];
		vi = mesh_tem_d->cm.vert.begin();
		vi+=flag;
		col[vi] = j;
	}

	vector<Triplet<float>> trips;
	j = 0;

	mesh_tem_d->updateDataMask(MeshModel::MM_VERTFACETOPO);
	mesh_tem_d->updateDataMask(MeshModel::MM_FACEFACETOPO);

	for(i=0;i<crownBoundaryT.size();++i,++j)//构造矩阵
	{
		vi = mesh_tem_d->cm.vert.begin();
		flag = crownBoundaryT[i];
		vi+=flag;
		if(!(*vi).IsD()){
			CVertexO *vii = &(*vi);
			vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
			CVertexO* firstV = pos.VFlip();
			CVertexO* tempV= NULL;//访问V的一环邻域
			degree = 0;
			do 
			{
				pos.NextE();
				tempV = pos.VFlip();
				if(col[tempV]!=-1){
					k = col[tempV];
					trips.push_back(Triplet<float>(j, k,-1));
					degree++;
				}
			} while (firstV != tempV);
			trips.push_back(Triplet<float>(j, j,degree));
		}
	}
	for(i=0;i<rootT.size();++i,++j)//构造矩阵
	{
		vi = mesh_tem_d->cm.vert.begin();
		flag = rootT[i];
		vi+=flag;
		if(!(*vi).IsD()){
			CVertexO *vii = &(*vi);
			vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
			CVertexO* firstV = pos.VFlip();
			CVertexO* tempV= NULL;//访问V的一环邻域
			degree = 0;
			do 
			{
				pos.NextE();
				tempV = pos.VFlip();
				if(col[tempV]!=-1){
					k = col[tempV];
					trips.push_back(Triplet<float>(j, k,-1));
					degree++;
				}
			} while (firstV != tempV);
			trips.push_back(Triplet<float>(j, j,degree));
		}
	}
	for(i=0;i<crownInnerT.size();++i,++j)//构造矩阵
	{
		vi = mesh_tem_d->cm.vert.begin();
		flag = crownInnerT[i];
		vi+=flag;
		if(!(*vi).IsD()){
			CVertexO *vii = &(*vi);
			vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
			CVertexO* firstV = pos.VFlip();
			CVertexO* tempV= NULL;//访问V的一环邻域
			degree = 0;
			do 
			{
				pos.NextE();
				tempV = pos.VFlip();
				if(col[tempV]!=-1){
					k = col[tempV];
					trips.push_back(Triplet<float>(j, k,-1));
					degree++;
				}
			} while (firstV != tempV);
			trips.push_back(Triplet<float>(j, j,degree));
		}
	}
	SparseMatrix<float> A0;//Laplacian方阵
	A0.resize(VN,VN);
	A0.setFromTriplets(trips.begin(), trips.end());

	//构造初始坐标向量
	VectorXf X0;
	X0.resize(VN);
	VectorXf Y0;
	Y0.resize(VN);
	VectorXf Z0;
	Z0.resize(VN);
	j = 0;
	for(i =0;i<crownBoundaryT.size();i++,j++)
	{
		X0(j)=mesh_tem->cm.vert[crownBoundaryT[i]].P().X();
		Y0(j)=mesh_tem->cm.vert[crownBoundaryT[i]].P().Y();
		Z0(j)=mesh_tem->cm.vert[crownBoundaryT[i]].P().Z();
	}
	for(i =0;i<rootT.size();i++,j++)
	{
		X0(j)=mesh_tem->cm.vert[rootT[i]].P().X();
		Y0(j)=mesh_tem->cm.vert[rootT[i]].P().Y();
		Z0(j)=mesh_tem->cm.vert[rootT[i]].P().Z();
	}
	for(i =0;i<crownInnerT.size();i++,j++)
	{
		X0(j)=mesh_tem->cm.vert[crownInnerT[i]].P().X();
		Y0(j)=mesh_tem->cm.vert[crownInnerT[i]].P().Y();
		Z0(j)=mesh_tem->cm.vert[crownInnerT[i]].P().Z();
	}

	VectorXf BX0;
	BX0.resize(VN);
	BX0 = A0*X0;
	VectorXf BY0;
	BY0.resize(VN);
	BY0 = A0*Y0;
	VectorXf BZ0;
	BZ0.resize(VN);
	BZ0 = A0*Z0;

	int VN1=VN+crownBoundaryT.size();
	VectorXf BX1;
	BX1.resize(VN1);
	VectorXf BY1;
	BY1.resize(VN1);
	VectorXf BZ1;
	BZ1.resize(VN1);
	for (i = 0;i<VN;i++)
	{
		BX1(i)=BX0(i);
		BY1(i)=BY0(i);
		BZ1(i)=BZ0(i);
	}
	int temp = VN;

	vcg::Matrix33f re3;

	vi=mesh_tem_d->cm.vert.begin();
	for(i = 0;i < crownBoundaryT.size();i++,temp++)
	{
		flag = crownBoundaryT[i];		
		trips.push_back(Triplet<float>(temp, col[vi+flag], v_weight));
		Point3f p3 =( mesh_tem->cm.vert[crownBoundaryT[i]].P()+vcg::Point3f(0.0,0.0,2.0)) * v_weight;
		BX1(temp) = p3.X();
		BY1(temp) = p3.Y();
		BZ1(temp) = p3.Z();
	}
	// 	for(i = 0;i < crownInnerT.size();i++,temp++)
	// 	{
	// 		flag = crownInnerT[i];		
	// 		trips.push_back(Triplet<float>(temp, col[vi+flag], v_weight));
	// 		Point3f p3 =  mesh_tem_d->cm.vert[crownInnerT[i]].P() * v_weight;
	// 		BX1(temp) = p3.X();
	// 		BY1(temp) = p3.Y();
	// 		BZ1(temp) = p3.Z();
	// 	}

	SparseMatrix<float> A;//构造Laplacian增广矩阵
	A.resize(VN1,VN);
	A.setFromTriplets(trips.begin(), trips.end());
	SparseMatrix<float> AT = A.transpose();
	SparseMatrix<float> M = AT*A;

	VectorXf BX2 = AT*BX1;
	VectorXf BY2 = AT*BY1;
	VectorXf BZ2 = AT*BZ1;

	SimplicialLLT<SparseMatrix<float>, Lower> solverA;
	solverA.compute(M);
	VectorXf X;
	VectorXf Y;
	VectorXf Z;
	X = solverA.solve(BX2);
	Y = solverA.solve(BY2);
	Z = solverA.solve(BZ2);

	j = crownBoundaryT.size();
	// 	for (i = 0;i<crownBoundaryT.size();i++,j++)
	// 	{
	// 		mesh_tem_d->cm.vert[crownBoundaryT[i]] = 
	// 		CMN->handleVex[i].cv.X() = X[j];
	// 		CMN->handleVex[i].cv.Y() = Y[j];
	// 		CMN->handleVex[i].cv.Z() = Z[j];
	// 	}
	for (i = 0;i<rootT.size();i++,j++)
	{
		mesh_tem_d->cm.vert[rootT[i]].P().X() = X[j];
		mesh_tem_d->cm.vert[rootT[i]].P().X() = X[j];
		mesh_tem_d->cm.vert[rootT[i]].P().X() = X[j];
	}
	// 	for (i = 0;i<CMN->originstaticVex.size();i++,j++)
	// 	{
	// 		CMN->staticVex[i].cv.X() = X[j];
	// 		CMN->staticVex[i].cv.Y() = Y[j];
	// 		CMN->staticVex[i].cv.Z() = Z[j];
	// 	}
	tri::Allocator<CMeshO>::DeletePerVertexAttribute(mesh_tem_d->cm ,std::string("v_col"));
}

void EditRepairPlugin::slotGetBoundaryTem(){
	CMeshO::VertexIterator vit;
	crownBoundaryT.clear();
	crownInnerT.clear();
	rootT.clear();

	mesh_tem->updateDataMask(MeshModel::MM_VERTFACETOPO);
	mesh_tem->updateDataMask(MeshModel::MM_FACEFACETOPO);

	int i=0;
	for(vit = mesh_tem->cm.vert.begin(); vit != mesh_tem->cm.vert.end(); ++vit){
		bool findCrownBoundaryT = false;
		vcg::Color4f c1 = Color4f::Construct(vit->C());
		if ( c1 == Color4<float>(Color4<float>::Red))//if white color?
		{
			if(!(*vit).IsD()){//判断该点是否存在
				CVertexO *vii = &(*vit);
				vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
				CVertexO* firstV = pos.VFlip();//取和vitd相邻的第一个点
				CVertexO* tempV= NULL;
				do
				{
					pos.NextE();
					tempV = pos.VFlip();
					if (Color4f::Construct(tempV->C())!= Color4<float>(Color4<float>::Red))
					{
						findCrownBoundaryT = true;
						tempV->SetS();
						crownBoundaryT.push_back(i);
						break;
					}
				}while (firstV != tempV);//从第一个点开始按顺序遍历，知道回到第一个点为止
				if (!findCrownBoundaryT)
				{
					crownInnerT.push_back(i);
				}
			}
		}else{
			rootT.push_back(i);
		}
		i++;
	}

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	rootBoundaryT.clear();
	i=0;
	for(vit = mesh_tem->cm.vert.begin(); vit != mesh_tem->cm.vert.end(); ++vit){
		if(!(*vit).IsD()){//判断该点是否存在
			if ((*vit).IsS() && Color4f::Construct(vit->C())!=Color4<float>(Color4<float>::Red))
			{
				rootBoundaryT.push_back(i);
			}
		}
		i++;
	}
	tri::UpdateSelection<CMeshO>::Clear(mesh_tem->cm);
}

void EditRepairPlugin::slotSmoothRoot(){
	CMeshO::VertexIterator vitd=mesh_tem->cm.vert.begin();
	int stepSmoothNum;

	for(vitd = mesh_tem->cm.vert.begin()+mesh->cm.vert.size(); vitd != mesh_tem->cm.vert.end(); ++vitd){
		vcg::Color4f c1 = Color4f::Construct(vitd->C());
		if ( c1 != Color4<float>(Color4<float>::Red))//if white color?
		{
			if(!(*vitd).IsD()){//判断该点是否存在
				(*vitd).SetS();
			}
		}
	}
	stepSmoothNum = 3;

	tri::UpdateFlags<CMeshO>::FaceBorderFromNone(mesh_tem->cm);
	tri::Smooth<CMeshO>::VertexCoordLaplacian(mesh_tem->cm,stepSmoothNum,true,false);
	tri::UpdateBounding<CMeshO>::Box(mesh_tem->cm);
	if(mesh_tem->cm.fn>0) {
		tri::UpdateNormals<CMeshO>::PerFaceNormalized(mesh_tem->cm);
		tri::UpdateNormals<CMeshO>::PerVertexAngleWeighted(mesh_tem->cm);
	}

	tri::UpdateSelection<CMeshO>::ClearVertex(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::ClearFace(mesh_tem->cm);

}

void EditRepairPlugin::slotSmoothGap(){
	int stepSmoothNum;

	CVertexO* vitd;

	for (int i = 0; i < borders.size(); i++)
	{
		vitd = &repaired_mesh->cm.vert[borders[i]];
		if(!vitd->IsD()){//判断该点是否存在
			vitd->SetS();
		}
	}

	int v1 = mesh->cm.vert.size();
	for (int i = 0; i < borders_root.size(); i++)
	{
		vitd = &repaired_mesh->cm.vert[borders_root[i]+v1];
		if(!vitd->IsD()){//判断该点是否存在
			vitd->SetS();
		}
	}

	stepSmoothNum = 1;

	tri::UpdateFlags<CMeshO>::FaceBorderFromNone(repaired_mesh->cm);
	tri::Smooth<CMeshO>::VertexCoordLaplacian(repaired_mesh->cm,stepSmoothNum,true,false);
	tri::UpdateBounding<CMeshO>::Box(repaired_mesh->cm);
	if(repaired_mesh->cm.fn>0) {
		tri::UpdateNormals<CMeshO>::PerFaceNormalized(repaired_mesh->cm);
		tri::UpdateNormals<CMeshO>::PerVertexAngleWeighted(repaired_mesh->cm);
	}

	tri::UpdateSelection<CMeshO>::ClearVertex(repaired_mesh->cm);
	tri::UpdateSelection<CMeshO>::ClearFace(repaired_mesh->cm);

	gla->update();
}

void EditRepairPlugin::slotSmoothGap2(){
	int stepSmoothNum;

	CVertexO* vitd;

	for (int i = 0; i < borders.size(); i++)
	{
		vitd = &mesh_tem->cm.vert[borders[i]];
		if(!vitd->IsD()){//判断该点是否存在
			vitd->SetS();
		}
	}

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	// 	int v1 = mesh->cm.vert.size();
	// 	for (int i = 0; i < borders_root.size(); i++)
	// 	{
	// 		vitd = &mesh_tem->cm.vert[borders_root[i]+v1];
	// 		if(!vitd->IsD()){//判断该点是否存在
	// 			vitd->SetS();
	// 		}
	// 	}

	stepSmoothNum = 3;

	tri::UpdateFlags<CMeshO>::FaceBorderFromNone(mesh_tem->cm);
	tri::Smooth<CMeshO>::VertexCoordLaplacian(mesh_tem->cm,stepSmoothNum,true,false);
	tri::UpdateBounding<CMeshO>::Box(mesh_tem->cm);
	if(mesh_tem->cm.fn>0) {
		tri::UpdateNormals<CMeshO>::PerFaceNormalized(mesh_tem->cm);
		tri::UpdateNormals<CMeshO>::PerVertexAngleWeighted(mesh_tem->cm);
	}

	tri::UpdateSelection<CMeshO>::ClearVertex(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::ClearFace(mesh_tem->cm);

	gla->update();
}

void EditRepairPlugin::slotDefomationTem(){
	clock_t start,finish;
	double totaltime;
	start=clock();

	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;

	AIndex gIndex;
	gIndex.Set(mesh->cm.face.begin(), mesh->cm.face.end());

	for(vi = mesh_tem->cm.vert.begin(); vi != mesh_tem->cm.vert.end(); ++vi){
		if (!vi->IsD())
		{
			ProLine pl;
			pl.lp = vi->P();
			pl.ldir = (vi->P()-barycentre)/Distance(vi->P(),barycentre);

			const bool TEST_BACK_FACES = true;
			vcg::RayTriangleIntersectionFunctor<TEST_BACK_FACES> rayIntersector;
			const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();
			vcg::Ray3<AIndex::ScalarType, false> ray(barycentre, pl.ldir);

			AIndex::ObjPtr isectFace;
			AIndex::ScalarType rayT;
			vcg::EmptyClass a;
			isectFace = gIndex.DoRay(rayIntersector, a , ray, maxDist, rayT);

			if (isectFace != 0) {
				vi->P() = barycentre + pl.ldir*rayT;
				vi->C() = Color4b::Red;
				continue;
			}
		}
	}
	mesh_tem->updateDataMask(MeshModel::MM_VERTCOLOR);

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n投影形变的时间为"<<totaltime<<endl;


	//*************划分模板上的区域**************
	CMeshO::VertexIterator vit;
	crownBoundaryT.clear();
	crownInnerT.clear();
	rootT.clear();
	rootBoundaryT.clear();

	int i=0;
	for(vit = mesh_tem->cm.vert.begin(); vit != mesh_tem->cm.vert.end(); ++vit){
		bool findCrownBoundaryT = false;
		vcg::Color4f c1 = Color4f::Construct(vit->C());
		if ( c1 == Color4<float>(Color4<float>::Red))//if white color?
		{
			if(!(*vit).IsD()){//判断该点是否存在
				CVertexO *vii = &(*vit);
				vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
				CVertexO* firstV = pos.VFlip();//取和vitd相邻的第一个点
				CVertexO* tempV= NULL;
				do
				{
					pos.NextE();
					tempV = pos.VFlip();
					if (Color4f::Construct(tempV->C())!= Color4<float>(Color4<float>::Red))
					{
						findCrownBoundaryT = true;
						tempV->SetS();
						crownBoundaryT.push_back(i);
						break;
					}
				}while (firstV != tempV);//从第一个点开始按顺序遍历，知道回到第一个点为止
				if (!findCrownBoundaryT)
				{
					crownInnerT.push_back(i);
				}
			}
		}else{
			rootT.push_back(i);
		}
		i++;
	}

	//六环领域
	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);

	tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(mesh_tem->cm);


	i=0;
	for(vit = mesh_tem->cm.vert.begin(); vit != mesh_tem->cm.vert.end(); ++vit){
		if(!(*vit).IsD()){//判断该点是否存在
			if ((*vit).IsS() && Color4f::Construct(vit->C())!=Color4<float>(Color4<float>::Red))
			{
				rootBoundaryT.push_back(i);
			}
		}
		i++;
	}
	tri::UpdateSelection<CMeshO>::Clear(mesh_tem->cm);

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n区域划分的时间为"<<totaltime<<endl;

	//*****************非投影区域形变**********************
	//CMeshO::VertexIterator vit;
	int j = 0;
	DIF = 1000;
	while (DIF > 0.1)
	{
		DIF = 0.0;
		CVertexO* vitd;
		std::vector <vcg::Point3f> crownBoundary;
		for (int i = 0; i < rootBoundaryT.size(); i++)
		{
			vcg::Point3f tran = vcg::Point3f(0.0,0.0,0.0);
			int n = 0;
			vitd = &mesh_tem->cm.vert[rootBoundaryT[i]];
			if(!(*vitd).IsD()){//判断该点是否存在
				CVertexO *vii = &(*vitd);
				vcg::face::JumpingPos<CFaceO> pos(vii->VFp(),vii);
				CVertexO* firstV = pos.VFlip();//取和vitd相邻的第一个点
				CVertexO* tempV= NULL;
				crownBoundary.clear();
				do
				{
					pos.NextE();
					tempV = pos.VFlip();
					if (Color4f::Construct(tempV->C())== Color4<float>(Color4<float>::Red))
					{
						crownBoundary.push_back(tempV->P());
					}
					int num = tempV - &*(mesh_tem->cm.vert.begin());
					tran += (tempV->P() - (mesh_tem_d->cm.vert.begin() + num)->P());
					n++;
				}while (firstV != tempV);//从第一个点开始按顺序遍历，知道回到第一个点为止
				int num2 = &*vitd - &*(mesh_tem->cm.vert.begin());
				vcg::Point3f temp = (mesh_tem_d->cm.vert.begin() + num2)->P() + tran/n;
				if (crownBoundary.size()>0)
				{
					for ( int i = 0; i < crownBoundary.size();i++)
					{
						temp = temp + crownBoundary[i];
					}
					temp = temp / (crownBoundary.size()+1);
				}

				DIF += Distance(temp,vitd->P());
				vitd->P() = temp;
			}
		}
		cout<<"slotDeformationNonPro "<<j<<" :"<<DIF<<endl;
		j++;
	}

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n非投影区域的形变时间为"<<totaltime<<endl;

	gla->update();
}

void EditRepairPlugin::slotRepairByTem(){
	clock_t start,finish;
	double totaltime;
	start=clock();


	//mesh_tem_d = mesh_tem;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	for(fi = mesh_tem->cm.face.begin(); fi != mesh_tem->cm.face.end(); ++fi){
		if(!(*fi).IsD()){//判断该点是否存在
			for (int i = 0;i<3;i++)
			{
				if(Color4f::Construct((fi->V(i))->C())==Color4<float>(Color4<float>::Red)){
					fi->SetS();
				} 
			}
		}
	}
	tri::UpdateSelection<CMeshO>::ClearVertex(mesh_tem->cm);
	tri::UpdateSelection<CMeshO>::VertexFromFaceStrict(mesh_tem->cm);
	for(fi=mesh_tem->cm.face.begin();fi!=mesh_tem->cm.face.end();++fi)
		if(!(*fi).IsD() && (*fi).IsS() )
			tri::Allocator<CMeshO>::DeleteFace(mesh_tem->cm,*fi);
	for(vi=mesh_tem->cm.vert.begin();vi!=mesh_tem->cm.vert.end();++vi)
		if(!(*vi).IsD() && (*vi).IsS() )
			tri::Allocator<CMeshO>::DeleteVertex(mesh_tem->cm,*vi);
	mesh_tem->clearDataMask(MeshModel::MM_FACEFACETOPO);
	mesh_tem->clearDataMask(MeshModel::MM_FACEFLAGBORDER);

	tri::UpdateSelection<CMeshO>::Clear(mesh_tem->cm);

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n删除投影区域的时间为"<<totaltime<<endl;

	//新建替换网格复制底部区域
	// 	CMeshO::VertexIterator vi;
	// 	CMeshO::FaceIterator fi;
	std::vector<CVertexO> vv_mm;
	std::vector<CFaceO> vf_mm;
	std::vector<int> VertexId(mesh_tem->cm.vert.size());
	int numvert = 0;

	for (vi=mesh_tem->cm.vert.begin();vi!=mesh_tem->cm.vert.end();vi++){
		if (!(*vi).IsD())
		{
			VertexId[vi-mesh_tem->cm.vert.begin()]=numvert;
			vv_mm.push_back(*vi);
			numvert++;
		}
	}
	for(fi=mesh_tem->cm.face.begin(); fi!=mesh_tem->cm.face.end(); ++fi)
	{
		if( !(*fi).IsD() ){
			vf_mm.push_back(*fi);
		}
	}

	MeshModel *newm;
	newm =  md->addNewMesh(mesh_tem->fullName(),"",false);

	newm->cm.Clear();
	int VN,FN;
	VN=vv_mm.size();
	FN=vf_mm.size();
	vcg::tri::Allocator<CMeshO>::AddVertices(newm->cm,VN);
	vcg::tri::Allocator<CMeshO>::AddFaces(newm->cm,FN);
	CMeshO::VertexPointer *ivp = new CMeshO::VertexPointer[VN];
	int cnt=0;
	vi= newm->cm.vert.begin();
	std::vector<CVertexO>::const_iterator itv;
	for(itv=vv_mm.begin();itv!=vv_mm.end();itv++){
		ivp[cnt]=&*vi;
		(*vi).P()= (*itv).P();
		++vi;
		++cnt;
	}
	fi=newm->cm.face.begin();
	std::vector<CFaceO>::const_iterator itf;
	for (itf = vf_mm.begin(); itf!= vf_mm.end(); itf++)
	{
		for(int k=0;k<(*itf).VN();k++)
		{
			int vInd = VertexId[(*itf).V(k)-&*(mesh_tem->cm.vert.begin())];
			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}

	vcg::tri::UpdateBounding<CMeshO>::Box(newm->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(newm->cm);
	md->delMesh(mesh_tem);

	newm->clearDataMask(MeshModel::MM_FACEMARK);
	newm->updateDataMask(MeshModel::MM_FACEMARK);
	newm->updateDataMask(MeshModel::MM_FACEFACETOPO);

	mesh_tem = newm;

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n新建替换网格的时间为"<<totaltime<<endl;



	mesh->updateDataMask(MeshModel::MM_FACEFLAGBORDER);
	borders.clear();
	//CMeshO::VertexIterator vi;
	int i = 0;
	for(vi = mesh->cm.vert.begin(); vi != mesh->cm.vert.end(); ++vi){
		if( !(*vi).IsD() )
		{
			if((*vi).IsB() )
			{
				borders.push_back(i);
			}
		}
		i++;
	}
	//mesh_tem_d = mesh_tem;//test
	mesh_tem->updateDataMask(MeshModel::MM_FACEFLAGBORDER);
	borders_root.clear();
	i = 0;
	for(vi = mesh_tem->cm.vert.begin(); vi != mesh_tem->cm.vert.end(); ++vi){
		if( !(*vi).IsD() )
		{
			if((*vi).IsB() )
			{
				borders_root.push_back(i);
			}
		}
		i++;
	}

	slotSortBorder2();

	slotFillGap2();
	cout<<"slotFillGap2"<<endl;

	slotSmoothGap2();
	cout<<"slotSmoothGap2"<<endl;

	slotSmoothRoot();

	gla->update();
}

bool EditRepairPlugin::StartEdit(MeshDocument &_md, GLArea *_gla )
{
	this->md = &_md;
	if (md->mm() == NULL)
		return false;//没有模型无法启动

	this->gla = _gla;
	getSE = false;//Qiainjiahong Test
	importset = false;
	mesh = md->mm();
	barycentre = vcg::Point3f(0.0, 0.0, 0.0); 
	xcoord = vcg::Point3f(1.0, 0.0, 0.0);
	ycoord = vcg::Point3f(0.0, 1.0, 0.0);
	zcoord = vcg::Point3f(0.0, 0.0, 1.0);

	// 	mesh = md->meshList[0];
	// 	mesh2 = md->meshList[1];

// 	QString qst1 = "LL1.obj";
// 	QString qst2 = "LL1_htr.obj";
// 	for(int i = 0; i < md->meshList.size(); i++){
// 		QString qstr = md->meshList[i]->shortName();
// 		char*  ch;
// 		QByteArray ba = qstr.toLatin1();
// 		ch=ba.data();
// 		if (ch[3] == 'r')
// 		{
// 			mesh2 = md->meshList[i];
// 		}else if (ch[3] == '.')
// 		{
// 			mesh = md->meshList[i];
// 		}
// 	}

	gla->setCursor(QCursor(QPixmap(":/images/cur_info.png"),1,1));
	if(stagesDialog==0){
		stagesDialog = new RepairDialog(gla->window(),this);
		connect(stagesDialog->ui.projectButton,SIGNAL(clicked()),this,SLOT(slotProject()));
		connect(stagesDialog->ui.deformationButton,SIGNAL(clicked()),this,SLOT(slotDeformationPro()));
		connect(stagesDialog->ui.deformatinoNonProButton,SIGNAL(clicked()),this,SLOT(slotDeformationNonPro()));
		connect(stagesDialog->ui.deformationAutoButton,SIGNAL(clicked()),this,SLOT(slotDeformationAuto()));
		connect(stagesDialog->ui.deleteProButton,SIGNAL(clicked()),this,SLOT(slotDeletePro()));
		connect(stagesDialog->ui.importSetButton,SIGNAL(clicked()),this,SLOT(slotImportSet()));
		connect(stagesDialog->ui.fillGapButton,SIGNAL(clicked()),this,SLOT(slotFillGap()));
		connect(stagesDialog->ui.findBoundaryButton,SIGNAL(clicked()),this,SLOT(slotFindBoundary()));
		connect(stagesDialog->ui.sortBorderButton,SIGNAL(clicked()),this,SLOT(slotSortBorder()));
		connect(stagesDialog->ui.copyMeshButton,SIGNAL(clicked()),this,SLOT(slotCopyMesh()));
		connect(stagesDialog->ui.importSetTemButton,SIGNAL(clicked()),this,SLOT(slotImportSetTem()));
		connect(stagesDialog->ui.alignTemButton,SIGNAL(clicked()),this,SLOT(slotAlignTem()));
		connect(stagesDialog->ui.xAdjustdoubleSpinBox,SIGNAL(valueChanged(double)), this,SLOT(slotAdjustX(double)));
		connect(stagesDialog->ui.yAdjustdoubleSpinBox,SIGNAL(valueChanged(double)), this,SLOT(slotAdjustY(double)));
		connect(stagesDialog->ui.zAdjustdoubleSpinBox,SIGNAL(valueChanged(double)), this,SLOT(slotAdjustZ(double)));
		connect(stagesDialog->ui.scaledoubleSpinBox,SIGNAL(valueChanged(double)),this,SLOT(slotAdjustScale(double)));
		connect(stagesDialog->ui.saveAdjustButton,SIGNAL(clicked()),this,SLOT(slotSaveAdjust()));
		connect(stagesDialog->ui.projectBarycentreButton,SIGNAL(clicked()),this,SLOT(slotProjectBarycentre()));
		connect(stagesDialog->ui.deformationButton2,SIGNAL(clicked()),this,SLOT(slotDeformationPro2()));
		connect(stagesDialog->ui.deformationAutoButton2,SIGNAL(clicked()),this,SLOT(slotDeformationAuto2()));
		connect(stagesDialog->ui.exportMeshButton,SIGNAL(clicked()),this,SLOT(slotExportMesh()));
		connect(stagesDialog->ui.getBoundaryTemButton,SIGNAL(clicked()),this,SLOT(slotGetBoundaryTem()));
		connect(stagesDialog->ui.smoothRootButton,SIGNAL(clicked()),this,SLOT(slotSmoothRoot()));
		connect(stagesDialog->ui.smoothGapButton,SIGNAL(clicked()),this,SLOT(slotSmoothGap()));
		connect(stagesDialog->ui.defomationTemButton,SIGNAL(clicked()),this,SLOT(slotDefomationTem()));
		connect(stagesDialog->ui.repairButton,SIGNAL(clicked()),this,SLOT(slotRepairByTem()));

		//New
		connect(stagesDialog->ui.importSetButton_2,SIGNAL(clicked()),this,SLOT(slotImportSet()));
		connect(stagesDialog->ui.cutButton,SIGNAL(clicked()),this,SLOT(slotCut()));

		//Qianjiahong
		connect(stagesDialog->ui.mergeButton,SIGNAL(clicked()),this,SLOT(slotMerge()));
		connect(stagesDialog->ui.normalizeButton,SIGNAL(clicked()),this,SLOT(slotNormalize()));
		connect(stagesDialog->ui.exportPButton,SIGNAL(clicked()),this,SLOT(slotExportPoints()));
		connect(stagesDialog->ui.exportPButton2,SIGNAL(clicked()),this,SLOT(slotExportPoints2()));
		connect(stagesDialog->ui.errorMapButton,SIGNAL(clicked()),this,SLOT(slotErrorMap()));

		connect(stagesDialog->ui.invertNormalButton,SIGNAL(clicked()),this,SLOT(slotInvertNormal()));
		connect(stagesDialog->ui.invertNormalFaceButton,SIGNAL(clicked()),this,SLOT(slotInvertNormalFace()));

		connect(stagesDialog->ui.findNeastestPointsButton,SIGNAL(clicked()),this,SLOT(slotFindNeastestPoints()));
		connect(stagesDialog->ui.extendBorderButton,SIGNAL(clicked()),this,SLOT(slotExtendBorder()));
		connect(stagesDialog->ui.linkPointsButton,SIGNAL(clicked()),this,SLOT(slotLinkPoints()));
		connect(stagesDialog->ui.importInsectionPButton,SIGNAL(clicked()),this,SLOT(slotImportInsectionP()));
	}

	stagesDialog->show();
	SetMode(TranMode::NONE);
	gla->update();
	return true;
}

bool EditRepairPlugin::IntersectionPlaneMesh(CMeshO & m, Plane3f pl,CMeshO & em)
{
	std::vector<Point3f> ptVec;
	std::vector<Point3f> nmVec;

	CMeshO::PerVertexAttributeHandle <float> qH =tri::Allocator<CMeshO>::AddPerVertexAttribute <float>(m,"TemporaryPlaneDistance");

	CMeshO::VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi) if(!(*vi).IsD())
		qH[vi] = vcg::Distance(pl,(*vi).cP());//这个距离是有方向的，+在上，-在下

	vector<int> ids;
	for(size_t i=0;i<m.face.size();i++){

		if(!m.face[i].IsD())
		{
			ptVec.clear();
			nmVec.clear();
			for(int j=0;j<3;++j)
			{
				if((qH[m.face[i].V0(j)] * qH[m.face[i].V1(j)])<0)
				{
					const Point3<float> &p0 = m.face[i].V0(j)->cP();
					const Point3<float> &p1 = m.face[i].V1(j)->cP();
					const Point3<float> &n0 = m.face[i].V0(j)->cN();
					const Point3<float> &n1 = m.face[i].V1(j)->cN();
					float q0 = qH[m.face[i].V0(j)];
					float q1 = qH[m.face[i].V1(j)];
					//         printf("Intersection ( %3.2f %3.2f %3.2f )-( %3.2f %3.2f %3.2f )\n",p0[0],p0[1],p0[2],p1[0],p1[1],p1[2]);
					Point3<float> pp;
					Segment3<float> seg(p0,p1);
					IntersectionPlaneSegment(pl,seg,pp);
					ptVec.push_back(pp);
					Point3<float> nn =(n0*fabs(q1) + n1*fabs(q0))/fabs(q0-q1);
					nmVec.push_back(nn);
				}
				if(qH[m.face[i].V(j)]==0)  ptVec.push_back(m.face[i].V(j)->cP());
			}
			if(ptVec.size()>=2)
			{
				CMeshO::VertexIterator vi;
				vcg::tri::Allocator<CMeshO>::AddEdges(em,1);
				vi = vcg::tri::Allocator<CMeshO>::AddVertices(em,2);
				(*vi).P() = ptVec[0];
				/*		(*vi).P().Y() = ptVec[0].Y();
				(*vi).P().Z() = ptVec[0].Z();*/
				(*vi).N() = nmVec[0];
				ids.push_back(i);
				em.edge.back().V(0) = &(*vi);
				vi++;
				(*vi).P() = ptVec[1];
				//(*vi).P().Y() = ptVec[1].Y();
				//(*vi).P().Z() = ptVec[1].Z()*0.01+i;
				(*vi).N() = nmVec[1];
				ids.push_back(i);
				em.edge.back().V(1) = &(*vi);
			}
		}
	}
	faceids.push_back(ids);
	tri::Allocator<CMeshO>::DeletePerVertexAttribute < float >(m,qH);

	return true;
}

void CapEdgeMesh(CMeshO &em, CMeshO &cm, bool revertFlag=false)
{
	// 	typedef CMeshO::edge EdgeType;
	// 	std::vector< std::vector<Point3f> > outlines;
	// 	std::vector<Point3f> outline;
	// 	//UpdateFlags<CMeshO>::EdgeClearV(em);
	// 	//UpdateTopology<CMeshO>::EdgeEdge(em);
	// 	int nv=0;
	// 	for(size_t i=0;i<em.edge.size();i++) if(!em.edge[i].IsD())
	// 	{
	// 		if (!em.edge[i].IsV())
	// 		{
	// 			//edge::Pos<EdgeType> startE(&em.edge[i],0);
	// 			//edge::Pos<EdgeType> curE=startE;
	// 			do
	// 			{
	// 				curE.E()->SetV();
	// 				outline.push_back(curE.V()->P());
	// 				curE.NextE();
	// 				nv++;
	// 			}
	// 			while(curE != startE);
	// 			if(revertFlag) std::reverse(outline.begin(),outline.end());
	// 			outlines.push_back(outline);
	// 			outline.clear();
	// 		}
	// 	}
	// 	if (nv<2) return;
	// 	//  printf("Found %i outlines for a total of %i vertices",outlines.size(),nv);
	// 
	// 	CMeshO::VertexIterator vi=vcg::tri::Allocator<CMeshO>::AddVertices(cm,nv);
	// 	for (size_t i=0;i<outlines.size();i++)
	// 	{
	// 		for(size_t j=0;j<outlines[i].size();++j,++vi)
	// 			(&*vi)->P()=outlines[i][j];
	// 	}
	// 
	// 	std::vector<int> indices;
	// 	glu_tesselator::tesselate(outlines, indices);
	// 	std::vector<Point3f> points;
	// 	glu_tesselator::unroll(outlines, points);
	// 	//typename MeshType::FaceIterator fi=tri::Allocator<MeshType>::AddFaces(cm,nv-2);
	// 	CMeshO::FaceIterator fi=tri::Allocator<CMeshO>::AddFaces(cm,indices.size()/3);
	// 	for (size_t i=0; i<indices.size(); i+=3,++fi)
	// 	{
	// 		(*&fi)->V(0)=&cm.vert[ indices[i+0] ];
	// 		(*&fi)->V(1)=&cm.vert[ indices[i+1] ];
	// 		(*&fi)->V(2)=&cm.vert[ indices[i+2] ];
	// 	}
	// //	tri::Clean<CMeshO>::RemoveDuplicateVertex(cm);
	// 	UpdateBounding<CMeshO>::Box(cm);
}

void EditRepairPlugin::slotCut(){
	if (!importset)
	{
		QMessageBox msgBoxWarn;
		msgBoxWarn.setText("Please import local coordinate!!!");
		msgBoxWarn.exec();
		return ;
	}
	for (std::vector<Circle *>::iterator it = circles.begin(); it != circles.end(); it ++){
		if (NULL != *it) 
		{
			md->delMesh((*it)->m);
			delete *it; 
			*it = NULL;
		}
	}
	circles.clear();
	vcg::Point3f planeAxis = zcoord;
	vcg::Point3f planeCenter = barycentre;
	vcg::Plane3f slicingPlane;
	//确定mesh最高点和最低点
	vector<Point3f> allpoint;
	CMeshO::VertexIterator vi;
	for(vi = mesh->cm.vert.begin();vi!= mesh->cm.vert.end();vi++){
		CVertexO *vii = &(*vi);
		vcg::Point3f c = vii->P();
		allpoint.push_back(c);
	}
	//先对图形进行平移
	for (int i=0;i<allpoint.size();i++)
	{
		allpoint[i].Z()=allpoint[i].X()*zcoord.X()+allpoint[i].Y()*zcoord.Y()+allpoint[i].Z()*zcoord.Z()-(barycentre.X()*zcoord.X()+barycentre.Y()*zcoord.Y()+barycentre.Z()*zcoord.Z());
	}
	//找出最大点和最小点
	float minfloat=allpoint[0].Z();
	float maxfloat=allpoint[0].Z();
	for (int i=1;i<allpoint.size();i++)
	{
		if (minfloat>allpoint[i].Z())
		{
			minfloat=allpoint[i].Z();
		}
		if (maxfloat<allpoint[i].Z())
		{
			maxfloat=allpoint[i].Z();
		}
	}
	float perStep = stagesDialog->ui.doubleSpinBoxPerStep->value();//确定步长
	float planeOffset = stagesDialog->ui.doubleSpinBoxPlaneOffset->value();//确定初始偏移量
	int cut_times = stagesDialog->ui.spinBoxCutTimes->value();//确定切割次数

	//自动确定步长，褚
	planeOffset=maxfloat;
	perStep=(minfloat-maxfloat)/cut_times;


	for (int i = 0; i<cut_times; i++)
	{
		planeCenter = barycentre+planeAxis*(planeOffset + i*perStep);//得到平面过的点

		slicingPlane.Init(planeCenter,planeAxis);//初始化平面
		MeshModel* orig = mesh;
		RenderMode rm;
		rm.drawMode = GLW::DMWire;
		MeshModel* cap= md->addNewMesh("","",false);
		//MeshModel* cap= new MeshModel(md, "", "");
		IntersectionPlaneMesh(orig->cm, slicingPlane, cap->cm );
		//vcg::IntersectionPlaneMesh<CMeshO, CMeshO, float>(orig->cm, slicingPlane, cap->cm );
		//tri::Clean<CMeshO>::RemoveDuplicateVertex(cap->cm);
		Circle* c = new Circle(cap,barycentre,xcoord,ycoord,zcoord,planeOffset + i*perStep);
		circles.push_back(c);

	}
	cout<<"*******************Cut Finished*******************"<<endl;
	QMessageBox msgBox;
	msgBox.setText("Finished!");
	msgBox.exec();
	if(CurrentMode() != TranMode::CUT)
		SetMode(TranMode::CUT);
	gla->update();
}

void EditRepairPlugin::slotCut1(){
	vcg::Point3f planeAxis = zcoord;
	vcg::Point3f planeCenter = barycentre;
	vcg::Plane3f slicingPlane;
	slicingPlane.Init(planeCenter,planeAxis);
	MeshModel* orig = mesh;

	MeshModel* cap= md->addNewMesh("","",false);
	IntersectionPlaneMesh(orig->cm, slicingPlane, cap->cm );
	//	tri::Clean<CMeshO>::RemoveDuplicateVertex(cap->cm);

	// 	MeshModel* cap2= md.addNewMesh("","",false);
	// 	tri::CapEdgeMesh(cap->cm, cap2->cm);
	// 	cap2->UpdateBoxAndNormals();


	cout<<"slotCut"<<endl;
	gla->update();
}

void EditRepairPlugin::slotMerge(){
	mesh = md->meshList[0];
	mesh2 = md->meshList[1];
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	std::vector<CVertexO> vv_mm;
	std::vector<CFaceO> vf_mm;
	std::vector<CVertexO> vv_mm2;
	std::vector<CFaceO> vf_mm2;
	std::vector<int> VertexId(mesh->cm.vert.size());//mesh_origin为原始已有的mesh，VertexId在这里作为一个索引，利用它可以找到编辑后的mesh所对应的新的索引，有点绕，不过不难理解
	int numvert = 0;
	std::vector<int> VertexId2(mesh2->cm.vert.size());//mesh_origin为原始已有的mesh，VertexId在这里作为一个索引，利用它可以找到编辑后的mesh所对应的新的索引，有点绕，不过不难理解
	int numvert2 = 0;

	for (vi=mesh->cm.vert.begin();vi!=mesh->cm.vert.end();vi++){//保存mesh_origin的顶点
		if (!(*vi).IsD())
		{
			VertexId[vi-mesh->cm.vert.begin()]=numvert;
			vv_mm.push_back(*vi);
			numvert++;
		}
	}
	for(fi=mesh->cm.face.begin(); fi!=mesh->cm.face.end(); ++fi)//保存mesh_origin的面片
	{
		if( !(*fi).IsD() ){
			vf_mm.push_back(*fi);
		}
	}

	for (vi=mesh2->cm.vert.begin();vi!=mesh2->cm.vert.end();vi++){//保存mesh_origin的顶点
		if (!(*vi).IsD())
		{
			VertexId2[vi-mesh2->cm.vert.begin()]=numvert2;
			vv_mm2.push_back(*vi);
			numvert2++;
		}
	}
	for(fi=mesh2->cm.face.begin(); fi!=mesh2->cm.face.end(); ++fi)//保存mesh_origin的面片
	{
		if( !(*fi).IsD() ){
			vf_mm2.push_back(*fi);
		}
	}

	MeshModel *mesh_copy = md->addNewMesh("mesh_copy","",false);//新建mesh_copy 
	mesh_copy->cm.Clear();
	int VN,FN;
	VN=vv_mm.size() + vv_mm2.size();
	FN=vf_mm.size() + vf_mm2.size();
	vcg::tri::Allocator<CMeshO>::AddVertices(mesh_copy->cm,VN);
	vcg::tri::Allocator<CMeshO>::AddFaces(mesh_copy->cm,FN);
	CMeshO::VertexPointer *ivp = new CMeshO::VertexPointer[VN];
	int cnt=0;
	vi= mesh_copy->cm.vert.begin();
	std::vector<CVertexO>::const_iterator itv;
	for(itv=vv_mm.begin();itv!=vv_mm.end();itv++){//复制顶点
		ivp[cnt]=&*vi;
		(*vi).P()= (*itv).P();
		++vi;
		++cnt;
	}

	for(itv=vv_mm2.begin();itv!=vv_mm2.end();itv++){//复制顶点
		ivp[cnt]=&*vi;
		(*vi).P()= (*itv).P();
		++vi;
		++cnt;
	}

	fi=mesh_copy->cm.face.begin();
	std::vector<CFaceO>::const_iterator itf;//130
	//CMeshO::FaceIterator itf;//133
	for (itf = vf_mm.begin(); itf!= vf_mm.end(); itf++)//复制面片
	{
		for(int k=0;k<(*itf).VN();k++)
		{
			int vInd = VertexId[(*itf).V(k)-&*(mesh->cm.vert.begin())];//确定顶点的索引值
			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}
	for (itf = vf_mm2.begin(); itf!= vf_mm2.end(); itf++)//复制面片
	{
		for(int k=0;k<(*itf).VN();k++)
		{
			int vInd = VertexId2[(*itf).V(k)-&*(mesh2->cm.vert.begin())] + vv_mm.size();//确定顶点的索引值
			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}

	vcg::tri::UpdateBounding<CMeshO>::Box(mesh_copy->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_copy->cm);
	cout<<"slotMerge"<<endl;
	gla->update();
}

void EditRepairPlugin::slotNormalize(){
	CMeshO::VertexIterator vi;

	for (vi=mesh->cm.vert.begin();vi!=mesh->cm.vert.end();vi++){
		if (!(*vi).IsD())
		{
			vcg::Point3f newp = vcg::Point3f((vi->P()-barycentre)*xcoord, (vi->P()-barycentre)*ycoord, (vi->P()-barycentre)*zcoord);
			vi->P() = newp;
		}
	}

	vcg::tri::UpdateBounding<CMeshO>::Box(mesh->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh->cm);
	cout<<"slotNormalize"<<endl;
	gla->update();
}

void EditRepairPlugin::slotExportPoints()
{
	QString filetype = ".txt";
	QString filename = QFileDialog::getSaveFileName(this->gla,tr("Save file"),"circles","*"+filetype);
	if (filename != "")
	{
		QFile f(filename);
		if(!f.open(QIODevice::WriteOnly | QIODevice::Text))
		{
			cout << "Open failed." << endl;
			return;
		}

		QTextStream txtOutput(&f);

		for (int i = 0; i < circles.size(); i++)
		{
			txtOutput <<"#circle "<<i<<endl;
			vcg::Point3f new_barycenter = circles[i]->barycenter - barycentre;
			//txtOutput <<"#center "<<circles[i]->barycenter.X()<<" "<<circles[i]->barycenter.Y()<<" "<<circles[i]->barycenter.Z()<<" " <<endl;
			txtOutput <<"#center "<<new_barycenter * xcoord<<" "<<new_barycenter * ycoord<<" "<<new_barycenter * zcoord<<" " <<endl;
			/*		txtOutput <<"#offset "<<circles[i].offset<<endl;*/
			vector<int> ids=faceids[i];
			for (int k = 0; k < circles[i]->m->cm.vert.size(); k++)
			{
				//int faceid=(int)circles[i]->m->cm.vert[k].P().Z();
				//circles[i]->m->cm.vert[k].P().Z()=(circles[i]->m->cm.vert[k].P().Z()-(int)(circles[i]->m->cm.vert[k].P().Z()))*100;
				vcg::Point3f new_circlepoint = circles[i]->m->cm.vert[k].P() - barycentre;
				//txtOutput <<"#circlepoint "<<circles[i]->m->cm.vert[k].P().X()<<" "<<circles[i]->m->cm.vert[k].P().Y()<<" "<<circles[i]->m->cm.vert[k].P().Z()<<endl;
				txtOutput <<"#circlepoint "<<" "<<new_circlepoint * xcoord<<" "<<new_circlepoint * ycoord<<" "<<new_circlepoint * zcoord<<" "<<ids[k]<<endl;
			}
		}
		f.close();
	}
}

void EditRepairPlugin::slotExportPoints2()
{
	QString filetype = ".txt";
	QString filename = QFileDialog::getSaveFileName(this->gla,tr("Save file"),"circles","*"+filetype);
	if (filename != "")
	{
		QFile f(filename);
		if(!f.open(QIODevice::WriteOnly | QIODevice::Text))
		{
			cout << "Open failed." << endl;
			return;
		}

		QTextStream txtOutput(&f);

		for (int i = 0; i < circles.size(); i++)
		{
			//txtOutput <<"#circle "<<i<<endl;
			vcg::Point3f new_barycenter = circles[i]->barycenter - barycentre;
			//txtOutput <<"#center "<<circles[i]->barycenter.X()<<" "<<circles[i]->barycenter.Y()<<" "<<circles[i]->barycenter.Z()<<" " <<endl;
			txtOutput <<new_barycenter * xcoord<<" "<<new_barycenter * ycoord<<" "<<new_barycenter * zcoord<<" " <<endl;
			/*		txtOutput <<"#offset "<<circles[i].offset<<endl;*/
			for (int k = 0; k < circles[i]->m->cm.vert.size(); k++)
			{
				vcg::Point3f new_circlepoint = circles[i]->m->cm.vert[k].P() - barycentre;
				//txtOutput <<"#circlepoint "<<circles[i]->m->cm.vert[k].P().X()<<" "<<circles[i]->m->cm.vert[k].P().Y()<<" "<<circles[i]->m->cm.vert[k].P().Z()<<endl;
				txtOutput <<new_circlepoint * xcoord<<" "<<new_circlepoint * ycoord<<" "<<new_circlepoint * zcoord<<endl;
			}
		}
		f.close();
	}
}

void EditRepairPlugin::slotErrorMap()
{
	for(int i = 0; i < md->meshList.size(); i++){
		QString qstr = md->meshList[i]->shortName();
		char*  ch;
		QByteArray ba = qstr.toLatin1();
		ch=ba.data();
		if (ch[3] == 'r')
		{
			mesh2 = md->meshList[i];
		}else if (ch[3] == '.')
		{
			mesh = md->meshList[i];
		}
	}

	clock_t start,finish;
	double totaltime;
	start=clock();
	// 	float threshold = stagesDialog->ui.doubleSpinBox->value();
	// 	float threshold2 = stagesDialog->ui.doubleSpinBox_2->value();
	// 	float threshold3 = stagesDialog->ui.doubleSpinBox_3->value();
	float threshold4 = stagesDialog->ui.doubleSpinBox_4->value();
	stagesDialog->ui.doubleSpinBox->setValue(threshold4/4*1);
	stagesDialog->ui.doubleSpinBox_2->setValue(threshold4/4*2);
	stagesDialog->ui.doubleSpinBox_3->setValue(threshold4/4*3);

	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	float t;
	float ei;
	float dsquare1 = 0.0;
	float dsquare2 = 0.0;
	float ei2;
	float dsquare12 = 0.0;
	float dsquare22 = 0.0;
	int numsnumzero = 0;

	AIndex gIndex;
	gIndex.Set(mesh->cm.face.begin(), mesh->cm.face.end());

	mesh_tem_d = mesh2;//Qianjiahong
	float maxerror = -1.0;
	for(vit = mesh_tem_d->cm.vert.begin(); vit != mesh_tem_d->cm.vert.end(); ++vit){
		if (!vit->IsD())
		{
			ProLine pl;
			pl.lp = vit->P();
			pl.ldir = vit->N();//(vit->P()-barycentre)/Distance(vit->P(),barycentre);

			const bool TEST_BACK_FACES = true;
			vcg::RayTriangleIntersectionFunctor<TEST_BACK_FACES> rayIntersector;
			const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();
			vcg::Ray3<AIndex::ScalarType, false> ray(pl.lp, pl.ldir);

			AIndex::ObjPtr isectFace;
			AIndex::ScalarType rayT;
			vcg::EmptyClass a;
			isectFace = gIndex.DoRay(rayIntersector, a , ray, maxDist, rayT);

			if (fabs(rayT) > 100)
			{
				rayT = 0.0;
			}
			// 			if ( < threshold)
			// 			{
			// 				vit->C() = Color4b::Blue;
			// 			}else if (fabs(rayT) < threshold2)
			// 			{
			// 				vit->C() = Color4b::Cyan;
			// 			}else if (fabs(rayT) < threshold3)
			// 			{
			// 				vit->C() = Color4b::Green;
			// 			}else if (fabs(rayT) < threshold4)
			// 			{
			// 				vit->C() = Color4b::Yellow;
			// 			}else
			// 			{
			// 				vit->C() = Color4b::Red;
			// 				cout<< "Red "<<fabs(rayT)<<endl;
			// 			}

			vcg::Color4f c;
			if (fabs(rayT)/threshold4 > 1.0)
			{
				c.SetHSVColor(0,1,1);
			}else{
				c.SetHSVColor((1.0-fabs(rayT)/threshold4)*4/6,1,1);
			}

			(*vit).C().X() = c[0];
			(*vit).C().Y() = c[1];
			(*vit).C().Z() = c[2];

			if (fabs(rayT) > maxerror)
			{
				maxerror = fabs(rayT);
			}
			if (fabs(rayT) > 0.0001)
			{
				dsquare12 += rayT*rayT;
				numsnumzero++;
			}
			dsquare1 += rayT*rayT;
		}
	}
	stagesDialog->ui.errorinfo->setText(QString("MaxError: ")+QString("%1").arg(maxerror));
	stagesDialog->ui.errorinfo->adjustSize();

	maxerror = -1.0;
	AIndex gIndex2;
	gIndex2.Set(mesh_tem_d->cm.face.begin(), mesh_tem_d->cm.face.end());
	for(vit = mesh->cm.vert.begin(); vit != mesh->cm.vert.end(); ++vit){
		if (!vit->IsD())
		{
			ProLine pl;
			pl.lp = vit->P();
			pl.ldir = vit->N();//(vit->P()-barycentre)/Distance(vit->P(),barycentre);

			const bool TEST_BACK_FACES = true;
			vcg::RayTriangleIntersectionFunctor<TEST_BACK_FACES> rayIntersector;
			const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();
			vcg::Ray3<AIndex::ScalarType, false> ray(pl.lp, pl.ldir);

			AIndex::ObjPtr isectFace;
			AIndex::ScalarType rayT;
			vcg::EmptyClass a;
			isectFace = gIndex2.DoRay(rayIntersector, a , ray, maxDist, rayT);

			if (fabs(rayT) > 100)
			{
				rayT = 0.0;
			}
			if (fabs(rayT) > maxerror)
			{
				maxerror = fabs(rayT);
			}
			if (fabs(rayT) > 0.0001)
			{
				dsquare22 += rayT*rayT;
				numsnumzero++;
			}
			dsquare2 += rayT*rayT;
		}
	}
	cout<<"dsquare1 "<<dsquare1<<" dsquare2 "<<dsquare2<<" maxerror "<<maxerror<<endl;
	cout<<"dsquare12 "<<dsquare12<<" dsquare22 "<<dsquare22<<" numsnumzero "<<numsnumzero<<endl;
	ei = (dsquare1 + dsquare2)/(mesh->cm.vert.size() + mesh_tem_d->cm.vert.size());
	ei2 = (dsquare12 + dsquare22)/(numsnumzero);
	stagesDialog->ui.errorinfo_2->setText(QString("Approximation Error(Ei): ")+QString("%1").arg(ei)+QString("\nEliminate Zero Point Ei: ")+QString("%1").arg(ei2));
	stagesDialog->ui.errorinfo_2->adjustSize();

	mesh_tem_d->updateDataMask(MeshModel::MM_VERTCOLOR);

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n此程序的运行时间为"<<totaltime<<endl;

	if(CurrentMode() != TranMode::ERRORMAP)
		SetMode(TranMode::ERRORMAP);

	gla->update();
}

void EditRepairPlugin::slotErrorMap2()
{
	clock_t start,finish;
	double totaltime;
	start=clock();
	float threshold = stagesDialog->ui.doubleSpinBox->value();
	float threshold2 = stagesDialog->ui.doubleSpinBox_2->value();
	float threshold3 = stagesDialog->ui.doubleSpinBox_3->value();
	float threshold4 = stagesDialog->ui.doubleSpinBox_4->value();

	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	float t;
	float ei;
	float dsquare1 = 0.0;
	float dsquare2 = 0.0;
	float ei2;
	float dsquare12 = 0.0;
	float dsquare22 = 0.0;
	int numsnumzero = 0;

	AIndex gIndex;
	gIndex.Set(mesh->cm.face.begin(), mesh->cm.face.end());

	mesh_tem_d = mesh2;//Qianjiahong
	float maxerror = -1.0;
	for(vit = mesh_tem_d->cm.vert.begin(); vit != mesh_tem_d->cm.vert.end(); ++vit){
		if (!vit->IsD())
		{
			ProLine pl;
			pl.lp = vit->P();
			pl.ldir = vit->N();//(vit->P()-barycentre)/Distance(vit->P(),barycentre);

			const bool TEST_BACK_FACES = true;
			vcg::RayTriangleIntersectionFunctor<TEST_BACK_FACES> rayIntersector;
			const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();
			vcg::Ray3<AIndex::ScalarType, false> ray(pl.lp, pl.ldir);

			AIndex::ObjPtr isectFace;
			AIndex::ScalarType rayT;
			vcg::EmptyClass a;
			isectFace = gIndex.DoRay(rayIntersector, a , ray, maxDist, rayT);

			if (fabs(rayT) > 100)
			{
				rayT = 0.0;
			}
			if (fabs(rayT) < threshold)
			{
				vit->C() = Color4b::Blue;
			}else if (fabs(rayT) < threshold2)
			{
				vit->C() = Color4b::Cyan;
			}else if (fabs(rayT) < threshold3)
			{
				vit->C() = Color4b::Green;
			}else if (fabs(rayT) < threshold4)
			{
				vit->C() = Color4b::Yellow;
			}else
			{
				vit->C() = Color4b::Red;
				cout<< "Red "<<fabs(rayT)<<endl;
			}

			if (fabs(rayT) > maxerror)
			{
				maxerror = fabs(rayT);
			}
			if (fabs(rayT) > 0.0001)
			{
				dsquare12 += rayT*rayT;
				numsnumzero++;
			}
			dsquare1 += rayT*rayT;
		}
	}
	stagesDialog->ui.errorinfo->setText(QString("MaxError: ")+QString("%1").arg(maxerror));
	stagesDialog->ui.errorinfo->adjustSize();

	maxerror = -1.0;
	AIndex gIndex2;
	gIndex2.Set(mesh_tem_d->cm.face.begin(), mesh_tem_d->cm.face.end());
	for(vit = mesh->cm.vert.begin(); vit != mesh->cm.vert.end(); ++vit){
		if (!vit->IsD())
		{
			ProLine pl;
			pl.lp = vit->P();
			pl.ldir = vit->N();//(vit->P()-barycentre)/Distance(vit->P(),barycentre);

			const bool TEST_BACK_FACES = true;
			vcg::RayTriangleIntersectionFunctor<TEST_BACK_FACES> rayIntersector;
			const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();
			vcg::Ray3<AIndex::ScalarType, false> ray(pl.lp, pl.ldir);

			AIndex::ObjPtr isectFace;
			AIndex::ScalarType rayT;
			vcg::EmptyClass a;
			isectFace = gIndex2.DoRay(rayIntersector, a , ray, maxDist, rayT);

			if (fabs(rayT) > 100)
			{
				rayT = 0.0;
			}
			if (fabs(rayT) > maxerror)
			{
				maxerror = fabs(rayT);
			}
			if (fabs(rayT) > 0.0001)
			{
				dsquare22 += rayT*rayT;
				numsnumzero++;
			}
			dsquare2 += rayT*rayT;
		}
	}
	cout<<"dsquare1 "<<dsquare1<<" dsquare2 "<<dsquare2<<" maxerror "<<maxerror<<endl;
	cout<<"dsquare12 "<<dsquare12<<" dsquare22 "<<dsquare22<<" numsnumzero "<<numsnumzero<<endl;
	ei = (dsquare1 + dsquare2)/(mesh->cm.vert.size() + mesh_tem_d->cm.vert.size());
	ei2 = (dsquare12 + dsquare22)/(numsnumzero);
	stagesDialog->ui.errorinfo_2->setText(QString("Approximation Error(Ei): ")+QString("%1").arg(ei)+QString("\nEliminate Zero Point Ei: ")+QString("%1").arg(ei2));
	stagesDialog->ui.errorinfo_2->adjustSize();

	mesh_tem_d->updateDataMask(MeshModel::MM_VERTCOLOR);

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n此程序的运行时间为"<<totaltime<<endl;

	if(CurrentMode() != TranMode::ERRORMAP)
		SetMode(TranMode::ERRORMAP);

	gla->update();
}

void EditRepairPlugin::slotErrorMap1()
{
	clock_t start,finish;
	double totaltime;
	start=clock();

	CMeshO::VertexIterator vit;
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	float t;

	AIndex gIndex;
	gIndex.Set(mesh->cm.face.begin(), mesh->cm.face.end());

	mesh_tem_d = mesh2;//Qianjiahong
	for(vit = mesh_tem_d->cm.vert.begin(); vit != mesh_tem_d->cm.vert.end(); ++vit){
		if (!vit->IsD())
		{
			ProLine pl;
			pl.lp = vit->P();
			pl.ldir = (vit->P()-barycentre)/Distance(vit->P(),barycentre);

			const bool TEST_BACK_FACES = true;
			vcg::RayTriangleIntersectionFunctor<TEST_BACK_FACES> rayIntersector;
			const AIndex::ScalarType maxDist = std::numeric_limits<AIndex::ScalarType>::max();
			vcg::Ray3<AIndex::ScalarType, false> ray(barycentre, pl.ldir);

			AIndex::ObjPtr isectFace;
			AIndex::ScalarType rayT;
			vcg::EmptyClass a;
			isectFace = gIndex.DoRay(rayIntersector, a , ray, maxDist, rayT);

			if (fabs(fabs(rayT) - Distance(vit->P(),barycentre)) < 0.04)
			{
				vit->C() = Color4b::Blue;
			}else if (fabs(fabs(rayT) - Distance(vit->P(),barycentre)) < 0.08)
			{
				vit->C() = Color4b::Cyan;
			}else if (fabs(fabs(rayT) - Distance(vit->P(),barycentre)) < 0.12)
			{
				vit->C() = Color4b::Green;
			}else if (fabs(fabs(rayT) - Distance(vit->P(),barycentre)) < 0.17)
			{
				vit->C() = Color4b::Yellow;
			}else if(fabs(fabs(rayT) - Distance(vit->P(),barycentre)) < 4)
			{
				vit->C() = Color4b::Red;
				cout<< "Red "<<fabs(fabs(rayT) - Distance(vit->P(),barycentre))<<endl;
			}else {
				vit->C() = Color4b::Blue;
			}
			if (isectFace != 0) {

				//vit->P() = barycentre + pl.ldir*rayT;
				//vit->C() = Color4b::Red;
				//cout<< "fabs(fabs(rayT) - Distance(vit->P(),barycentre)) "<<fabs(fabs(rayT) - Distance(vit->P(),barycentre))<<endl;
				continue;
			}
		}
	}
	mesh_tem_d->updateDataMask(MeshModel::MM_VERTCOLOR);

	finish=clock();
	totaltime=(double)(finish-start);
	cout<<"\n此程序的运行时间为"<<totaltime<<endl;

	gla->update();
}

void EditRepairPlugin::DrawColorBar(QPainter *p)
{
	p->endNativePainting();
	p->save();
	QColor color;  
	QRect section;

	float colorBarLength=800;//设置颜色条的长度  
	int x = 300;
	int y = 50;
	int barwidth = 30;

	for(int i=0;i<=colorBarLength/6*4;i++)// hsv  
	{  
		color.setHsvF(i/colorBarLength,1,1);  
		//section.setRect(x+50,colorBarLength+y-i*1,20,1);  
		section.setRect(x+50,y+i*1,barwidth,1);  
		p->fillRect(section,color);  
	}  

	//---------设置边框--------------//  
	//刻度值的绘制可以自己设计，使用drawText函数即可,刻度的绘制可以使用drawLine函数  
	p->drawLine(x+40, y + colorBarLength/6*0, x+60, y + colorBarLength/6*0); 
	p->drawLine(x+40, y + colorBarLength/6*1, x+60, y + colorBarLength/6*1);
	p->drawLine(x+40, y + colorBarLength/6*2, x+60, y + colorBarLength/6*2);
	p->drawLine(x+40, y + colorBarLength/6*3, x+60, y + colorBarLength/6*3);
	p->drawLine(x+40, y + colorBarLength/6*4, x+60, y + colorBarLength/6*4);
	p->drawRect(x+50,y,barwidth,colorBarLength/6*4);

	p->setRenderHint(QPainter::TextAntialiasing);
	p->setPen(Qt::black);
	QFont qFont;
	qFont.setStyleStrategy(QFont::NoAntialias);
	qFont.setFamily("Helvetica");
	qFont.setPixelSize(18);
	p->setFont(qFont);
	//float barHeight = qFont.pixelSize()*5;
	//QFontMetrics metrics = p->fontMetrics();
	//int border = qMax(4, metrics.leading());

	QString col1Text =" ";
	QRect Column_1(x+50,40,barwidth,10);
	p->drawText(Column_1, Qt::AlignLeft | Qt::TextWordWrap, col1Text);

	p->drawText(x+50-45,y + colorBarLength/6.0 * 0.0 + 6,QString("%1").arg(stagesDialog->ui.doubleSpinBox_4->value()));  
	p->drawText(x+50-45,y + colorBarLength/6.0 * 1.0 + 6,QString("%1").arg(stagesDialog->ui.doubleSpinBox_3->value()));  
	p->drawText(x+50-45,y + colorBarLength/6.0 * 2.0 + 6,QString("%1").arg(stagesDialog->ui.doubleSpinBox_2->value()));  
	p->drawText(x+50-45,y + colorBarLength/6.0 * 3.0 + 6,QString("%1").arg(stagesDialog->ui.doubleSpinBox->value()));  
	p->drawText(x+50-45,y + colorBarLength/6.0 * 4.0 + 6,QString("0.00")); 
	p->restore();
	p->beginNativePainting(); 
}

void EditRepairPlugin::slotInvertNormal()
{
	CMeshO::VertexIterator vi;
	for(vi = mesh->cm.vert.begin(); vi != mesh->cm.vert.end(); ++vi){
		if (!vi->IsD()){
			vi->N() = -vi->N();
		}
	}
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh->cm);
	gla->update();
}

void EditRepairPlugin::slotInvertNormalFace()
{
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	std::vector<CVertexO> vv_mm;
	std::vector<CFaceO> vf_mm;
	std::vector<int> VertexId(mesh->cm.vert.size());
	int numvert = 0;

	for (vi=mesh->cm.vert.begin();vi!=mesh->cm.vert.end();vi++){
		if (!(*vi).IsD())
		{
			VertexId[vi-mesh->cm.vert.begin()]=numvert;
			vv_mm.push_back(*vi);
			numvert++;
		}
	}
	for(fi=mesh->cm.face.begin(); fi!=mesh->cm.face.end(); ++fi)
	{
		if( !(*fi).IsD() ){
			vf_mm.push_back(*fi);
		}
	}

	mesh_copy = md->addNewMesh("mesh_copy","",false);

	mesh_copy->cm.Clear();
	int VN,FN;
	VN=vv_mm.size();
	FN=vf_mm.size();
	vcg::tri::Allocator<CMeshO>::AddVertices(mesh_copy->cm,VN);
	vcg::tri::Allocator<CMeshO>::AddFaces(mesh_copy->cm,FN);
	CMeshO::VertexPointer *ivp = new CMeshO::VertexPointer[VN];
	int cnt=0;
	vi= mesh_copy->cm.vert.begin();
	std::vector<CVertexO>::const_iterator itv;
	for(itv=vv_mm.begin();itv!=vv_mm.end();itv++){
		ivp[cnt]=&*vi;
		(*vi).P()= (*itv).P();
		++vi;
		++cnt;
	}
	fi=mesh_copy->cm.face.begin();
	std::vector<CFaceO>::const_iterator itf;
	for (itf = vf_mm.begin(); itf!= vf_mm.end(); itf++)
	{
		for(int k=0;k<(*itf).VN();k++)
		{
			int vInd = VertexId[(*itf).V(k)-&*(mesh->cm.vert.begin())];
			(*fi).V((*itf).VN()-k-1)= ivp[vInd];
		}
		++fi;
	}

	vcg::tri::UpdateBounding<CMeshO>::Box(mesh_copy->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh_copy->cm);
	gla->update();
}

void EditRepairPlugin::slotFindNeastestPoints()
{
	for(int i = 0; i < md->meshList.size(); i++){
		QString qstr = md->meshList[i]->shortName();
		char*  ch;
		QByteArray ba = qstr.toLatin1();
		ch=ba.data();
		if (ch[3] == '.')
		{
			mesh = md->meshList[i];
		}else{
			mesh2 = md->meshList[i];
		}
	}

	mesh->updateDataMask(MeshModel::MM_FACEFLAGBORDER);
	mesh2->updateDataMask(MeshModel::MM_FACEFLAGBORDER);

// 	angleV1 = vcg::Point3f(24.8964,6.2205,-1.64945);//UL1 hole2
// 	angleV2 = vcg::Point3f(21.7318,4.25893,-1.12725);
// 	angleV1 = vcg::Point3f(23.8603,9.44466,-2.537);//UL1 hole1
// 	angleV2 = vcg::Point3f(19.8984,7.74382,-2.32539);

// 	angleV1 = vcg::Point3f(-4.12969, 2.98214, -1.43422);//UL7 hole1
// 	angleV2 = vcg::Point3f(3.60245, 5.00985, -1.32964);

// 	angleV1 = vcg::Point3f(-3.29796, -4.51897, -1.14542);//UL7 hole2
// 	angleV2 = vcg::Point3f(4.3167, -2.9806, -1.23865);

	QString dir = QFileDialog::getOpenFileName(this->gla, tr("Open Set File..."),
		currentfilepath,
		tr("Setting (*.txt)"));

	int index = dir.lastIndexOf("/");
	currentfilepath = dir;
	currentfilepath.remove(index,dir.length()-index);

	std::string str = dir.toStdString();
	ifstream inputs;
	std::string oneline;
	inputs.open(str.c_str());
	int i = 0;
	if (inputs.is_open())
	{
		int i = -1;
		while(getline(inputs, oneline)){
			istringstream strStream(oneline);

			double ax,ay,az;
			strStream >> ax >> ay >> az;
			
			if (i == 0)
			{
				angleV1.X() = ax;
				angleV1.Y() = ay;
				angleV1.Z() = az;
			}else if(i == 1){
				angleV2.X() = ax;
				angleV2.Y() = ay;
				angleV2.Z() = az;
			}
			
			i++;
		}
		inputs.close();

		float max = -100000;
		float min = 100000;
		for (int i = 0 ; i < mesh->cm.vert.size(); i++)
		{
			float zpro = (mesh->cm.vert[i].P() - barycentre)* zcoord;
			if (zpro > max)
			{
				max = zpro;
			}
			if (zpro < min)
			{
				min = zpro;
			}
		}
		stagesDialog->ui.toothinfo->setText(QString("Max: ")+QString("%1").arg(max)+QString("; Min: ") + QString("%1").arg(min));
		stagesDialog->ui.toothinfo->adjustSize();


		if(CurrentMode() != TranMode::IMPORTCOORD)
			SetMode(TranMode::IMPORTCOORD);
		importset = true;
	}


	float minDis1 = Distance(angleV1, mesh->cm.vert[0].P());
	flag1 = 0;
	for(int i = 1; i< mesh->cm.vert.size(); i++){
		if(mesh->cm.vert[i].IsB()){
			if(Distance(angleV1, mesh->cm.vert[i].P())<minDis1){
				minDis1 = Distance(angleV1, mesh->cm.vert[i].P());
				flag1 = i;
			}
		}
	}
	cout<<"minDis1 "<<minDis1<<endl;
	cout<<"flag1 "<<flag1<<endl;

	minDis1 = Distance(angleV1, mesh2->cm.vert[0].P());
	flagH1 = 0;
	for(int i = 1; i< mesh2->cm.vert.size(); i++){
		if(mesh2->cm.vert[i].IsB()){
			if(Distance(angleV1, mesh2->cm.vert[i].P())<minDis1){
				minDis1 = Distance(angleV1, mesh2->cm.vert[i].P());
				flagH1 = i;
			}
		}
	}
	cout<<"minDis1 "<<minDis1<<endl;
	cout<<"flagH1 "<<flagH1<<endl;

	float minDis2 = Distance(angleV2, mesh->cm.vert[0].P());
	flag2 = 0;
	for(int i = 1; i< mesh->cm.vert.size(); i++){
		if(mesh->cm.vert[i].IsB()){
			if(Distance(angleV2, mesh->cm.vert[i].P())<minDis2){
				minDis2 = Distance(angleV2, mesh->cm.vert[i].P());
				flag2 = i;
			}
		}
	}
	cout<<"minDis2 "<<minDis2<<endl;
	cout<<"flag2 "<<flag2<<endl;
	minDis2 = Distance(angleV2, mesh2->cm.vert[0].P());
	flagH2 = 0;
	for(int i = 1; i< mesh2->cm.vert.size(); i++){
		if(mesh2->cm.vert[i].IsB()){
			if(Distance(angleV2, mesh2->cm.vert[i].P())<minDis2){
				minDis2 = Distance(angleV2, mesh2->cm.vert[i].P());
				flagH2 = i;
			}
		}
	}
	cout<<"minDis2 "<<minDis2<<endl;
	cout<<"flagH2 "<<flagH2<<endl;
	if(CurrentMode() != TranMode::FINDNESTP)
		SetMode(TranMode::FINDNESTP);
	gla->update();
}


void EditRepairPlugin::slotExtendBorder()
{
	if (1)
	{
		mesh->updateDataMask(MeshModel::MM_VERTFACETOPO);
		mesh->updateDataMask(MeshModel::MM_FACEFACETOPO);

		teethBordersInOrder.clear();
		teethBordersInOrder.push_back(flag1);

		int findtheFirst = 0;
		CVertexO* V1= NULL;

		CVertexO *vi = &(mesh->cm.vert[flag1]);
		vcg::face::JumpingPos<CFaceO> pos(vi->VFp(),vi);
		CVertexO* firstV = pos.VFlip();
		CVertexO* tempV= NULL;
		CVertexO* tempVLast= NULL;
		CVertexO* tempVLastLast= NULL;
		Point3f p1,p2;
		do 
		{
			pos.NextE();
			tempV = pos.VFlip();
			// 			cout<<"vi:"<< vi->P().X()<< vi->P().Y()<< vi->P().Z()<<endl;
			// 			cout<<"firstV:"<< firstV->P().X()<< firstV->P().Y()<< firstV->P().Z()<<endl;
			// 			cout<<"tempV:"<< tempV->P().X()<< tempV->P().Y()<< tempV->P().Z()<<endl;
			if (tempVLast != NULL && tempVLastLast != NULL && tempVLast->IsB()){
				Point3f ptemp = tempVLast->P()-vi->P();
				if (ptemp*(zcoord)>0)//沿Y轴负方向(待改进，比较大的那一个)
				{
					if (tempV == tempVLastLast)
					{
						break;
						//cout<<"slotSortBorder2"<<endl;
					}
				}
			}
			tempVLastLast = tempVLast;
			tempVLast = tempV;
		} while (1);

		cout<<"slotSortBorder3"<<endl;

		CVertexO* Vp= vi;

		int num = 1;
		while (tempVLast != &mesh->cm.vert[flag2])//到另一个边界点
		{
			teethBordersInOrder.push_back(tempVLast-&*(mesh->cm.vert.begin()));
			CVertexO *vi0 = tempVLast;
			vcg::face::JumpingPos<CFaceO> pos0(vi0->VFp(),vi0);
			CVertexO* firstV0 = pos0.VFlip();
			CVertexO* tempV0= NULL;//访问V的一环邻域
			CVertexO* tempV0Last= NULL;
			CVertexO* tempV0LastLast= NULL;
			do 
			{

				pos0.NextE();
				tempV0 = pos0.VFlip();
				if (tempV0Last != NULL && tempV0LastLast != NULL && tempV0Last->IsB() && tempV0Last != Vp)
				{
					if (tempV0 == tempV0LastLast)
					{
						break;
					}
				}
				tempV0LastLast = tempV0Last;
				tempV0Last = tempV0;
			} while (1);
			Vp = tempVLast;
			tempVLast = tempV0Last;
			num++;
			if (num>1000)
			{
				break;
			}
		}

		teethBordersInOrder.push_back(flag2);
	}
	
	if (1)
	{
		mesh2->updateDataMask(MeshModel::MM_VERTFACETOPO);
		mesh2->updateDataMask(MeshModel::MM_FACEFACETOPO);

		teethBordersInOrder2.clear();
		teethBordersInOrder2.push_back(flagH1);

		int findtheFirst = 0;
		CVertexO* V1= NULL;

		CVertexO *vi = &(mesh2->cm.vert[flagH1]);
		vcg::face::JumpingPos<CFaceO> pos(vi->VFp(),vi);
		CVertexO* firstV = pos.VFlip();
		CVertexO* tempV= NULL;
		CVertexO* tempVLast= NULL;
		CVertexO* tempVLastLast= NULL;
		Point3f p1,p2;
		float dotproduct = 1000;
		do 
		{
			pos.NextE();
			tempV = pos.VFlip();
			// 			cout<<"vi:"<< vi->P().X()<< vi->P().Y()<< vi->P().Z()<<endl;
			// 			cout<<"firstV:"<< firstV->P().X()<< firstV->P().Y()<< firstV->P().Z()<<endl;
			// 			cout<<"tempV:"<< tempV->P().X()<< tempV->P().Y()<< tempV->P().Z()<<endl;
			if (tempVLast != NULL && tempVLastLast != NULL && tempVLast->IsB()){
				Point3f ptemp = tempVLast->P()-vi->P();
				if (ptemp*(zcoord)>dotproduct)//沿Y轴负方向(待改进，比较大的那一个)
				{
					if (tempV == tempVLastLast)
					{
						break;
						//cout<<"slotSortBorder2"<<endl;
					}
				}else{
					dotproduct = ptemp*(zcoord);
				}
			}
			tempVLastLast = tempVLast;
			tempVLast = tempV;
		} while (1);

		cout<<"slotSortBorder3"<<endl;

		CVertexO* Vp= vi;

		int num = 1;
		while (tempVLast != &mesh2->cm.vert[flagH2])//到另一个边界点
		{
			teethBordersInOrder2.push_back(tempVLast-&*(mesh2->cm.vert.begin()));
			CVertexO *vi0 = tempVLast;
			vcg::face::JumpingPos<CFaceO> pos0(vi0->VFp(),vi0);
			CVertexO* firstV0 = pos0.VFlip();
			CVertexO* tempV0= NULL;//访问V的一环邻域
			CVertexO* tempV0Last= NULL;
			CVertexO* tempV0LastLast= NULL;
			do 
			{

				pos0.NextE();
				tempV0 = pos0.VFlip();
				if (tempV0Last != NULL && tempV0LastLast != NULL && tempV0Last->IsB() && tempV0Last != Vp)
				{
					if (tempV0 == tempV0LastLast)
					{
						break;
					}
				}
				tempV0LastLast = tempV0Last;
				tempV0Last = tempV0;
			} while (1);
			Vp = tempVLast;
			tempVLast = tempV0Last;
			num++;
			if (num>1000)
			{
				break;
			}
		}

		teethBordersInOrder2.push_back(flagH2);
	}
	
	cout<<"slotSortBorder4"<<endl;

	if(CurrentMode() != TranMode::EXTENDB)
		SetMode(TranMode::EXTENDB);
	gla->update();
}

void EditRepairPlugin::slotLinkPoints()
{
	CMeshO::VertexIterator vi;
	CMeshO::FaceIterator fi;
	std::vector<CVertexO> vv_mm;
	std::vector<CFaceO> vf_mm;

	for (vi=mesh->cm.vert.begin();vi!=mesh->cm.vert.end();vi++){
		vv_mm.push_back(*vi);
	}
	for (vi=mesh2->cm.vert.begin();vi!=mesh2->cm.vert.end();vi++){
		vv_mm.push_back(*vi);
	}
	for (int i = 0; i<mesh->cm.face.size();i++){
		vf_mm.push_back(mesh->cm.face[i]);
	}
	for (int i = 0; i<mesh2->cm.face.size();i++){
		vf_mm.push_back(mesh2->cm.face[i]);
	}
	repaired_mesh = md->addNewMesh("repaired_mesh","",false);

	repaired_mesh->cm.Clear();
	int VN,FN;
	VN=vv_mm.size();
	FN=vf_mm.size()+teethBordersInOrder.size()-1 + teethBordersInOrder2.size()-1;

	//FN=vf_mm.size()+teethBordersInOrder.size()+rootBordersInOrder.size();
	vcg::tri::Allocator<CMeshO>::AddVertices(repaired_mesh->cm,VN);
	vcg::tri::Allocator<CMeshO>::AddFaces(repaired_mesh->cm,FN);
	CMeshO::VertexPointer *ivp = new CMeshO::VertexPointer[VN];
	int cnt=0;
	vi= repaired_mesh->cm.vert.begin();
	std::vector<CVertexO>::const_iterator itv;
	for(itv=vv_mm.begin();itv!=vv_mm.end();itv++){
		ivp[cnt]=&*vi;
		(*vi).P()= (*itv).P();
		++vi;
		++cnt;
	}
	fi=repaired_mesh->cm.face.begin();
	std::vector<CFaceO>::const_iterator itf;
	for(int i = 0; i<mesh->cm.face.size(); i++){
		for (int k = 0 ;k<3;k++)
		{
			int vInd = vf_mm[i].V(k)-&*(mesh->cm.vert.begin());

			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}
	int f1 = mesh->cm.face.size();
	int v1 = mesh->cm.vert.size();
	for(int i = f1; i<vf_mm.size(); i++){
		for (int k = 0 ;k<3;k++)
		{
			int vInd = vf_mm[i].V(k)-&*(mesh2->cm.vert.begin())+v1;
			(*fi).V(k)= ivp[vInd];
		}
		++fi;
	}
	int times = vf_mm.size();
	int remain;
	for (int i = 0; i < teethBordersInOrder.size()-1; i++)
	{
		int id1 = i;
		int id2 = i+1;
		float dis1 = Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh2->cm.vert[teethBordersInOrder2[0]].P());
		float dis2 = Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh2->cm.vert[teethBordersInOrder2[0]].P());
		int r1 = 0;
		int r2 = 0;
		for (int j = 1; j< teethBordersInOrder2.size()-1; j++)
		{
			if (Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh2->cm.vert[teethBordersInOrder2[j]].P())<dis1)
			{
				dis1 = Distance(mesh->cm.vert[teethBordersInOrder[id1]].P(),mesh2->cm.vert[teethBordersInOrder2[j]].P());
				r1 = j;
			}
			if (Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh2->cm.vert[teethBordersInOrder2[j]].P())<dis2)
			{
				dis2 = Distance(mesh->cm.vert[teethBordersInOrder[id2]].P(),mesh2->cm.vert[teethBordersInOrder2[j]].P());
				r2 = j;
			}
		}

		int vInd;
		vInd = teethBordersInOrder[id1];
		(*fi).V(0)= ivp[vInd];
		vInd = teethBordersInOrder[id2];
		(*fi).V(1)= ivp[vInd];
		vInd = teethBordersInOrder2[r1]+v1;
		(*fi).V(2)= ivp[vInd];
		if ((((*fi).V(1)->P()-(*fi).V(0)->P())^((*fi).V(2)->P()-(*fi).V(0)->P()))*((*fi).V(0)->P()-barycentre) < 0)
		{
			vInd = teethBordersInOrder[id1];
			(*fi).V(0)= ivp[vInd];
			vInd = teethBordersInOrder2[r1]+v1;		
			(*fi).V(1)= ivp[vInd];
			vInd = teethBordersInOrder[id2];			;
			(*fi).V(2)= ivp[vInd];
		}
		++fi;
		times++;
		remain = r2;

		if (r1 < r2)
		{
			for (int j = r1; j<r2; j++)
			{
				vInd = teethBordersInOrder2[j]+v1;
				(*fi).V(0)= ivp[vInd];
				vInd = teethBordersInOrder[id2];
				(*fi).V(1)= ivp[vInd];
				vInd = teethBordersInOrder2[j+1]+v1;
				(*fi).V(2)= ivp[vInd];
				if ((((*fi).V(1)->P()-(*fi).V(0)->P())^((*fi).V(2)->P()-(*fi).V(0)->P()))*((*fi).V(0)->P()-barycentre) < 0)
				{
					vInd = teethBordersInOrder2[j]+v1;
					(*fi).V(0)= ivp[vInd];
					vInd = teethBordersInOrder2[j+1]+v1;
					(*fi).V(1)= ivp[vInd];
					vInd = teethBordersInOrder[id2];
					(*fi).V(2)= ivp[vInd];
				}
				++fi;
				times++;
			}
			remain = r2;
		}else if (r1 > r2)
		{
			for (int j = r1; j<r2+teethBordersInOrder2.size(); j++)
			{
				vInd = teethBordersInOrder2[j]+v1;
				(*fi).V(0)= ivp[vInd];
				vInd = teethBordersInOrder[id2];
				(*fi).V(1)= ivp[vInd];
				vInd = teethBordersInOrder2[j+1]+v1;
				(*fi).V(2)= ivp[vInd];
				if ((((*fi).V(1)->P()-(*fi).V(0)->P())^((*fi).V(2)->P()-(*fi).V(0)->P()))*((*fi).V(0)->P()-barycentre) < 0)
				{
					vInd = teethBordersInOrder2[j]+v1;
					(*fi).V(0)= ivp[vInd];
					vInd = teethBordersInOrder2[j+1]+v1;
					(*fi).V(1)= ivp[vInd];
					vInd = teethBordersInOrder[id2];
					(*fi).V(2)= ivp[vInd];
				}
				++fi;
				times++;
			}
			remain = r1;
		}
	}
	while (remain < teethBordersInOrder2.size() - 1 )
	{
		int vInd = teethBordersInOrder2[remain]+v1;
		(*fi).V(0)= ivp[vInd];
		vInd = teethBordersInOrder[teethBordersInOrder.size() - 1];
		(*fi).V(1)= ivp[vInd];
		vInd = teethBordersInOrder2[remain+1]+v1;
		(*fi).V(2)= ivp[vInd];
		if ((((*fi).V(1)->P()-(*fi).V(0)->P())^((*fi).V(2)->P()-(*fi).V(0)->P()))*((*fi).V(0)->P()-barycentre) < 0)
		{
			vInd = teethBordersInOrder2[remain]+v1;
			(*fi).V(0)= ivp[vInd];
			vInd = teethBordersInOrder2[remain+1]+v1;
			(*fi).V(1)= ivp[vInd];
			vInd = teethBordersInOrder[teethBordersInOrder.size() - 1];
			(*fi).V(2)= ivp[vInd];
		}
		++fi;
		remain++; 
	}

	cout<<"Success!"<<times<<" "<<FN<<endl;

	vcg::tri::UpdateBounding<CMeshO>::Box(repaired_mesh->cm);
	vcg::tri::UpdateNormals<CMeshO>::PerVertexNormalizedPerFaceNormalized(repaired_mesh->cm);
	if(CurrentMode() != TranMode::LINKP)
		SetMode(TranMode::LINKP);
	gla->update();
}

void EditRepairPlugin::slotImportInsectionP()
{
	if(CurrentMode() != TranMode::LINKP)
		SetMode(TranMode::LINKP);
	gla->update();
}