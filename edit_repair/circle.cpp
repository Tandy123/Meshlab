#include "circle.h"
#include <vcg/simplex/edge/pos.h>
#include <vcg/space/intersection2.h>

vcg::Point3f Circle::RotationA(vcg::Point3f &srcv, double angle, vcg::Point3f &axis){
	double lenaxis = sqrt(axis.X()*axis.X()+axis.Y()*axis.Y()+axis.Z()*axis.Z());
	axis /= lenaxis;
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

	Eigen::Vector3f vf0(srcv.X(),srcv.Y(),srcv.Z());
	Eigen::MatrixXf mf(3,1);
	mf =matrix.transpose()*vf0;	
	vcg::Point3f dstv1 = vcg::Point3f(mf(0,0),mf(1,0),mf(2,0));
	return dstv1;
}

Circle::Circle(MeshModel* _m, vcg::Point3f _barycentrecoord, vcg::Point3f _xcoord, vcg::Point3f _ycoord, vcg::Point3f _zcoord, double _offset) {
	m = _m;
	barycentrecoord = _barycentrecoord;
	xcoord = _xcoord;
	ycoord = _ycoord;
	zcoord = _zcoord;
	offset = _offset;
	barycenter = barycentrecoord+_zcoord*offset;
	cutn = 10;
}

Circle::Circle(vcg::Point3f _barycentrecoord, vcg::Point3f _xcoord, vcg::Point3f _ycoord, vcg::Point3f _zcoord,double _offset,double r){
	barycentrecoord = _barycentrecoord;
	xcoord = _xcoord;
	ycoord = _ycoord;
	zcoord = _zcoord;
	offset = _offset;
	barycenter = barycentrecoord+_zcoord*offset;
	radius = r;
	circleType = BOTTOM;
	cutn = 10;
}

// void Circle::Classify(){
// 	outlines.clear();
// 	std::vector<vcg::Point3f> outline;
// 	vcg::tri::UpdateFlags<CMeshO>::EdgeClearV(m->cm);
// 	vcg::tri::UpdateTopology<CMeshO>::EdgeEdge(m->cm);
// 	int nv=0;
// 	for(size_t i=0;i<m->cm.edge.size();i++) if(!m->cm.edge[i].IsD())
// 	{
// 		if (!m->cm.edge[i].IsV())
// 		{
// 			vcg::edge::Pos<CMeshO::EdgeType> startE(&m->cm.edge[i],0);
// 			vcg::edge::Pos<CMeshO::EdgeType> curE=startE;
// 			vcg::edge::Pos<CMeshO::EdgeType> pre_curE;
// 			do
// 			{	
// 				curE.E()->SetV();
// 				pre_curE = curE;
// 				curE.NextE();
// 				if (pre_curE.E() == curE.E())
// 				{
// 					outline.push_back(pre_curE.V()->P());
// 					nv++;
// 				}
// 			}
// 			while(curE != startE);
// 			if (outline.size() > 0)
// 			{
// 				outlines.push_back(outline);
// 				outline.clear();
// 			}
// 		}
// 	}
// 	if (outlines.size() == 0)
// 	{
// 		circleType = FULL;
// 	}else if (outlines.size() == 1)
// 	{
// 		double cosZ = (outlines[0][0]-barycentrecoord) * zcoord;
// 		if (cosZ > 0)
// 		{
// 			vcg::Point3f v = (outlines[0][0] + outlines[0][1])/2-barycentrecoord;
// 			double cosA = v.Normalize() * ycoord;
// 			if (cosA >= 0)
// 			{
// 				if (cosA == 0)
// 				{
// 					qDebug() << "WARNING! Classify:cosA = 0" << "\n";
// 				}
// 				circleType = RIGHTGAP;
// 			}else{
// 				circleType = LEFTGAP;
// 			}
// 		}else{
// 			vcg::Point3f v = (outlines[0][0] + outlines[0][1])/2-barycentrecoord;
// 			double cosA = v.Normalize() * xcoord;
// 			if (cosA >= 0)
// 			{
// 				if (cosA == 0)
// 				{
// 					qDebug() << "WARNING! Classify:cosA = 0" << "\n";
// 				}
// 				circleType = ONLYFRONT;
// 			}else{
// 				circleType = ONLYBACK;
// 			}
// 		}
// 		
// 	}else if (outlines.size() == 2)
// 	{
// 		circleType = TWOGAP;
// 	}else{
// 		qDebug() << "ERROR! Classify:can not find type" << "\n";
// 	}
// }

void Circle::SortBorderPoints(){
	switch(circleType){
		case FULL:
			break;
		case LEFTGAP:
			{
				vcg::Point3f v = outlines[0][0] - outlines[0][1];
				double cosA = v * xcoord;
				if (cosA >= 0)
				{
					if (cosA == 0)
					{
						qDebug() << "WARNING! SortBorderPoints:cosA1 = 0" << "\n";
					}
				}else{
					vcg::Point3f ptemp = outlines[0][0];
					outlines[0][0] = outlines[0][1];
					outlines[0][1] = ptemp;
				}
				break;
			}
		case RIGHTGAP:
			{
				vcg::Point3f v = outlines[0][0] - outlines[0][1];
				double cosA = v * xcoord;
				if (cosA >= 0)
				{
					if (cosA == 0)
					{
						qDebug() << "WARNING! SortBorderPoints:cosA1 = 0" << "\n";
					}
				}else{
					vcg::Point3f ptemp = outlines[0][0];
					outlines[0][0] = outlines[0][1];
					outlines[0][1] = ptemp;
				}
				break;
			}
		case ONLYFRONT:
			{
				vcg::Point3f v = outlines[0][0] - outlines[0][1];
				double cosA = v * ycoord;
				if (cosA >= 0)
				{
					if (cosA == 0)
					{
						qDebug() << "WARNING! SortBorderPoints:cosA1 = 0" << "\n";
					}
				}else{
					vcg::Point3f ptemp = outlines[0][0];
					outlines[0][0] = outlines[0][1];
					outlines[0][1] = ptemp;
				}
				break;
			}
		case ONLYBACK:
			{
				vcg::Point3f v = outlines[0][1] - outlines[0][0];
				double cosA = v * ycoord;
				if (cosA >= 0)
				{
					if (cosA == 0)
					{
						qDebug() << "WARNING! SortBorderPoints:cosA1 = 0" << "\n";
					}
				}else{
					vcg::Point3f ptemp = outlines[0][0];
					outlines[0][0] = outlines[0][1];
					outlines[0][1] = ptemp;
				}
				break;
			}
		case TWOGAP:
			{
				vcg::Point3f p1 = (outlines[0][0] + outlines[0][1])/2.0;
				vcg::Point3f p2 = (outlines[1][0] + outlines[1][1])/2.0;
				vcg::Point3f v = p1 - p2;
				double cosA = v * xcoord;
				if (cosA >= 0)
				{
					if (cosA == 0)
					{
						qDebug() << "WARNING! SortBorderPoints:cosA = 0" << "\n";
					}
				}else{
					vcg::Point3f ptemp = outlines[0][0];
					outlines[0][0] = outlines[1][0];
					outlines[1][0] = ptemp;
					ptemp = outlines[0][1];
					outlines[0][1] = outlines[1][1];
					outlines[1][1] = ptemp;
				}
				v = outlines[0][0] - outlines[0][1];
				cosA = v * ycoord;
				if (cosA >= 0)
				{
					if (cosA == 0)
					{
						qDebug() << "WARNING! SortBorderPoints:cosA1 = 0" << "\n";
					}
				}else{
					vcg::Point3f ptemp = outlines[0][0];
					outlines[0][0] = outlines[0][1];
					outlines[0][1] = ptemp;
				}

				v = outlines[1][1] - outlines[1][0];
				cosA = v * ycoord;
				if (cosA >= 0)
				{
					if (cosA == 0)
					{
						qDebug() << "WARNING! SortBorderPoints:cosA2 = 0" << "\n";
					}
				}else{
					vcg::Point3f ptemp = outlines[1][1];
					outlines[1][1] = outlines[1][0];
					outlines[1][0] = ptemp;
				}
				break;
			}	
		default:
			break;
	}
}

std::vector< std::vector<vcg::Point3f> > Circle::GetOutLines(){
	return outlines;
}

std::vector<vcg::Point3f> Circle::GetPatchPoints(){
	return patchPoints;
}

std::vector<vcg::Point3f> Circle::GetPatchPointsR(){
	return patchPointsRight;
}

std::vector<vcg::Point3f> Circle::GetPatchPointsL(){
	return patchPointsLeft;
}

std::vector<vcg::Point3f> Circle::GetPatchPointsB(){
	return patchPointsBottom;
}


void Circle::LinkBorderPoints(){
	double PI = 3.1415926;
	patchPoints.clear();
	patchPointsLeft.clear();
	patchPointsRight.clear();
	patchPointsBottom.clear();
	switch(circleType){
		case FULL:
			break;
		case LEFTGAP:
			{
				double angle;
				vcg::Point3f v1 = outlines[0][0] - barycenter;
				vcg::Point3f v2 = outlines[0][1] - barycenter;
				double r1 = Distance(outlines[0][0],barycenter);
				double r2 = Distance(outlines[0][1],barycenter);
				v1.Normalize();
				v2.Normalize();
				angle = acos(v1*v2);
				double angleStep = angle / cutn;
				double radiusStep = (r2 - r1)/cutn;
				//vcg::Point3f rotation_axis = v1^v2;
				for (int i = 0; i < cutn +1; i++)
				{
					double deltAngle = angleStep * i;
					vcg::Point3f result;
					result = RotationA(v1, deltAngle, -zcoord);
					result = barycenter + result * (r1 + radiusStep * i);
					patchPoints.push_back(result);
					patchPointsLeft.push_back(result);
				}
				break;
			}
		case RIGHTGAP:
			{
				double angle;
				vcg::Point3f v1 = outlines[0][0] - barycenter;
				vcg::Point3f v2 = outlines[0][1] - barycenter;
				double r1 = Distance(outlines[0][0],barycenter);
				double r2 = Distance(outlines[0][1],barycenter);
				v1.Normalize();
				v2.Normalize();
				angle = acos(v1*v2);
				double angleStep = angle / cutn;
				double radiusStep = (r2 - r1)/cutn;
				//vcg::Point3f rotation_axis = v1^v2;
				for (int i = 0; i < cutn +1; i++)
				{
					double deltAngle = angleStep * i;
					vcg::Point3f result;
					result = RotationA(v1, deltAngle, zcoord);
					result = barycenter + result * (r1 + radiusStep * i);
					patchPoints.push_back(result);
					patchPointsRight.push_back(result);
				}
				break;
			}
		case ONLYFRONT:
			{
				double angle;
				vcg::Point3f v1 = outlines[0][0] - barycenter;
				vcg::Point3f v2 = outlines[0][1] - barycenter;
				double r1 = Distance(outlines[0][0],barycenter);
				double r2 = Distance(outlines[0][1],barycenter);
				v1.Normalize();
				v2.Normalize();
				angle = 2.0*PI - acos(v1*v2);
				double angleStep = angle / (2*cutn);
				double radiusStep = (r2 - r1)/(2*cutn);
				//vcg::Point3f rotation_axis = v2 ^ v1;
				for (int i = 0; i < cutn +1; i++)
				{
					double deltAngle = angleStep * i;
					vcg::Point3f result;
					result = RotationA(v1, deltAngle, zcoord);
					result = barycenter + result * (r1 + radiusStep * i);
					patchPoints.push_back(result);
					patchPointsRight.push_back(result);
				}
				for (int i = 2*cutn; i > cutn - 1; i--)
				{
					double deltAngle = angleStep * i;
					vcg::Point3f result;
					result = RotationA(v1, deltAngle, zcoord);
					result = barycenter + result * (r1 + radiusStep * i);
					patchPoints.push_back(result);
					patchPointsLeft.push_back(result);
				}
				break;
			}
		case ONLYBACK:
			{
				double angle;
				vcg::Point3f v1 = outlines[0][0] - barycenter;
				vcg::Point3f v2 = outlines[0][1] - barycenter;
				double r1 = Distance(outlines[0][0],barycenter);
				double r2 = Distance(outlines[0][1],barycenter);
				v1.Normalize();
				v2.Normalize();
				angle = 2 * PI - acos(v1*v2);
				double angleStep = angle / (2*cutn);
				double radiusStep = (r2 - r1)/(2*cutn);
				//vcg::Point3f rotation_axis = v2^v1;
				for (int i = cutn; i > -1; i--)
				{
					double deltAngle = angleStep * i;
					vcg::Point3f result;
					result = RotationA(v1, deltAngle, zcoord);
					result = barycenter + result * (r1 + radiusStep * i);
					patchPoints.push_back(result);
					patchPointsLeft.push_back(result);
				}
				for (int i = cutn; i < 2*cutn+1; i++)
				{
					double deltAngle = angleStep * i;
					vcg::Point3f result;
					result = RotationA(v1, deltAngle, zcoord);
					result = barycenter + result * (r1 + radiusStep * i);
					patchPoints.push_back(result);
					patchPointsRight.push_back(result);
				}
				break;
			}
		case TWOGAP:
			{
				double angle;
				vcg::Point3f v1 = outlines[0][0] - barycenter;
				vcg::Point3f v2 = outlines[1][1] - barycenter;
				double r1 = Distance(outlines[0][0],barycenter);
				double r2 = Distance(outlines[1][1],barycenter);
				v1.Normalize();
				v2.Normalize();
				angle = acos(v1*v2);
				double angleStep = angle / cutn;
				double radiusStep = (r2 - r1)/cutn;
				vcg::Point3f rotation_axis = v1^v2;
				for (int i = 0; i < cutn +1; i++)
				{
					double deltAngle = angleStep * i;
					vcg::Point3f result;
					result = RotationA(v1, deltAngle, zcoord);
					result = barycenter + result * (r1 + radiusStep * i);
					//result = result * (r1 + radiusStep * i);
					patchPoints.push_back(result);
					patchPointsRight.push_back(result);
				}

				v1 = outlines[0][1] - barycenter;
				v2 = outlines[1][0] - barycenter;
				r1 = Distance(outlines[0][1],barycenter);
				r2 = Distance(outlines[1][0],barycenter);
				v1.Normalize();
				v2.Normalize();
				angle = acos(v1*v2);
				angleStep = angle / cutn;
				radiusStep = (r2 - r1)/cutn;
				rotation_axis = v1^v2;
				for (int i = 0; i < cutn +1; i++)
				{
					double deltAngle = angleStep * i;
					vcg::Point3f result;
					result = RotationA(v1, deltAngle, -zcoord);
					result = barycenter + result * (r1 + radiusStep * i);
					patchPoints.push_back(result);
					patchPointsLeft.push_back(result);
				}
				break;
			}
		case BOTTOM:
			{
				double angle;
				vcg::Point3f v = xcoord;
				double r = radius;
				v.Normalize();
				angle = 2.0*PI;
				double angleStep = angle / (2*cutn);
				//vcg::Point3f rotation_axis = v2 ^ v1;
				for (int i = 0; i < 2*cutn+1; i++)
				{
					double deltAngle = angleStep * i;
					vcg::Point3f result;
					result = RotationA(v, deltAngle, zcoord);
					result = barycenter + result * r;
					patchPoints.push_back(result);
					patchPointsBottom.push_back(result);
/*					patchPointsRight.push_back(result);*/
				}
// 				for (int i = 2*cutn; i > cutn - 1; i--)
// 				{
// 					double deltAngle = angleStep * i;
// 					vcg::Point3f result;
// 					result = RotationA(v1, deltAngle, zcoord);
// 					result = barycenter + result * (r1 + radiusStep * i);
// 					patchPoints.push_back(result);
// 					patchPointsLeft.push_back(result);
// 				}
				break;
			}
		default:
			break;
	}
}