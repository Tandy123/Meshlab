#ifndef TEETH_TREE_H
#define TEETH_TREE_H

#include <common/interfaces.h>

class ToothNode
{
public:
	ToothNode(MeshModel *_m, int _id,std::string _nodename)
	{
		m=_m;
		id = _id;
		nodename = _nodename;
		qnodename = QString::fromLocal8Bit(_nodename.c_str());
		glued=false;
		toothaxisvex.clear();
		isgum = false;
		isattachment = false;
		nodetr.SetIdentity();
		origntr.SetIdentity();
		current_save_path.SetIdentity();
	}
	ToothNode() { m=0;id=-1;nodename="0";qnodename = "0";}
	bool glued;
	int id;
	std::string nodename;
	QString qnodename;
	//vcg::Matrix33f pca;
	MeshModel *m;
	vcg::Matrix44f &tr() {return m->cm.Tr;}
	const vcg::Box3f &bbox() const {return m->cm.bbox;}
	std::vector<vcg::Point3f> toothaxisvex;
//	vcg::TeethPath path;

	vcg::Point3f translation,rotation;
	void UpdateOriginPCA();
	void UpdatefeaturePoint(vcg::Matrix44f r);
	void UpdateTr(){m->cm.Tr = current_save_path;}
	bool isgum;
	bool isattachment;
	int move_state[60];
	std::vector<vcg::Matrix44f> move_state_matrix;
	void CalmatrixByTranData();
	bool iseditingkeyframe;//标记当前牙齿是否修改过。
	std::vector<int> edit_state;
	vcg::Matrix44f nodetr;
	vcg::Matrix44f origntr;
	vcg::Matrix44f current_save_path;//保存当前牙齿已经保存的path
};

class TeethTree
{
public:
	TeethTree(){};

	QList<ToothNode *> nodeList;
	vcg::CallBackPos * cb;

	MeshModel *MM(unsigned int i) {return nodeList.value(i)->m;}

	void clear()
	{
		foreach(ToothNode *mp, nodeList) 
			delete mp;
		nodeList.clear();
	}

	void resetID();

	ToothNode *find(int id)
	{
		foreach(ToothNode *mp, nodeList) 
			if(mp->id==id) return mp;
		assert("You are trying to find an unexistent mesh"==0);
		return 0;
	}

	ToothNode *find(MeshModel *m)
	{
		foreach(ToothNode *mp, nodeList) 
			if(mp->m==m) return mp;
		assert("You are trying to find an unexistent mesh"==0);
		return 0;
	}
	int gluedNum();

	inline vcg::Box3f bbox() {
		vcg::Box3f FullBBox;
		foreach(ToothNode *mp, nodeList) 
			FullBBox.Add(vcg::Matrix44f::Construct(mp->tr()),mp->bbox());
		return FullBBox;
	}

	inline vcg::Box3f gluedBBox() {
		vcg::Box3f FullBBox;
		foreach(ToothNode *mp, nodeList) 
			if(mp->glued)
				FullBBox.Add(vcg::Matrix44f::Construct(mp->tr()),mp->bbox());
		return FullBBox;
	}
};

class IdMesh{//自定义数据结果，用于存储网格和序号
public:
	int id;
	MeshModel *mm;
	std::vector<vcg::Point3f> local_coordinate;
	vcg::Point3f _CMPCA[3];
	vcg::Point3f _repairax[3];
	vcg::Point3f _barycentric_coord;
	vcg::Point3f _initbarycentric_coord;
	vcg::Point3f _rootpoint;
	bool _isrepairing_axis;
	int  _rotationaxis;//记录绕哪个轴转动，用于绘制调整轴时的坐标
	float _height;
	bool haveBottom;
	int bottomFlag;
	std::vector<CVertexO> vHole;
	int deleteFn;
};

#endif