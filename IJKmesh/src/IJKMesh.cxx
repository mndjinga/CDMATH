/*
 * IJKmesh.cxx
 *
 *  Created on: 24 March 2019
 *      Authors: CDMATH
 */

#include "IJKMesh.hxx"
#include "IJKNode.hxx"
#include "IJKCell.hxx"
#include "IJKFace.hxx"

#include "MEDFileMesh.hxx"
#include "MEDLoader.hxx"
#include "MEDCouplingIMesh.hxx"

#include "CdmathException.hxx"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iterator>
#include <algorithm> 

using namespace MEDCoupling;
using namespace std;

//----------------------------------------------------------------------
Mesh::Mesh( void )
//----------------------------------------------------------------------
{
	_mesh=NULL;
	_spaceDim = 0 ;
	_meshDim  = 0 ;
	_numberOfNodes = 0;
	_numberOfFaces = 0;
	_numberOfCells = 0;
	_numberOfEdges = 0;
	_xMin=0.;
	_xMax=0.;
	_yMin=0.;
	_yMax=0.;
	_zMin=0.;
	_zMax=0.;
    _nxyz.resize(0);
    _dxyz.resize(0.);
	_faceGroupNames.resize(0);
	_faceGroups.resize(0);
	_nodeGroupNames.resize(0);
	_nodeGroups.resize(0);
    _indexFacePeriodicSet=false;
    _name="";
}

//----------------------------------------------------------------------
Mesh::~Mesh( void )
//----------------------------------------------------------------------
{
}

std::string 
Mesh::getName( void ) const
{
    return _name;
}

Mesh::Mesh( const MEDCoupling::MEDCouplingIMesh* mesh )
{
	_spaceDim=mesh->getSpaceDimension();
	_meshDim=mesh->getMeshDimension();
	_dxyz=mesh->getDXYZ();
	_nxyz=mesh->getCellGridStructure();
	double* Box0=new double[2*_spaceDim];
	mesh->getBoundingBox(Box0);
    _name=mesh->getName();
    _indexFacePeriodicSet=false;
    
	_xMin=Box0[0];
	_xMax=Box0[1];
	if (_spaceDim>=2)
	{
		_yMin=Box0[2];
		_yMax=Box0[3];
	}
	if (_spaceDim>=3)
	{
		_zMin=Box0[4];
		_zMax=Box0[5];
	}

	double *originPtr = new double[_spaceDim];
	double *dxyzPtr = new double[_spaceDim];
	int *nodeStrctPtr = new int[_spaceDim];

	for(int i=0;i<_spaceDim;i++)
	{
		originPtr[i]=Box0[2*i];
		nodeStrctPtr[i]=_nxyz[i]+1;
		dxyzPtr[i]=_dxyz[i];
	}
	_mesh=MEDCouplingIMesh::New(_name,
			_spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+_spaceDim,
			originPtr,
			originPtr+_spaceDim,
			dxyzPtr,
			dxyzPtr+_spaceDim);
	delete [] originPtr;
	delete [] dxyzPtr;
	delete [] nodeStrctPtr;
	delete [] Box0 ;
	setMesh();
}

//----------------------------------------------------------------------
Mesh::Mesh( const IJKMesh& m )
//----------------------------------------------------------------------
{
	_spaceDim = m.getSpaceDimension() ;
	_meshDim = m.getMeshDimension() ;
    _name=m.getName();
    _xMax=m.getXMax();
    _yMin=m.getYMin();
    _yMax=m.getYMax();
    _zMin=m.getZMin();
    _zMax=m.getZMax();
    _nxyz = m.getCellGridStructure() ;
    _dxyz=m.getDXYZ();

	_numberOfNodes = m.getNumberOfNodes();
	_numberOfFaces = m.getNumberOfFaces();
	_numberOfCells = m.getNumberOfCells();
	_numberOfEdges = m.getNumberOfEdges();
	_faceGroupNames = m.getNameOfFaceGroups() ;
	_faceGroups = m.getFaceGroups() ;
	_nodeGroupNames = m.getNameOfNodeGroups() ;
	_nodeGroups = m.getNodeGroups() ;
    _indexFacePeriodicSet= m.isIndexFacePeriodicSet();
    if(_indexFacePeriodicSet)
        _indexFacePeriodicMap=m.getIndexFacePeriodic();
    
	MCAuto<MEDCouplingIMesh> m1=m.getMEDCouplingIMesh()->deepCopy();
	_mesh=m1;
}

//----------------------------------------------------------------------
Mesh::Mesh( const std::string filename, int meshLevel )
//----------------------------------------------------------------------
{
	readMeshMed(filename, meshLevel);
}

//----------------------------------------------------------------------
void
Mesh::readMeshMed( const std::string filename, const int meshLevel)
//----------------------------------------------------------------------
{
	MEDFileIMesh *m=MEDFileIMesh::New(filename.c_str());//reads the first mesh encountered in the file, otherwise call New (const char *fileName, const char *mName, int dt=-1, int it=-1)
	_mesh=m->getMeshAtLevel(meshLevel);
    _mesh->checkConsistencyLight();
	_mesh->setName(_mesh->getName());
	_meshDim=_mesh->getMeshDimension();
	_spaceDim=_mesh->getSpaceDimension();
    _name=_mesh->getName();
    _indexFacePeriodicSet=false;
    MEDCoupling::MEDCouplingIMesh* structuredMesh = dynamic_cast<MEDCoupling::MEDCouplingIMesh*> (_mesh.retn());
    if(structuredMesh)
    {
        _dxyz=structuredMesh->getDXYZ();
        _nxyz=structuredMesh->getCellGridStructure();
        double* Box0=new double[2*_spaceDim];
        structuredMesh->getBoundingBox(Box0);
    
        _xMin=Box0[0];
        _xMax=Box0[1];
        std::cout<<"nx= "<<_nxyz[0];
        if (_spaceDim>=2)
        {
            _yMin=Box0[2];
            _yMax=Box0[3];
            std::cout<<", "<<"ny= "<<_nxyz[1];
        }
        if (_spaceDim>=3)
        {
            _zMin=Box0[4];
            _zMax=Box0[5];
            std::cout<<", "<<"nz= "<<_nxyz[2];
        }
    }
    else
        throw CdmathException("Mesh::readMeshMed med file does not contain a structured MedcouplingIMesh mesh");
    
	MEDCouplingIMesh*  mu = setMesh();
	setGroups(m, mu);
	cout<<endl<< "Loaded file "<< filename<<endl;
    cout<<"Structured Mesh name= "<<m->getName()<<", mesh dim="<< _meshDim<< ", space dim="<< _spaceDim<< ", nb cells= "<<getNumberOfCells()<< ", nb nodes= "<<getNumberOfNodes()<<endl;

	m->decrRef();
	mu->decrRef();
}

void
Mesh::setGroupAtFaceByCoords(double x, double y, double z, double eps, std::string groupName)
{
	int nbFace=getNumberOfFaces();
	bool flag=false;
	for (int iface=0;iface<nbFace;iface++)
	{
		double FX=_faces[iface].x();
		double FY=_faces[iface].y();
		double FZ=_faces[iface].z();
		if (abs(FX-x)<eps && abs(FY-y)<eps && abs(FZ-z)<eps)
		{
			_faces[iface].setGroupName(groupName);
			IntTab nodesID= _faces[iface].getNodesId();
			int nbNodes = _faces[iface].getNumberOfNodes();
			for(int inode=0 ; inode<nbNodes ; inode++)
				_nodes[nodesID[inode]].setGroupName(groupName);

			flag=true;
		}
	}
	if (flag)
    {
		_faceGroupNames.push_back(groupName);
		_nodeGroupNames.push_back(groupName);
        //To do : update _faceGroups and _nodeGroups
    }
}

void
Mesh::setGroupAtPlan(double value, int direction, double eps, std::string groupName)
{
	int nbFace=getNumberOfFaces();
	bool flag=false;
	for (int iface=0;iface<nbFace;iface++)
	{
		double cord=_faces[iface].getBarryCenter()[direction];
		if (abs(cord-value)<eps)
		{
			_faces[iface].setGroupName(groupName);
			IntTab nodesID= _faces[iface].getNodesId();
			int nbNodes = _faces[iface].getNumberOfNodes();
			for(int inode=0 ; inode<nbNodes ; inode++)
                {
				_nodes[nodesID[inode]].setGroupName(groupName);
                }

			flag=true;
		}
	}
	if (flag)
    {
		_faceGroupNames.push_back(groupName);
		_nodeGroupNames.push_back(groupName);
        //To do : update _faceGroups, _nodeGroups
    }
}

std::map<int,int>
Mesh::getIndexFacePeriodic( void ) const
{
    return _indexFacePeriodicMap;
}

void
Mesh::setPeriodicFaces(bool check_groups, bool use_central_inversion)
{
    if(_indexFacePeriodicSet)
        return;
        
    double eps=1.E-10;

    for (int indexFace=0;indexFace<_boundaryFaceIds.size() ; indexFace++)
    {
        Face my_face=_faces[_boundaryFaceIds[indexFace]];
        int iface_perio=-1;
        if(_meshDim==1)
        {
            for (int iface=0;iface<_boundaryFaceIds.size() ; iface++)
                if(iface!=indexFace)
                {
                    iface_perio=_boundaryFaceIds[iface];
                    break;
                }
        }
        else if(_meshDim==2)
        {
            double x=my_face.x();
            double y=my_face.y();
            
            for (int iface=0;iface<_boundaryFaceIds.size() ; iface++)
            {
                Face face_i=_faces[_boundaryFaceIds[iface]];
                double xi=face_i.x();
                double yi=face_i.y();
                if (   (abs(y-yi)<eps || abs(x-xi)<eps )// Case of a square geometry
                    && ( !check_groups || my_face.getGroupName()!=face_i.getGroupName()) //In case groups need to be checked
                    && ( !use_central_inversion || abs(y+yi) + abs(x+xi)<eps ) // Case of a central inversion
                    && fabs(my_face.getMeasure() - face_i.getMeasure())<eps
                    && fabs(my_face.getXN()      + face_i.getXN())<eps
                    && fabs(my_face.getYN()      + face_i.getYN())<eps
                    && fabs(my_face.getZN()      + face_i.getZN())<eps )
                {
                    iface_perio=_boundaryFaceIds[iface];
                    break;
                }
            }
        }
        else if(_meshDim==3)
        {
            double x=my_face.x();
            double y=my_face.y();
            double z=my_face.z();
        
            for (int iface=0;iface<_boundaryFaceIds.size() ; iface++)
            {
                Face face_i=_faces[_boundaryFaceIds[iface]];
                double xi=face_i.x();
                double yi=face_i.y();
                double zi=face_i.z();
                if ( ((abs(y-yi)<eps && abs(x-xi)<eps) || (abs(x-xi)<eps && abs(z-zi)<eps) || (abs(y-yi)<eps && abs(z-zi)<eps))// Case of a cube geometry
                    && ( !check_groups || my_face.getGroupName()!=face_i.getGroupName()) //In case groups need to be checked
                    && ( !use_central_inversion || abs(y+yi) + abs(x+xi) + abs(z+zi)<eps )// Case of a central inversion
                    && fabs(my_face.getMeasure() - face_i.getMeasure())<eps
                    && fabs(my_face.getXN()      + face_i.getXN())<eps
                    && fabs(my_face.getYN()      + face_i.getYN())<eps
                    && fabs(my_face.getZN()      + face_i.getZN())<eps )
                {
                    iface_perio=_boundaryFaceIds[iface];
                    break;
                }
            }  
        }
        else
            throw CdmathException("Mesh::setPeriodicFaces: Mesh dimension should be 1, 2 or 3");
        
        if (iface_perio==-1)
            throw CdmathException("Mesh::setPeriodicFaces: periodic face not found, iface_perio==-1 " );
        else
            _indexFacePeriodicMap[_boundaryFaceIds[indexFace]]=iface_perio;
    }
    _indexFacePeriodicSet=true;    
}

int
Mesh::getIndexFacePeriodic(int indexFace, bool check_groups, bool use_central_inversion)
{
	if (!_faces[indexFace].isBorder())
        {
            cout<<"Pb with indexFace= "<<indexFace<<endl;
            throw CdmathException("Mesh::getIndexFacePeriodic: not a border face" );
        }
        
    if(!_indexFacePeriodicSet)
        setPeriodicFaces(check_groups, use_central_inversion);

    std::map<int,int>::const_iterator  it = _indexFacePeriodicMap.find(indexFace);
    if( it != _indexFacePeriodicMap.end() )
        return it->second;
    else
    {
        cout<<"Pb with indexFace= "<<indexFace<<endl;
        throw CdmathException("Mesh::getIndexFacePeriodic: not a periodic face" );
    }
}

bool
Mesh::isBorderNode(int nodeid) const
{
	return getNode(nodeid).isBorder();
}

bool
Mesh::isBorderFace(int faceid) const
{
	return getFace(faceid).isBorder();
}

bool
Mesh::isQuadrangular() const
{
	return _eltsTypes.size()==1 && _eltsTypes[0]==INTERP_KERNEL::NORM_QUAD4;
}
bool
Mesh::isHexahedral() const
{
	return _eltsTypes.size()==1 && _eltsTypes[0]==INTERP_KERNEL::NORM_HEXA8;
}
bool
Mesh::isStructured() const
{
	return true;
}

std::string 
Mesh::getElementTypes() const
{
    std::string result;    
    for(int i=0; i< _eltsTypes.size(); i++)
    {
        if( _eltsTypes[i]==INTERP_KERNEL::NORM_POINT1)
            result += "Points ";
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_SEG2)
            result += "Segments ";
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_TRI3)
            result += "Triangles ";
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_QUAD4)
            result += "Quadrangles ";
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_POLYGON)
            result += "Polygons ";
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_TETRA4)
            result += "Tetrahedra ";
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_PYRA5)
            result += "Pyramids ";
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_PENTA6)
            result += "Pentahedra ";
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_HEXA8)
            result += "Hexahedra ";
        else if( _eltsTypes[i]==INTERP_KERNEL::NORM_POLYHED)
            result += "Polyhedrons ";
        else
		{
			cout<< "Mesh " + _name + " contains an element of type " <<endl;
			cout<< _eltsTypes[i]<<endl;
			throw CdmathException("Mesh::getElementTypes : recognised cell med types are NORM_POINT1, NORM_SEG2, NORM_TRI3, NORM_QUAD4, NORM_TETRA4, NORM_PYRA5, NORM_PENTA6, NORM_HEXA8, NORM_POLYGON, NORM_POLYHED");
        }
    }
    return result;
}

void
Mesh::setGroups( const MEDFileIMesh* medmesh, MEDCouplingIMesh*  mu)
{
	//Searching for face groups
	vector<string> faceGroups=medmesh->getGroupsNames() ;

	for (unsigned int i=0;i<faceGroups.size();i++ )
	{
		string groupName=faceGroups[i];
		vector<int> nonEmptyGrp(medmesh->getGrpNonEmptyLevels(groupName));
		//We check if the group has a relative dimension equal to -1 
		//before call to the function getGroup(-1,groupName.c_str())
		vector<int>::iterator it = find(nonEmptyGrp.begin(), nonEmptyGrp.end(), -1);
		if (it != nonEmptyGrp.end())
		{
			cout<<"Boundary face group named "<< groupName << " found"<<endl;
			MEDCouplingIMesh *m=medmesh->getGroup(-1,groupName.c_str());
			_faceGroups.push_back(m);
			_faceGroupNames.push_back(groupName);
			DataArrayDouble *baryCell = m->computeCellCenterOfMass() ;
			const double *coorBary=baryCell->getConstPointer();

			int nbCellsSubMesh=m->getNumberOfCells();
			for (int ic(0), k(0); ic<nbCellsSubMesh; ic++, k+=_spaceDim)
			{
				vector<double> coorBaryXyz(3,0);
				for (int d=0; d<_spaceDim; d++)
					coorBaryXyz[d] = coorBary[k+d];
				Point p1(coorBaryXyz[0],coorBaryXyz[1],coorBaryXyz[2]) ;

				int flag=0;
				for (int iface=0;iface<_numberOfFaces;iface++ )
				{
					Point p2=_faces[iface].getBarryCenter();
					if(p1.distance(p2)<1.E-10)
					{
						_faces[iface].setGroupName(groupName);
						flag=1;
						break;
					}
				}
				if (flag==0)
					throw CdmathException("No face belonging to group " + groupName + " found");
			}
			baryCell->decrRef();
			//m->decrRef();
		}
	}

	//Searching for node groups
	vector<string> nodeGroups=medmesh->getGroupsOnSpecifiedLev(1) ;

	for (unsigned int i=0;i<nodeGroups.size();i++ )
	{
		string groupName=nodeGroups[i];
		DataArrayInt * nodeGroup=medmesh->getNodeGroupArr( groupName );
		const int *nodeids=nodeGroup->getConstPointer();

		if(nodeids!=NULL)
		{
			cout<<"Boundary node group named "<< groupName << " found"<<endl;

			_nodeGroups.push_back(nodeGroup);
			_nodeGroupNames.push_back(groupName);

			int nbNodesSubMesh=nodeGroup->getNumberOfTuples();//nodeGroup->getNbOfElems();

			DataArrayDouble *coo = mu->getCoords() ;
			const double *cood=coo->getConstPointer();

			for (int ic(0); ic<nbNodesSubMesh; ic++)
			{
				vector<double> coorP(3,0);
				for (int d=0; d<_spaceDim; d++)
					coorP[d] = cood[nodeids[ic]*_spaceDim+d];
				Point p1(coorP[0],coorP[1],coorP[2]) ;

				int flag=0;
				for (int inode=0;inode<_numberOfNodes;inode++ )
				{
					Point p2=_nodes[inode].getPoint();
					if(p1.distance(p2)<1.E-10)
					{
						_nodes[inode].setGroupName(groupName);
						flag=1;
						break;
					}
				}
				if (flag==0)
					throw CdmathException("No node belonging to group " + groupName + " found");
			}
		}
	}
}

//----------------------------------------------------------------------
MEDCouplingIMesh* 
Mesh::setMesh( void )
//----------------------------------------------------------------------
{

	// _cells, _nodes and _faces (with normal) initialization:
	
    return mu;
}

//----------------------------------------------------------------------
Mesh::Mesh( double xmin, double xmax, int nx, std::string meshName )
//----------------------------------------------------------------------
{
	if(nx<=0)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx) : nx <= 0");
	if(xmin>=xmax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx) : xmin >= xmax");

	double dx = (xmax - xmin)/nx ;

	_spaceDim = 1 ;
	_meshDim  = 1 ;
    _name=meshName;
    _indexFacePeriodicSet=false;

	_xMin=xmin;
	_xMax=xmax;
	_yMin=0.;
	_yMax=0.;
	_zMin=0.;
	_zMax=0.;

	_dxyz.resize(_spaceDim);
	_dxyz[0]=dx;
	_nxyz.resize(_spaceDim);
	_nxyz[0]=nx;

	double *originPtr = new double[_spaceDim];
	double *dxyzPtr = new double[_spaceDim];
	int *nodeStrctPtr = new int[_spaceDim];

	originPtr[0]=xmin;
	nodeStrctPtr[0]=nx+1;
	dxyzPtr[0]=dx;

	_mesh=MEDCouplingIMesh::New(meshName,
			_spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+_spaceDim,
			originPtr,
			originPtr+_spaceDim,
			dxyzPtr,
			dxyzPtr+_spaceDim);
	delete [] originPtr;
	delete [] dxyzPtr;
	delete [] nodeStrctPtr;

	_numberOfCells = _mesh->getNumberOfCells() ;

	_numberOfNodes = mu->getNumberOfNodes() ;

	_numberOfFaces = _numberOfNodes;

    _numberOfEdges = _numberOfCells;
    

		Point p(coorBary[id],0.0,0.0) ;
		Cell ci( nbVertices, nbVertices, lon[id], p ) ;
        double xn = (cood[nodeIdsOfCell[nbVertices-1]] - cood[nodeIdsOfCell[0]] > 0.0) ? -1.0 : 1.0;

        int nbFaces=tmpI[id+1]-tmpI[id];
        const int *work=tmp+tmpI[id];
		
        for( int el=0;el<nbFaces;el++ )
		{
			ci.addNormalVector(el,xn,0.0,0.0) ;
			ci.addFaceId(el,work[el]) ;
			xn = - xn;
		}

        
		Node vi( nbCells, nbFaces, nbNeighbourNodes, p ) ;
        for( int el=0;el<nbCells;el++ )
			vi.addCellId(el,workc[el]) ;
		for( int el=0;el<nbFaces;el++ )
			vi.addFaceId(el,id) ;
        for( int el=0;el<nbNeighbourNodes;el++ )
			vi.addNeighbourNodeId(el,workn[el]) ;


		int nbVertices=1;
		Face fi( nbVertices, nbCells, 1.0, p, 1., 0., 0. ) ;
        for( int el=0;el<nbVertices;el++ )
			fi.addNodeId(el,id) ;

		for( int el=0;el<nbCells;el++ )
			fi.addCellId(el,workc[el]) ;
}

//----------------------------------------------------------------------
Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, std::string meshName)
//----------------------------------------------------------------------
{
	if(nx<=0 || ny<=0)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : nx <= 0 or ny <= 0");
	if(xmin>=xmax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : xmin >= xmax");
	if(ymin>=ymax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny) : ymin >= ymax");

	_xMin=xmin;
	_xMax=xmax;
	_yMin=ymin;
	_yMax=ymax;
	_zMin=0.;
	_zMax=0.;


	double dx = (xmax - xmin)/nx ;
	double dy = (ymax - ymin)/ny ;

	_spaceDim = 2 ;
	_meshDim  = 2 ;
    _name=meshName;
    _indexFacePeriodicSet=false;
	_nxyz.resize(_spaceDim);
	_nxyz[0]=nx;
	_nxyz[1]=ny;

	_dxyz.resize(_spaceDim);
	_dxyz[0]=dx;
	_dxyz[1]=dy;

	double *originPtr = new double[_spaceDim];
	double *dxyzPtr = new double[_spaceDim];
	int *nodeStrctPtr = new int[_spaceDim];

	originPtr[0]=xmin;
	originPtr[1]=ymin;
	nodeStrctPtr[0]=nx+1;
	nodeStrctPtr[1]=ny+1;
	dxyzPtr[0]=dx;
	dxyzPtr[1]=dy;

	_mesh=MEDCouplingIMesh::New(meshName,
			_spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+_spaceDim,
			originPtr,
			originPtr+_spaceDim,
			dxyzPtr,
			dxyzPtr+_spaceDim);

	delete [] originPtr;
	delete [] dxyzPtr;
	delete [] nodeStrctPtr;
    
	setMesh();
}

//----------------------------------------------------------------------
Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz, std::string meshName)
//----------------------------------------------------------------------
{
	if(nx<=0 || ny<=0 || nz<=0)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : nx <= 0 or ny <= 0 or nz <= 0");
	if(xmin>=xmax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : xmin >= xmax");
	if(ymin>=ymax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : ymin >= ymax");
	if(zmin>=zmax)
		throw CdmathException("Mesh::Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz) : zmin >= zmax");

	_spaceDim = 3;
	_meshDim  = 3;
    _name=meshName;
    _indexFacePeriodicSet=false;
	_xMin=xmin;
	_xMax=xmax;
	_yMin=ymin;
	_yMax=ymax;
	_zMin=zmin;
	_zMax=zmax;

	double dx = (xmax - xmin)/nx ;
	double dy = (ymax - ymin)/ny ;
	double dz = (zmax - zmin)/nz ;

	_dxyz.resize(_spaceDim);
	_dxyz[0]=dx;
	_dxyz[1]=dy;
	_dxyz[2]=dz;

	_nxyz.resize(_spaceDim);
	_nxyz[0]=nx;
	_nxyz[1]=ny;
	_nxyz[2]=nz;

	double *originPtr = new double[_spaceDim];
	double *dxyzPtr = new double[_spaceDim];
	int *nodeStrctPtr = new int[_spaceDim];

	originPtr[0]=xmin;
	originPtr[1]=ymin;
	originPtr[2]=zmin;
	nodeStrctPtr[0]=nx+1;
	nodeStrctPtr[1]=ny+1;
	nodeStrctPtr[2]=nz+1;
	dxyzPtr[0]=dx;
	dxyzPtr[1]=dy;
	dxyzPtr[2]=dz;

	_mesh=MEDCouplingIMesh::New(meshName,
			_spaceDim,
			nodeStrctPtr,
			nodeStrctPtr+_spaceDim,
			originPtr,
			originPtr+_spaceDim,
			dxyzPtr,
			dxyzPtr+_spaceDim);

	delete [] originPtr;
	delete [] dxyzPtr;
	delete [] nodeStrctPtr;
    
	setMesh();
}

//----------------------------------------------------------------------
int
Mesh::getSpaceDimension( void )  const
//----------------------------------------------------------------------
{
	return _spaceDim ;
}

//----------------------------------------------------------------------
int
Mesh::getMeshDimension( void )  const
//----------------------------------------------------------------------
{
	return _meshDim ;
}

std::vector<double>
Mesh::getDXYZ() const
{
	return _dxyz;
}

std::vector<int>
Mesh::getCellGridStructure() const
{
	return _nxyz;
}

//----------------------------------------------------------------------
int
Mesh::getNx( void )  const
//----------------------------------------------------------------------
{
	return _nxyz[0];
}

//----------------------------------------------------------------------
int
Mesh::getNy( void )  const
//----------------------------------------------------------------------
{
	if(_spaceDim < 2)
		throw CdmathException("int double& Field::operator[ielem] : Ny is not defined in dimension < 2!");

	return _nxyz[1];
}

//----------------------------------------------------------------------
int
Mesh::getNz( void )  const
//----------------------------------------------------------------------
{
	if(_spaceDim < 3)
		throw CdmathException("int Mesh::getNz( void ) : Nz is not defined in dimension < 3!");

	return _nxyz[2];
}

//----------------------------------------------------------------------
double
Mesh::getXMin( void )  const
//----------------------------------------------------------------------
{        
	return _xMin ;
}

//----------------------------------------------------------------------
double
Mesh::getXMax( void )  const
//----------------------------------------------------------------------
{
	return _xMax ;
}

//----------------------------------------------------------------------
double
Mesh::getYMin( void )  const
//----------------------------------------------------------------------
{
	return _yMin ;
}

//----------------------------------------------------------------------
double
Mesh::getYMax( void )  const
//----------------------------------------------------------------------
{
	return _yMax ;
}

//----------------------------------------------------------------------
double
Mesh::getZMin( void )  const
//----------------------------------------------------------------------
{
	return _zMin ;
}

//----------------------------------------------------------------------
double
Mesh::getZMax( void )  const
//----------------------------------------------------------------------
{
	return _zMax ;
}

//----------------------------------------------------------------------
MCAuto<MEDCouplingIMesh>
Mesh::getMEDCouplingIMesh( void )  const
//----------------------------------------------------------------------
{
	return _mesh ;
}

//----------------------------------------------------------------------
int
Mesh::getNumberOfNodes ( void ) const
//----------------------------------------------------------------------
{
	return _numberOfNodes ;
}

//----------------------------------------------------------------------
int
Mesh::getNumberOfCells ( void ) const
//----------------------------------------------------------------------
{
	return _numberOfCells ;
}

//----------------------------------------------------------------------
int
Mesh::getNumberOfFaces ( void ) const
//----------------------------------------------------------------------
{
	return _numberOfFaces ;
}

//----------------------------------------------------------------------
int
Mesh::getNumberOfEdges ( void ) const
//----------------------------------------------------------------------
{
	return _numberOfEdges ;
}

//----------------------------------------------------------------------
Cell&
Mesh::getCell ( int i ) const
//----------------------------------------------------------------------
{
	return _cells[i] ;
}

//----------------------------------------------------------------------
Face&
Mesh::getFace ( int i ) const
//----------------------------------------------------------------------
{
	return _faces[i] ;
}

//----------------------------------------------------------------------
Node&
Mesh::getNode ( int i ) const
//----------------------------------------------------------------------
{
	return _nodes[i] ;
}

vector<string>
Mesh::getNameOfFaceGroups( void )  const
{
	return _faceGroupNames;
}

vector<MEDCoupling::MEDCouplingIMesh *>
Mesh::getFaceGroups( void )  const
{
	return _faceGroups;
}

vector<string>
Mesh::getNameOfNodeGroups( void )  const
{
	return _nodeGroupNames;
}

vector<MEDCoupling::DataArrayInt *>
Mesh::getNodeGroups( void )  const
{
	return _nodeGroups;
}

//----------------------------------------------------------------------
const IJKMesh&
Mesh::operator= ( const IJKMesh& mesh )
//----------------------------------------------------------------------
{
	_spaceDim = mesh.getSpaceDimension() ;
	_meshDim  = mesh.getMeshDimension() ;
    _name = mesh.getName();
	_numberOfNodes = mesh.getNumberOfNodes();
	_numberOfFaces = mesh.getNumberOfFaces();
	_numberOfCells = mesh.getNumberOfCells();
	_numberOfEdges = mesh.getNumberOfEdges();
    _indexFacePeriodicSet= mesh.isIndexFacePeriodicSet();
    if(_indexFacePeriodicSet)
        _indexFacePeriodicMap=mesh.getIndexFacePeriodic();
    
        _nxyz = mesh.getCellGridStructure() ;
        _dxyz=mesh.getDXYZ();
        _xMin=mesh.getXMin();
        _xMax=mesh.getXMax();
        _yMin=mesh.getYMin();
        _yMax=mesh.getYMax();
        _zMin=mesh.getZMin();
        _zMax=mesh.getZMax();

	_faceGroupNames = mesh.getNameOfFaceGroups() ;
	_faceGroups = mesh.getFaceGroups() ;
	_nodeGroupNames = mesh.getNameOfNodeGroups() ;
	_nodeGroups = mesh.getNodeGroups() ;

	MCAuto<MEDCouplingIMesh> m1=mesh.getMEDCouplingIMesh()->deepCopy();
	_mesh=m1;
	return *this;
}

bool Mesh::isIndexFacePeriodicSet() const
{
 return    _indexFacePeriodicSet;
}
//----------------------------------------------------------------------
double 
Mesh::minRatioVolSurf()
{
    return dx_min;
}
int 
Mesh::getMaxNbNeighbours(EntityType type) const
{
    double result=0;
    
    if (type==CELLS)
	{
        for(int i=0; i<_numberOfCells; i++)
            if(result < _cells[i].getNumberOfFaces())
                result=_cells[i].getNumberOfFaces();
	}
    else if(type==NODES)
	{
        for(int i=0; i<_numberOfNodes; i++)
            if(result < _nodes[i].getNumberOfEdges())
                result=_nodes[i].getNumberOfEdges();
	}
    else
		throw CdmathException("Mesh::getMaxNbNeighbours : entity type is not accepted. Should be CELLS or NODES");

    return result;
}
//----------------------------------------------------------------------
void
Mesh::writeVTK ( const std::string fileName ) const
//----------------------------------------------------------------------
{
	string fname=fileName+".vtu";
	_mesh->writeVTK(fname.c_str()) ;
}

//----------------------------------------------------------------------
void
Mesh::writeMED ( const std::string fileName ) const
//----------------------------------------------------------------------
{
	string fname=fileName+".med";
	MEDCoupling::WriteCMesh(fname.c_str(),_mesh,true);

	mu->decrRef();
}
