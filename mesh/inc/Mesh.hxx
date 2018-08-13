/*
 * mesh.hxx
 *
 *  Created on: 22 janv. 2012
 *      Authors: CDMAT
 */

#ifndef MESH_HXX_
#define MESH_HXX_

/**
 * Mesh class is defined by
 * - case 1: file name of mesh med file
 * - case 2: 1D : xmin and xmax and number of cells
 * - case 3: 2D : xmin, xmax, ymin and ymax and numbers of cells in x direction and y direction
 * - case 4 (not yet): 3D : xmin, xmax, ymin and ymax and numbers of cells in x direction, y direction and z direction
 */

namespace MEDCoupling
{
class MEDFileUMesh;
class MEDCouplingMesh;
class MEDCouplingIMesh;
class MEDCouplingUMesh;
}
#include <MCAuto.hxx>
#include "NormalizedGeometricTypes"

class Node;
class Cell;
class Face;

#include <vector>
#include <string>
#include <map>

class Mesh
{

public: //----------------------------------------------------------------
	/**
	 * default constructor
	 */
	Mesh ( void ) ;

	/**
	 * constructor with data to load a general unstructured mesh
	 * @param filename name of mesh file
	 * @param meshLevel : relative mesh dimension : 0->cells, 1->Faces etc
	 */
	Mesh ( const std::string filename, int meshLevel=0 ) ;

	/**
	 * constructor with data for a regular 1D grid 
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param nx : Number of cells in x direction
	 */
	Mesh( double xmin, double xmax, int nx, std::string meshName="MESH1D_Regular_Grid" ) ;

	/**
	 * constructor with data for an unstructured 1D mesh
	 * @param points : abscissas of the mesh nodes
	 */
	Mesh( std::vector<double> points, std::string meshName="MESH1D_unstructured" ) ;

	/**
	 * constructor with data for a regular 2D grid 
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param ymin : minimum y
	 * @param ymax : maximum y
	 * @param nx : Number of cells in x direction
	 * @param ny : Number of cells in y direction
	 */
	Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, std::string meshName="MESH2D_Regular_Rectangle_Grid") ;

	/**
	 * constructor with data for a regular 3D grid 
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param ymin : minimum y
	 * @param ymax : maximum y
	 * @param zmin : minimum z
	 * @param zmax : maximum z
	 * @param nx : Number of cells in x direction
	 * @param ny : Number of cells in y direction
	 * @param nz : Number of cells in z direction
	 */
	Mesh( double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz, std::string meshName="MESH3D_Regular_Cuboid_Grid") ;

	/**
	 * constructor with data for a regular 2D triangular mesh
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param ymin : minimum y
	 * @param ymax : maximum y
	 * @param nx : Number of square cells in x direction
	 * @param ny : Number of square cells in y direction
	 */
	Mesh( int policy, double xmin, double xmax, int nx, double ymin, double ymax, int ny, std::string meshName="MESH2D_Regular_Triangle") ;

	/**
	 * constructor with data for a regular 3D hexahedral mesh 
	 * @param xmin : minimum x
	 * @param xmax : maximum x
	 * @param ymin : minimum y
	 * @param ymax : maximum y
	 * @param zmin : minimum z
	 * @param zmax : maximum z
	 * @param nx : Number of cells in x direction
	 * @param ny : Number of cells in y direction
	 * @param nz : Number of cells in z direction
	 */
	Mesh( int policy, double xmin, double xmax, int nx, double ymin, double ymax, int ny, double zmin, double zmax, int nz, std::string meshName="MESH3D_Regular_Tetrahedra") ;

	Mesh( const MEDCoupling::MEDCouplingIMesh* mesh ) ;

	/**
	 * constructor with data
	 * @param filename : file name of mesh med file
	 * @param meshLevel : relative mesh dimension : 0->cells, 1->Faces etc
	 */
	void readMeshMed( const std::string filename, int meshLevel=0 ) ;

	/**
	 * constructor by copy
	 * @param mesh : The Mesh object to be copied
	 */
	Mesh ( const Mesh & mesh ) ;

	/**
	 * destructor
	 */
	~Mesh( void ) ;

	/**
	 * return mesh name
	 * @return _name
	 */
	std::string getName( void ) const ;

	/**
	 * return Space dimension
	 * @return _spaceDim
	 */
	int getSpaceDimension( void ) const ;

	/**
	 * return Mesh dimension
	 * @return _meshDim
	 */
	int getMeshDimension( void ) const ;

	/**
	 * return The nodes in this mesh
	 * @return _nodes
	 */
	Node* getNodes ( void ) const ;

	/**
	 * return The cells in this mesh
	 * @return _vertices
	 */
	Cell* getCells ( void ) const ;

	/**
	 * return the faces in this mesh
	 * @return _vertices
	 */
	Face* getFaces ( void ) const ;

	/**
	 * return the number of nodes in this mesh
	 * @return _numberOfNodes
	 */
	int getNumberOfNodes ( void )  const ;

	/**
	 * return the number of faces in this mesh
	 * @return _numberOfFaces
	 */
	int getNumberOfFaces ( void )  const ;

	/**
	 * return the number of cells in this mesh
	 * @return _numberOfCells
	 */
	int getNumberOfCells ( void )  const ;

	/**
	 * return the number of edges in this mesh
	 * @return _numberOfEdges
	 */
	int getNumberOfEdges ( void )  const ;

	/**
	 * return The cell i in this mesh
	 * @return _cells[i]
	 */
	Cell& getCell ( int i )  const ;

	/**
	 * return The face i in this mesh
	 * @return _faces[i]
	 */
	Face& getFace ( int i )  const ;

	/**
	 * return The node i in this mesh
	 * @return _nodes[i]
	 */
	Node& getNode ( int i )  const ;

	/**
	 * return number of cell in x direction (structured mesh)
	 * return _nX
	 */
	int getNx( void )  const ;

	/**
	 * return number of cell in y direction (structured mesh)
	 * return _nY
	 */
	int getNy( void )  const ;

	/**
	 * return number of cell in z direction (structured mesh)
	 * return _nZ
	 */
	int getNz( void )  const ;

	double getXMin( void )  const ;

	double getXMax( void )  const ;

	double getYMin( void )  const ;

	double getYMax( void )  const ;

	double getZMin( void )  const ;

	double getZMax( void )  const ;

	std::vector<double> getDXYZ() const ;// for structured meshes

	std::vector<int> getCellGridStructure() const;// for structured meshes

	/**
	 * surcharge operator =
	 * @param mesh : The Mesh object to be copied
	 */
	const Mesh& operator= ( const Mesh& mesh ) ;

	/**
	 * return the mesh MEDCoupling
	 * return _mesh
	 */
	MEDCoupling::MCAuto<MEDCoupling::MEDCouplingMesh> getMEDCouplingMesh ( void )  const ;

	/**
	 * return the list of group names
	 * return _groupNames
	 */
	std::vector<std::string> getNamesOfGroups( void )  const ;

	/**
	 * return the list of groups
	 * return _groups
	 */
	std::vector<MEDCoupling::MEDCouplingUMesh *> getGroups( void )  const ;

	/**
	 * write mesh in the VTK format
	 */
	void writeVTK ( const std::string fileName ) const ;

	/**
	 * write mesh in the MED format
	 */
	void writeMED ( const std::string fileName ) const ;

	void setGroupAtPlan(double value, int direction, double eps, std::string groupName) ;

	void setGroupAtFaceByCoords(double x, double y, double z, double eps, std::string groupName) ;

	/*
     * Functions to manage periodic boundary condition in square/cubic geometries 
     */
    int getIndexFacePeriodic(int indexFace) const ;
    void setPeriodicFaces();
    std::map<int,int> getIndexFacePeriodic( void ) const;
    bool isIndexFacePeriodicSet() const ;
    
	bool isBorderNode(int nodeid) const ;
	bool isBorderFace(int faceid) const ;
	
	bool isTriangular() const ;
	bool isTetrahedral() const ;
	bool isQuadrangular() const ;
	bool isHexahedral() const ;
    bool isStructured() const ;
    
	/**
	 * Compute the minimum value over all cells of the ratio cell perimeter/cell vaolume
	 */
    double minRatioVolSurf();
    
private: //----------------------------------------------------------------

	MEDCoupling::MEDCouplingUMesh*  setMesh( void ) ;

	void setGroups( const MEDCoupling::MEDFileUMesh* medmesh, MEDCoupling::MEDCouplingUMesh*  mu) ;

    std::string _name;
    
	/**
	 * Space dimension
	 */
	int _spaceDim ;

	/**
	 * Mesh dimension
	 */
	int _meshDim ;
    
    /*
     * Structured mes parameters
     */

    bool _isStructured;
    
	double _xMin;

	double _xMax;

	double _yMin;

	double _yMax;

	double _zMin;

	double _zMax;

	std::vector<int> _nxyz;

	std::vector<double> _dxyz;
	/*
	 * The nodes in this mesh.
	 */
	Node *_nodes;

	/*
	 * The number of nodes in this mesh.
	 */
	int _numberOfNodes;

	/*
	 * The faces in this mesh.
	 */
	Face *_faces;

	/*
	 * The numbers of faces in this mesh.
	 */
	int _numberOfFaces;

	/*
	 * The cells in this mesh.
	 */
	Cell *_cells;

	/*
	 * The number of cells in this mesh.
	 */
	int _numberOfCells;

	/*
	 * The number of edges in this mesh.
	 */
	int _numberOfEdges;

	/*
	 * The names of groups.
	 */
	std::vector<std::string> _groupNames;

	/*
	 * The list of groups.
	 */
	std::vector<MEDCoupling::MEDCouplingUMesh *> _groups;
	/*
	 * The mesh MEDCoupling
	 */
	MEDCoupling::MCAuto<MEDCoupling::MEDCouplingMesh> _mesh;
	std::vector< INTERP_KERNEL::NormalizedCellType > _eltsTypes;//List of cell types contained in the mesh
    
    /*
     * Tools to manage periodic boundary conditions in square/cube geometries
     */
     bool _indexFacePeriodicSet;
     std::map<int,int> _indexFacePeriodicMap;
};

#endif /* MESH_HXX_ */
