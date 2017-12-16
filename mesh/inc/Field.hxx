/*
 * field.hxx
 *
 *  Created on: 07 fevrier. 2012
 *      Authors: CDMAT groups
 */

#ifndef FIELD_HXX_
#define FIELD_HXX_

namespace MEDCoupling
{
  class MEDCouplingFieldDouble;
  class DataArrayDouble;
}

typedef enum
  {
    CELLS = 0,
    NODES = 1,
    FACES = 2
  } TypeField;

#include "DoubleTab.hxx"
#include "Vector.hxx"
#include "Mesh.hxx"

#include <MCAuto.hxx>

/**
 * Field class is defined by
 * - ........
 */

class Field
{
    public: //----------------------------------------------------------------
    /**
     * default constructor
     */
    Field ( void ) ;

    /**
    * constructor with data:
    * @param fieldName : name of the field
    * @param type : type of the field
    * @param mesh : mesh of the field
    * @param numberOfComponents : number of the component
    * @param time : time of the field
    */
    Field(const std::string fieldName, TypeField type, const Mesh& mesh, int numberOfComponents, double time) ;

    /**
    * constructor with data:
    * @param fieldName : name of the field
    * @param type : type of the field
    * @param mesh : mesh of the field
    * @param numberOfComponents : number of the component
    */
    Field( const std::string fieldName, TypeField type, const Mesh& mesh, int numberOfComponents) ;

    /**
    * constructor with data:
    * @param fieldName : name of the field
    * @param type : type of the field
    * @param mesh : mesh of the field
    */
    Field( const std::string fieldName, TypeField type, const Mesh& mesh) ;

    /**
    * destructor
    */
    ~Field ( void ) ;

    /**
    * constructor by copy
    * @param field : The Field object to be copied
    */
    Field ( const Field & field ) ;

    /**
     * constructor with data
     * @param filename : file name of field med file
     * @param fieldType: field type
     * @param fieldName: field name
     * @param iteration: iteration number (optional)
     * @param order:     order inside an iteration (optional)
     */
    Field( const std::string filename, TypeField fieldType,
           const std::string & fieldName = "",
           int iteration = -1, int order = -1);
  
    MEDCoupling::DataArrayDouble * getArray();

    void readFieldMed( const std::string & fileNameRadical,
                       TypeField type,
                       const std::string & fieldName = "",
                       int iteration = -1,
                       int order = -1) ;

    double& operator[] ( int ielem ) ;

    double operator[] ( int ielem ) const;

    double& operator() ( int ielem ) ;

    double operator() ( int ielem ) const;

    double& operator() ( int ielem, int jcomp ) ;

    double operator() ( int ielem, int jcomp ) const ;

    int getNumberOfComponents ( void ) const ;

    const double* getValues ( void ) const ;

    const std::string getName ( void ) const;

    const Mesh& getMesh ( void ) const ;

    int getNumberOfElements ( void ) const ;

    TypeField getTypeOfField ( void ) const ;

    /**
     * return the MEDCouplingField pointer
     * return _field
     */
    MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> getField ( void )  const ;

    void setFieldByMEDCouplingFieldDouble ( const MEDCoupling::MEDCouplingFieldDouble* field );

    void setFieldByDataArrayDouble ( const MEDCoupling::DataArrayDouble* array );

    DoubleTab getNormEuclidean( void ) const ;

    void setTime ( double time, int iter );

    Vector getValuesOnComponent(int compo) const ;

    Vector getValuesOnAllComponents(int elem) const ;

    int getSpaceDimension( void ) const;

    double getTime ( void ) const;

    void setName ( const std::string fieldName ) ;

    void setInfoOnComponent(int icomp, std::string nameCompo) ;

    std::string getInfoOnComponent(int icomp) const;

    /**
     * Computes all the components of the sum of values of the field multiplied by dual cell measures. In case of a field on cells, the dual mesh coincides with the underlying mesh
     * The field may be multicomponent so the result of the integral should be a vector
     * return the vector of numerical value of the integral of the field
     */
    Vector integral() const;

    /**
     * Computes the sum of values of a given component of the field multiplied by dual cell measures. In case of a field on cells, the dual mesh coincides with the underlying mesh
     * @param the index of the component of interest
     * return the numerical value of the integral of the field
     */
    double integral(int compId) const;

    /**
     * Computes for each component the sum of the absolute values of the field components multiplied by dual cell measures. In case of a field on cells, the dual mesh coincides with the underlying mesh.
     * The field may be multicomponent so the result of the integral should be a vector
     * return the vector of numerical value of the L1 norm of each component of the field
     */
    Vector normL1() const;

    /**
     * Computes all the components of the sum of squares of the values of the field components multiplied by dual cell measures. In case of a field on cells, the dual mesh coincides with the underlying mesh
     * The field may be multicomponent so the result of the integral should be a vector
     * return the vector of numerical value of the L2 norm of each component of the field
     */
    Vector normL2() const;

    const Field& operator= ( const Field& f ) ;

    Field operator+ ( const Field& f ) const ;

    Field operator- ( const Field& f ) const ;

    const Field& operator+= ( const Field& f ) ;

    const Field& operator-= ( const Field& f ) ;

    const Field& operator*= ( double s ) ;

    const Field& operator/= ( double s ) ;

    const Field& operator-= ( double s ) ;

    const Field& operator+= ( double s ) ;

    void writeVTK ( const std::string fileName, bool fromScratch=true ) const ;

    void writeMED ( const std::string fileName, bool fromScratch=true ) const ;

    void writeCSV ( const std::string fileName ) const ;

    friend Field operator* (double value , const Field& field ) ;

    friend Field operator* (const Field& field, double value ) ;

    friend Field operator/ (const Field& field, double value) ;

    friend std::ostream& operator<<(std::ostream& out, const Field& field ) ;

    protected: //----------------------------------------------------------------

    MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> _field;
    Mesh _mesh ;
    TypeField _typeField;

    private:

};

#endif /* Field_HXX_ */
