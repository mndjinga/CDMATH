#include "Vector.hxx"
#include "CdmathException.hxx"

#include <cmath>

//----------------------------------------------------------------------
Vector::Vector( void )
//----------------------------------------------------------------------
{

}

//----------------------------------------------------------------------
Vector::~Vector( void )
//----------------------------------------------------------------------
{

}

//----------------------------------------------------------------------
Vector::Vector( int numberOfRows ):
Matrix(numberOfRows,1)
//----------------------------------------------------------------------
{

}

int Vector::size()
{
	return getNumberOfRows();
}
double&
Vector::operator()(int i)
{
	return Matrix::operator()(i,0);
}

double
Vector::operator()(int i) const
{
	return Matrix::operator()(i,0);
}

double&
Vector::operator[](int i)
{
	return Matrix::operator()(i,0);
}
double
Vector::operator[](int i) const
{
	return Matrix::operator()(i,0);
}

double
Vector::operator* (const Vector& vector) const
{
	double res=0.;
	int numberOfRows=getNumberOfRows();
	if(numberOfRows!= vector.getNumberOfRows())
		throw CdmathException("Vector::operator* vectors should have the same dimension for scalar product");
	for(int i=0; i<numberOfRows; i++)
	{
		res=res+Matrix::operator()(i,0)*vector(i);
	}
	return res;
}

Vector
Vector::innerProduct (const Vector& vector) const
{
	double res=0.;
	int numberOfRows=getNumberOfRows();
	if(numberOfRows!= vector.getNumberOfRows())
		throw CdmathException("Vector::operator* vectors should have the same dimension for scalar product");
	for(int i=0; i<numberOfRows; i++)
	{
		res=res+Matrix::operator()(i,0)*vector(i);
	}
	return res;
}

Vector
Vector::crossProduct (const Vector& vector) const
{
	  int numberOfRows1 = getNumberOfRows();
	  int numberOfRows2 = vector.getNumberOfRows();
	  if(numberOfRows1!= 3 || numberOfRows2!= 3 )
			throw CdmathException("Vector::operator* vectors should have the dimension 3 for cross-product");

	  Vector res(3);
	  res(0) = Matrix::operator()(1,0) * vector(2) - Matrix::operator()(2,0) * vector(1);
	  res(1) = Matrix::operator()(2,0) * vector(0) - Matrix::operator()(0,0) * vector(2);
	  res(2) = Matrix::operator()(0,0) * vector(1) - Matrix::operator()(1,0) * vector(0);

	  return res;
}

Matrix Vector::tensProduct (const Vector& vector) const
{
	Matrix res(getNumberOfRows(),vector.getNumberOfRows());
	for(int i=0;i<getNumberOfRows();i++)
		for(int j=0;j<vector.getNumberOfRows();j++)
			res(i,j)=Matrix::operator()(i,0)*vector(j);
	return res;
}
double Vector::norm() const
{
	double norm = 0.0;
	int dim=getNumberOfRows();
	for(int i=0; i<dim; i++)
		norm += Matrix::operator()(i,0)*Matrix::operator()(i,0);
	return sqrt(norm);
}

Vector
operator+ (const Vector& vector1, const Vector& vector2)
{
  int numberOfRows = vector1.getNumberOfRows();
  if(numberOfRows!= vector2.getNumberOfRows())
		throw CdmathException("Vector::operator+ vectors should have the same dimension for addition");
  Vector res(numberOfRows);
  for (int i=0; i<numberOfRows; i++)
	  res(i) = vector1(i) + vector2(i);
   return res;
}

Vector
operator- (const Vector& vector1, const Vector& vector2)
{
  int numberOfRows = vector1.getNumberOfRows();
  if(numberOfRows!= vector2.getNumberOfRows())
		throw CdmathException("Vector::operator+ vectors should have the same dimension for substraction");
  Vector res(numberOfRows);
  for (int i=0; i<numberOfRows; i++)
	  res(i) = vector1(i) - vector2(i);
   return res;
}


Vector
operator* (double value , const Vector& vector )
{
	int numberOfRows = vector.getNumberOfRows();
	Vector res(numberOfRows);
	for (int i=0; i<numberOfRows; i++)
		  res(i) = vector(i)*value;
   return res;
}

Vector
operator* (const Vector& vector, double value )
{
	int numberOfRows = vector.getNumberOfRows();
	Vector res(numberOfRows);
	for (int i=0; i<numberOfRows; i++)
		  res(i) = vector(i)*value;
   return res;
}

Vector
operator/ (const Vector& vector, double value)
{
	int numberOfRows = vector.getNumberOfRows();
	Vector res(numberOfRows);
	for (int i=0; i<numberOfRows; i++)
		  res(i) = vector(i)/value;
   return res;
}

Matrix
operator^(const Vector& vector1, const Vector& vector2)
{
	Matrix res(vector1.getNumberOfRows(),vector2.getNumberOfRows());
	for(int i=0;i<vector1.getNumberOfRows();i++)
		for(int j=0;j<vector2.getNumberOfRows();j++)
			res(i,j)=vector1(i)*vector1(j);
	return res;
}

Vector
operator% (const Vector& vector1, const Vector& vector2)
{
  int numberOfRows1 = vector1.getNumberOfRows();
  int numberOfRows2 = vector2.getNumberOfRows();
  if(numberOfRows1!= 3 || numberOfRows2!= 3 )
		throw CdmathException("Vector::operator* vectors should have the dimension 3 for cross-product");

  Vector res(3);
  res(0) = vector1(1) * vector2(2) - vector1(2) * vector2(1);
  res(1) = vector1(2) * vector2(0) - vector1(0) * vector2(2);
  res(2) = vector1(0) * vector2(1) - vector1(1) * vector2(0);

  return res;
}
