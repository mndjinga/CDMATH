/*
 * SparseMatrixPetsc.cxx
 *
 *  Created on: 04/11/2017
 *      Author: mndjinga
 */

#include "SparseMatrixPetsc.hxx"
#include "CdmathException.hxx"

#include <cstring>

using namespace std;

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc()
//----------------------------------------------------------------------
{
	_numberOfColumns=0;
	_numberOfRows=0;
	_numberOfNonZeros=0;
	_isSparseMatrix=true;
	_mat=NULL;
	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);
}

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc( int numberOfRows, int numberOfColumns)
//----------------------------------------------------------------------
{
	_numberOfRows = numberOfRows;
	_numberOfColumns=numberOfColumns;
	_isSparseMatrix=true;
	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);
	MatCreateSeqAIJ(MPI_COMM_SELF,_numberOfRows,_numberOfColumns,PETSC_DEFAULT,NULL,&_mat);
}

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc( Mat mat )
//----------------------------------------------------------------------
{
	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);
	_isSparseMatrix=true;
	_mat=mat;
	//extract number of row and column
	MatGetSize(mat,&_numberOfRows,&_numberOfColumns);

	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;
}

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc( int numberOfRows, int numberOfColumns, int nnz )
//----------------------------------------------------------------------
{
	_numberOfRows = numberOfRows;
	_numberOfColumns=numberOfColumns;
	_numberOfNonZeros=nnz;
	_isSparseMatrix=true;
	_mat=NULL;
	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);
	MatCreateSeqAIJ(MPI_COMM_SELF,_numberOfRows,_numberOfColumns,_numberOfNonZeros,NULL,&_mat);
}

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc( int blockSize, int numberOfRows, int numberOfColumns, int nnz )
//----------------------------------------------------------------------
{
	_numberOfRows = numberOfRows;
	_numberOfColumns=numberOfColumns;
	_numberOfNonZeros=nnz;
	_isSparseMatrix=true;
	_mat=NULL;
	PetscInitialize(0, (char ***)"", PETSC_NULL, PETSC_NULL);
	MatCreateSeqBAIJ(MPI_COMM_SELF,blockSize, _numberOfRows,_numberOfColumns,_numberOfNonZeros,NULL,&_mat);
}

//----------------------------------------------------------------------
SparseMatrixPetsc::SparseMatrixPetsc(const SparseMatrixPetsc& matrix)
//----------------------------------------------------------------------
{
	_isSparseMatrix=matrix.isSparseMatrix();
	MatDuplicate(matrix.getPetscMatrix(), MAT_COPY_VALUES,&_mat);
	MatGetSize(_mat,&_numberOfRows,&_numberOfColumns);	
	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(_mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;
}

SparseMatrixPetsc
SparseMatrixPetsc::transpose() const
{
	Mat mattranspose;
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	MatTranspose(_mat,MAT_INITIAL_MATRIX, &mattranspose);
	return SparseMatrixPetsc(mattranspose);
}

void
SparseMatrixPetsc::setValue( int i, int j, double value )
{
	MatSetValues(_mat,1, &i, 1, &j, &value, INSERT_VALUES);
}

void
SparseMatrixPetsc::addValue( int i, int j, double value )
{
	MatSetValues(_mat,1, &i, 1, &j, &value, ADD_VALUES);
}

void
SparseMatrixPetsc::setValue( int i, int j, Matrix M  )
{
    int I,J;
    for (int k=0; k<M.getNumberOfRows(); k++)
        for (int l=0; l<M.getNumberOfColumns(); l++)
        {
            I=i+k;
            J=j+l;
            MatSetValues(_mat,1, &I, 1, &J, &M(k,l), INSERT_VALUES);
        }
}

void
SparseMatrixPetsc::addValue( int i, int j, Matrix M  )
{
    int I,J;
    for (int k=0; k<M.getNumberOfRows(); k++)
        for (int l=0; l<M.getNumberOfColumns(); l++)
        {
            I=i+k;
            J=j+l;
            MatSetValues(_mat,1, &I, 1, &J, &M(k,l), ADD_VALUES);
        }
}

void
SparseMatrixPetsc::setValuesBlocked( int i, int j, Matrix M  )
{
    int blockSize;
    MatGetBlockSize(_mat,&blockSize);
    if(blockSize!=M.getNumberOfRows() || blockSize!=M.getNumberOfColumns())
        throw CdmathException("SparseMatrixPetsc::setValuesBlocked : matrix size is different from sparse matrix block structure");
    double petscValues[blockSize*blockSize];
    for (int k=0; k<M.getNumberOfRows(); k++)
        for (int l=0; l<M.getNumberOfColumns(); l++)
            petscValues[k*blockSize+l]=M(k,l);
    MatSetValuesBlocked(_mat,1, &i, 1, &j, petscValues, INSERT_VALUES);
}

void
SparseMatrixPetsc::addValuesBlocked( int i, int j, Matrix M  )
{
    int blockSize;
    MatGetBlockSize(_mat,&blockSize);
    if(blockSize!=M.getNumberOfRows() || blockSize!=M.getNumberOfColumns())
        throw CdmathException("SparseMatrixPetsc::addValuesBlocked : matrix size is different from sparse matrix block structure");
    double petscValues[blockSize*blockSize];
    for (int k=0; k<M.getNumberOfRows(); k++)
        for (int l=0; l<M.getNumberOfColumns(); l++)
            petscValues[k*blockSize+l]=M(k,l);
    MatSetValuesBlocked(_mat,1, &i, 1, &j, petscValues, ADD_VALUES);
}

//----------------------------------------------------------------------
double
SparseMatrixPetsc::operator()( int i, int j ) const
//----------------------------------------------------------------------
{
	double res;
	int idxm=i,idxn=j;
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	MatGetValues(_mat,1,&idxm,1, &idxn,&res);
	return res;
}

//----------------------------------------------------------------------
SparseMatrixPetsc::~SparseMatrixPetsc()
//----------------------------------------------------------------------
{
	if(&_mat != NULL)
		MatDestroy(&_mat);
	//PetscFinalize();
}

Mat
SparseMatrixPetsc::getPetscMatrix() const
{
	return (_mat);
}

bool
SparseMatrixPetsc::containsPetscMatrix() const
{
	return true;
}
//----------------------------------------------------------------------
const SparseMatrixPetsc&
SparseMatrixPetsc::operator= ( const SparseMatrixPetsc& matrix )
//----------------------------------------------------------------------
{
	_isSparseMatrix=matrix.isSparseMatrix();
	MatDuplicate(matrix.getPetscMatrix(), MAT_COPY_VALUES,&_mat);
	MatGetSize(_mat,&_numberOfRows,&_numberOfColumns);	
	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(_mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;
	return (*this);
}

SparseMatrixPetsc
operator+ (const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2)
{
	int numberOfRows = matrix1.getNumberOfRows();
	int numberOfColumns = matrix1.getNumberOfColumns();
	int numberOfRows2 = matrix2.getNumberOfRows();
	int numberOfColumns2 = matrix2.getNumberOfColumns();

	if(numberOfRows2!=numberOfRows || numberOfColumns2!=numberOfColumns)
	{
		string msg="SparseMatrixPetsc::operator()+(const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2): number of rows or columns of the matrices is different!";
		throw CdmathException(msg);
	}

	Mat mat1=matrix1.getPetscMatrix();
	Mat mat2=matrix2.getPetscMatrix();
	Mat mat;
	MatDuplicate(mat1, MAT_COPY_VALUES,&mat);
	MatAXPY(mat,1,mat2,DIFFERENT_NONZERO_PATTERN);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc
operator- (const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2)
{
	int numberOfRows = matrix1.getNumberOfRows();
	int numberOfColumns = matrix1.getNumberOfColumns();
	int numberOfRows2 = matrix2.getNumberOfRows();
	int numberOfColumns2 = matrix2.getNumberOfColumns();

	if(numberOfRows2!=numberOfRows || numberOfColumns2!=numberOfColumns)
	{
		string msg="SparseMatrixPetsc::operator()-(const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2): number of rows or columns of the matrices is different!";
		throw CdmathException(msg);
	}
	Mat mat1=matrix1.getPetscMatrix();
	Mat mat2=matrix2.getPetscMatrix();
	Mat mat;
	MatDuplicate(mat1, MAT_COPY_VALUES,&mat);
	MatAXPY(mat,-1,mat2,DIFFERENT_NONZERO_PATTERN);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc
operator* (double value , const SparseMatrixPetsc& matrix )
{
	Mat mat;
	MatDuplicate(matrix.getPetscMatrix(), MAT_COPY_VALUES,&mat);
	MatScale(mat, value);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc
operator* (const SparseMatrixPetsc& matrix, double value )
{
	Mat mat;
	MatDuplicate(matrix.getPetscMatrix(), MAT_COPY_VALUES,&mat);
	MatScale(mat, value);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc
operator/ (const SparseMatrixPetsc& matrix, double value)
{
	if(value==0.)
	{
		string msg="SparseMatrixPetsc SparseMatrixPetsc::operator()/(const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2): division by zero";
		throw CdmathException(msg);
	}
	Mat mat;
	MatDuplicate(matrix.getPetscMatrix(), MAT_COPY_VALUES,&mat);
	MatScale(mat, 1/value);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc
operator*(const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2)
{
	Mat mat1=matrix1.getPetscMatrix();
	Mat mat2=matrix2.getPetscMatrix();
	Mat mat;
	MatMatMult(mat1, mat2, MAT_INITIAL_MATRIX,PETSC_DEFAULT,&mat);

	return SparseMatrixPetsc(mat);
}

SparseMatrixPetsc&
SparseMatrixPetsc::operator*= (const SparseMatrixPetsc& matrix)
{
	Mat mat1=matrix.getPetscMatrix();
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	Mat mat2;
	MatMatMult(_mat, mat1, MAT_INITIAL_MATRIX,PETSC_DEFAULT,&mat2);

	MatDestroy(&_mat);
	_mat=mat2;
	MatGetSize(_mat,&_numberOfRows,&_numberOfColumns);	
	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(_mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;
	return (*this);
}

Vector
SparseMatrixPetsc::operator* (const Vector& vec) const
{
	int numberOfRows=vec.getNumberOfRows();
	Vec X;

	VecCreate(PETSC_COMM_WORLD,&X);
	VecSetSizes(X,PETSC_DECIDE,numberOfRows);
	VecSetFromOptions(X);
	double value ;
	for (PetscInt i=0; i<numberOfRows; i++)
	{
		value = vec(i);
		VecSetValues(X,1,&i,&value,INSERT_VALUES);
	}

	VecAssemblyBegin(X);
	VecAssemblyEnd(X);

	Vec Y;
	VecDuplicate (X,&Y);
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
	MatMult(_mat,X,Y);

	Vector result(numberOfRows);
	for (PetscInt i=0; i<numberOfRows; i++)
	{
		VecGetValues(Y,1,&i,&value);
		result(i)=value;
	}
	return result;
}

SparseMatrixPetsc&
SparseMatrixPetsc::operator+= (const SparseMatrixPetsc& matrix)
{
	Mat mat1=matrix.getPetscMatrix();
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
	MatAXPY(_mat,1,mat1,DIFFERENT_NONZERO_PATTERN);

	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(_mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;

	return (*this);
}

SparseMatrixPetsc&
SparseMatrixPetsc::operator-= (const SparseMatrixPetsc& matrix)
{
	Mat mat1=matrix.getPetscMatrix();
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
	MatAXPY(_mat,-1,mat1,DIFFERENT_NONZERO_PATTERN);

	//extract an upper bound for the total number of non zero coefficients
	MatInfo info;
	MatGetInfo(_mat,MAT_LOCAL,&info);
	_numberOfNonZeros = info.nz_allocated;

	return (*this);
}

SparseMatrixPetsc&
SparseMatrixPetsc::operator*= (double value)
{
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
	MatScale(_mat, value);
	return (*this);
}

SparseMatrixPetsc&
SparseMatrixPetsc::operator/= (double value)
{
	if(value==0.)
	{
		string msg="SparseMatrixPetsc SparseMatrixPetsc::operator()/=(const SparseMatrixPetsc& matrix1, const SparseMatrixPetsc& matrix2): division by zero";
		throw CdmathException(msg);
	}
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
	MatScale(_mat, 1/value);
	return (*this);
}

void
SparseMatrixPetsc::viewMatrix() const 
{
    MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	MatView(_mat,PETSC_VIEWER_STDOUT_SELF);
}
double
SparseMatrixPetsc::getMatrixCoeff(int i, int j) const 
{
	double res;
	int idxm=i,idxn=j;
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	MatGetValues(_mat,1,&idxm,1, &idxn,&res);
	return res;
}

void 
SparseMatrixPetsc::diagonalShift(double lambda)
{
    MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

    MatShift(_mat, lambda);
}

void 
SparseMatrixPetsc::zeroEntries()
{
    MatZeroEntries(_mat);
}
