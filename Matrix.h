#ifndef MATRIX_CLAUDE
#define MATRIX_CLAUDE

#include <cassert>
#include <cstdlib>

class Matrix
{
public:
    Matrix(const unsigned r = 0, const unsigned c = 0);
    Matrix(const Matrix &m);
    Matrix(const double *arr, const unsigned r = 0, const unsigned c = 0);
    ~Matrix();
    
    bool empty() const {return e == NULL || r == 0 || c == 0;}
    Matrix sub(const unsigned r1, const unsigned r2, const unsigned c1, const unsigned c2); //sub matrix :r1 <= r < r2; c1 <= c < c2
    void dim(unsigned &r, unsigned &c) const { r = this->r; c = this->c; }
    void changeDim(const unsigned r, const unsigned c);
    
    Matrix& operator=(const Matrix &m);
    double*& operator[](const unsigned idx) const
    {
        assert(idx < r);
        return e[idx];
    }
    
    Matrix t() const;
    
    Matrix& operator+=(const Matrix &m);
    Matrix& operator+=(const double d);
    
    Matrix& operator-=(const Matrix &m);
    Matrix& operator-=(const double d);
    
    Matrix& operator*=(const Matrix &m);
    Matrix& operator*=(const double d);
    
    double det() const;
    Matrix inv() const;
    
    void show();
private:
    double **e;
    unsigned r, c;
    void copyMatrix(const Matrix &m);
    void rowOpt(const unsigned pv, const double d, const unsigned tg, const unsigned startC = 0);  //row tg += row pv * d
    void rowSwap(const unsigned r1, const unsigned r2);
};

Matrix operator+(const Matrix &lhs, const Matrix &rhs);
Matrix operator+(const Matrix &lhs, const double d);
Matrix operator+(const double d, const Matrix &rhs);

Matrix operator-(const Matrix &m);
Matrix operator-(const Matrix &lhs, const Matrix &rhs);
Matrix operator-(const Matrix &lhs, const double d);
Matrix operator-(const double d, const Matrix &rhs);

Matrix operator*(const Matrix &lhs, const double d);
Matrix operator*(const double d, const Matrix &rhs);
Matrix dot(const Matrix &lhs, const Matrix &rhs);   //Element-Wise Multiplication
Matrix operator*(const Matrix &lhs, const Matrix &rhs);

Matrix eye(const unsigned s);   //identity matrix

#endif