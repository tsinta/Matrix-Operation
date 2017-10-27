#ifndef MATRIX_CLAUDE
#define MATRIX_CLAUDE

#include <cassert>
#include <cstdlib>

enum {ROW, COL};

class Matrix
{
public:
    Matrix(const unsigned r = 0, const unsigned c = 0, const double d = 0.0);
    Matrix(const Matrix &m);
    Matrix(const double* const arr, const unsigned r = 0, const unsigned c = 0);
    ~Matrix();
    
    bool empty() const {return e == NULL || r == 0 || c == 0;}
    Matrix sub(const unsigned r1, const unsigned r2, const unsigned c1, const unsigned c2); //sub matrix : r1 <= r < r2; c1 <= c < c2
    void dim(unsigned &r, unsigned &c) const { r = this->r; c = this->c; }
    unsigned dim(const int type = ROW) const { return (type == COL) ? c : r; }
    void changeDim(const unsigned r, const unsigned c);
    void erase(const unsigned start, const unsigned end, const int mType = ROW);    //erase : start= < (r or c) < end
    void insert(const unsigned pos, const Matrix &m, const int mType = ROW);    //insert m before row(column) pos
    void opt(const unsigned pv, const double d, const unsigned tg, const int mType = ROW, const unsigned start = 0, const unsigned end = UINT_MAX);
    void swap(const unsigned s1, const unsigned s2, const int mType = ROW, const unsigned start = 0, const unsigned end = UINT_MAX);
    
    Matrix& operator<<(const double d);
    Matrix& operator<<(const double *d);
    void reset(unsigned p = 0)  //for operator<<
    {
        ip = p;
    }
    Matrix& operator=(const Matrix &m);
    double*& operator[](const unsigned idx) const
    {
        assert(idx < r);
        return e[idx];
    }
    
    Matrix t() const;
    
    const Matrix& operator+=(const Matrix &m);
    const Matrix& operator+=(const double d);
    
    const Matrix& operator-=(const Matrix &m);
    const Matrix& operator-=(const double d);
    
    const Matrix& operator*=(const Matrix &m);
    const Matrix& operator*=(const double d);
    
    const Matrix& operator/=(const double d);
    
    double det() const;
    Matrix inv() const;
    
    void show() const;
private:
    double **e;
    unsigned r, c;
    unsigned ip;    //for operator<<
    void copyMatrix(const Matrix &m);
    void rowOpt(const unsigned pv, const double d, const unsigned tg, const unsigned startC = 0, unsigned endC = UINT_MAX); //row tg += row pv * d
    void rowSwap(const unsigned r1, const unsigned r2, const unsigned startC = 0, unsigned endC = UINT_MAX);
    void colOpt(const unsigned pv, const double d, const unsigned tg, const unsigned startR = 0, unsigned endR = UINT_MAX); //col tg += col pv * d
    void colSwap(const unsigned c1, const unsigned c2, const unsigned startR = 0, unsigned endR = UINT_MAX);
};

//non-member Matrix functions
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

Matrix operator/(const Matrix &lhs, const double d);
Matrix operator/(const double d, const Matrix &rhs);
Matrix dotDiv(const Matrix &lhs, const Matrix &rhs);    //Element-Wise Division

Matrix eye(const unsigned s);   //identity matrix

//Solve linear equations
Matrix conjgrad(const Matrix &a, const Matrix &b, Matrix x = Matrix(), unsigned maxIter = 1000000); //solve ax = b; a must be symmetric and positive definite
Matrix graddesc(const Matrix &a, const Matrix &b, Matrix x = Matrix(), unsigned maxIter = 1000000); //solve ax = b; not promise to find the correct answer

#endif