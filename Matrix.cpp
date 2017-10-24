#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "Matrix.h"

#define MATMATOPT(A, O, B) for (unsigned i = 0; i < r; ++i) {  \
                            for (unsigned j = 0; j < c; ++j)    \
                                A[i][j] O B[i][j];   \
                           }

#define MATSCAOPT(A, O, B) for (unsigned i = 0; i < r; ++i) {  \
                            for (unsigned j = 0; j < c; ++j)    \
                                A[i][j] O B;   \
                           }

//#define mat_DEBUG
const double mat_TOL = 0.0000001;

Matrix::Matrix(const unsigned r, const unsigned c): r(r), c(c)
{
#ifdef mat_DEBUG
    std::cout << "con: " << this << "\n";
#endif
    if (r && c) {
        e = (double**)malloc(sizeof(double*) * r + sizeof(double) * r * c);
        double *p = (double*)&e[r];
        for (unsigned i = 0; i < r; ++i, p += c) {
            e[i] = p;
            memset(e[i], 0, sizeof(double) * c);
        }
    }
    else {
        e = NULL;
    }
}

void Matrix::copyMatrix(const Matrix &m)
{
    if ((r = m.r) && (c = m.c)) {
        const unsigned s = sizeof(double*) * r + sizeof(double) * r * c;
        e = (double**)malloc(s);
        memcpy(e, m.e, s);
        double *p = (double*)&e[r];
        for (unsigned i = 0; i < r; ++i, p += c)
            e[i] = p;
    }
    else
        e = NULL;
}

Matrix::Matrix(const Matrix &m)
{
#ifdef mat_DEBUG
    std::cout << "copy con: " << this << "\n";
#endif
    copyMatrix(m);
}

Matrix::Matrix(const double *arr, const unsigned r, const unsigned c): r(r), c(c)
{
#ifdef mat_DEBUG
    std::cout << "con double*: " << this << "\n";
#endif
    if ((this->r = r) && (this->c = c)) {
        const int s = sizeof(double*) * r + sizeof(double) * r * c;
        e = (double**)malloc(s);
        double *p = (double*)&e[r];
        for (unsigned i = 0; i < r; ++i, p += c) {
            e[i] = p;
            for (unsigned j = 0; j < c; ++j)
                e[i][j] = arr[i * c + j];
        }
    }
    else
        e = NULL;
}

Matrix::~Matrix()
{
#ifdef mat_DEBUG
    std::cout << "des: " << this << "\n";
#endif
    if (e)
        free(e);
}

Matrix Matrix::sub(const unsigned r1, const unsigned r2, const unsigned c1, const unsigned c2)
{
    assert(r1 <= r2 && r2 <= r && c1 <= c2 && c2 <= c);
    Matrix m(r2 - r1, c2 - c1);
    for (unsigned i = r1; i < r2; ++i)
        memcpy(m.e[i - r1], &e[i][c1], sizeof(double) * (c2 - c1));
    return m;
}

void Matrix::changeDim(const unsigned r, const unsigned c)
{
    if (this->r == r && this->c == c)
        return;
    if (r == 0 || c == 0) {
        this->r = this->c = 0;
        if (e)
            free(e);
        return;
    }
    const int s = sizeof(double*) * r + sizeof(double) * r * c;
    double **new_e = (double**)malloc(s);
    double *p = (double*)&new_e[r];
    for (unsigned i = 0; i < r; ++i, p += c) {
        new_e[i] = p;
        for (unsigned j = 0; j < c; ++j)
            new_e[i][j] = (i < this->r && j < this->c) ? e[i][j] : 0.0;
    }
    free(e);
    e = new_e;
    this->r = r;
    this->c = c;
}

Matrix& Matrix::operator=(const Matrix &m)
{
#ifdef mat_DEBUG
    std::cout << "operator =: " << this << "\n";
#endif
    if (e)
        free(e);
    copyMatrix(m);
    return *this;
}

Matrix Matrix::t() const
{
#ifdef mat_DEBUG
    std::cout << "transpose: " << this << "\n";
#endif
    Matrix m(c, r);
    for (unsigned i = 0; i < r; ++i) {
        for (unsigned j = 0; j < c; ++j)
            m.e[j][i] = e[i][j];
    }
    return m;
}

Matrix& Matrix::operator+=(const Matrix &m)
{
#ifdef mat_DEBUG
    std::cout << "operator +=: " << this << "\n";
#endif
    assert(r == m.r && c == m.c);
    MATMATOPT(e, +=, m.e);
    return *this;
}

Matrix& Matrix::operator+=(const double d)
{
#ifdef mat_DEBUG
    std::cout << "operator += double: " << this << "\n";
#endif
    MATSCAOPT(e, +=, d);
    return *this;
}

Matrix& Matrix::operator-=(const Matrix &m)
{
#ifdef mat_DEBUG
    std::cout << "operator -=: " << this << "\n";
#endif
    assert(r == m.r && c == m.c);
    MATMATOPT(e, -=, m.e);
    return *this;
}

Matrix& Matrix::operator-=(const double d)
{
#ifdef mat_DEBUG
    std::cout << "operator -= double: " << this << "\n";
#endif
    MATSCAOPT(e, -=, d);
    return *this;
}

Matrix& Matrix::operator*=(const Matrix &m)
{
#ifdef mat_DEBUG
    std::cout << "operator *=: " << this << "\n";
#endif
    assert(c == m.r);
    Matrix mulM(r, m.c);
    for (unsigned i = 0; i < r; ++i) {
        for (unsigned j = 0; j < m.c; ++j) {
            for (unsigned p = 0; p < c; ++p)
                mulM[i][j] += e[i][p] * m.e[p][j];
        }
    }
    *this = mulM;
    return *this;
}

Matrix& Matrix::operator*=(const double d)
{
#ifdef mat_DEBUG
    std::cout << "operator *= double: " << this << "\n";
#endif
    MATSCAOPT(e, *=, d);
    return *this;
}

void Matrix::rowOpt(const unsigned pv, const double d, const unsigned tg, const unsigned startC)
{
    assert(pv < r && tg < r);
    double *pvRow = (double*)malloc(sizeof(double) * c);
    memcpy(pvRow, e[pv], sizeof(double) * c);
    for (unsigned j = startC; j < c; ++j)
        e[tg][j] += pvRow[j] * d;
    free(pvRow);
}

void Matrix::rowSwap(const unsigned r1, const unsigned r2)
{
     assert(r2 < r && r1 < r);
     const unsigned sizeC = sizeof(double) * c;
     double *bufRow = (double*)malloc(sizeC);
     memcpy(bufRow, e[r1], sizeC);
     memcpy(e[r1], e[r2], sizeC);
     memcpy(e[r2], bufRow, sizeC);
     free(bufRow);
}

double Matrix::det() const
{
    assert(r == c);
    if (r == 0)
        return 0.0;
    Matrix m(*this);
    double prod = 1.0;
    for (unsigned i = 0; i < r; ++i) {
        unsigned check = i;
        for (; check < r; ++check) {
            if (fabs(m.e[check][i]) >= mat_TOL)
                break;
        }
        if (check == r)
            return 0.0;
        if (check != i) {
            m.rowSwap(i, check);
            prod = -prod;
        }
        for (unsigned j = i + 1; j < r; ++j)
            m.rowOpt(i, -m.e[j][i] / m.e[i][i], j, i);
        prod *= m.e[i][i];
    }
    return prod;
}

Matrix Matrix::inv(bool *nonSingular) const
{
    assert(r == c);
    if (r == 0) {
        if (nonSingular)
            *nonSingular = false;
        return Matrix();
    }
    Matrix lm = Matrix(*this), rm = eye(r);
    for (unsigned i = 0; i < r; ++i) {
        unsigned check = i;
        for (; check < r; ++check) {
            if (fabs(lm.e[check][i]) >= mat_TOL)
                break;
        }
        if (check == r) {
            if (nonSingular)
                *nonSingular = false;
            return Matrix();
        }
        if (check != i) {
            lm.rowSwap(i, check);
            rm.rowSwap(i, check);
        }
        for (unsigned j = 0; j < r; ++j) {
            const double d = ((i != j) ? -lm.e[j][i]: (1.0 - lm.e[i][i])) / lm.e[i][i];
            lm.rowOpt(i, d, j, i);
            rm.rowOpt(i, d, j);
        }
    }
    if (nonSingular)
        *nonSingular = true;
    return rm;
}

void Matrix::show()
{
    for (unsigned i = 0; i < r; ++i) {
        for (unsigned j = 0; j < c; ++j)
            printf("%g\t", e[i][j]);
        putchar('\n');
    }
}

//non-member Matrix functions

Matrix operator+(const Matrix &lhs, const Matrix &rhs)
{
#ifdef mat_DEBUG
    std::cout << "operator +: " << &lhs << ", " << &rhs << "\n";
#endif
    return Matrix(lhs) += rhs;
}

Matrix operator+(const Matrix &lhs, const double d)
{
#ifdef mat_DEBUG
    std::cout << "operator +: " << &lhs << ", " << d << "\n";
#endif
    return Matrix(lhs) += d;
}

Matrix operator+(const double d, const Matrix &rhs)
{
#ifdef mat_DEBUG
    std::cout << "operator +: " << d << ", " << &rhs << "\n";
#endif
    return Matrix(rhs) += d;
}

Matrix operator-(const Matrix &m)
{
#ifdef mat_DEBUG
    std::cout << "operator -: " << &m << "\n";
#endif
    Matrix negM(m);
    unsigned r, c;
    m.dim(r, c);
    MATMATOPT(negM, =, -m);
    return negM;
}

Matrix operator-(const Matrix &lhs, const Matrix &rhs)
{
#ifdef mat_DEBUG
    std::cout << "operator +: " << &lhs << ", " << &rhs << "\n";
#endif
    return Matrix(lhs) -= rhs;
}

Matrix operator-(const Matrix &lhs, const double d)
{
#ifdef mat_DEBUG
    std::cout << "operator -: " << &lhs << ", " << d << "\n";
#endif
    return Matrix(lhs) -= d;
}

Matrix operator-(const double d, const Matrix &rhs)
{
#ifdef mat_DEBUG
    std::cout << "operator -: " << d << ", " << &rhs << "\n";
#endif
    return -Matrix(rhs) += d;
}

Matrix operator*(const Matrix &lhs, const double d)
{
#ifdef mat_DEBUG
    std::cout << "operator *: " << &lhs << ", " << d << "\n";
#endif
    return Matrix(lhs) *= d;
}

Matrix operator*(const double d, const Matrix &rhs)
{
#ifdef mat_DEBUG
    std::cout << "operator *: " << d << ", " << &rhs << "\n";
#endif
    return Matrix(rhs) *= d;
}

Matrix dot(const Matrix &lhs, const Matrix &rhs)
{
#ifdef mat_DEBUG
    std::cout << "ewm: " << &lhs << ", " << &rhs << "\n";
#endif
    Matrix m(lhs);
    unsigned r, c, rr, cr;
    lhs.dim(r, c);
    rhs.dim(rr, cr);
    assert(r == rr && c == cr);
    MATMATOPT(m, *=, rhs);
    return m;
}

Matrix operator*(const Matrix &lhs, const Matrix &rhs)
{
#ifdef mat_DEBUG
    std::cout << "operator *: " << &lhs << ", " << &rhs << "\n";
#endif
    return Matrix(lhs) *= rhs;
}

Matrix eye(const unsigned s)
{
    Matrix m(s, s);
    for (unsigned i = 0; i < s; ++i)
        m[i][i] = 1.0;
    return m;
}
