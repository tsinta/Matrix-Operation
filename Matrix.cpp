#include <iostream>
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
const double mat_TOL = 0.0000000001;

Matrix::Matrix(const unsigned r, const unsigned c, const double d): r(r), c(c), ip(0)
{
#ifdef mat_DEBUG
    std::cout << "con: " << this << "\n";
#endif
    if (r && c) {
        e = (double**)malloc(sizeof(double*) * r + sizeof(double) * r * c);
        double *p = (double*)&e[r];
        for (unsigned i = 0; i < r; ++i, p += c) {
            e[i] = p;
            if (fabs(d) < mat_TOL)
                memset(e[i], 0, sizeof(double) * c);
            else {
                for (unsigned j = 0; j < c; ++j)
                    e[i][j] = d;
            }
        }
    }
    else {
        this->r = this->c = 0;
        e = NULL;
    }
}

void Matrix::copyMatrix(const Matrix &m)
{
    if (r && c) {
        const unsigned s = sizeof(double*) * r + sizeof(double) * r * c;
        e = (double**)malloc(s);
        memcpy(e, m.e, s);
        double *p = (double*)&e[r];
        for (unsigned i = 0; i < r; ++i, p += c)
            e[i] = p;
    }
    else {
        r = c = 0;
        e = NULL;
    }
}

Matrix::Matrix(const Matrix &m): r(m.r), c(m.c), ip(0)
{
#ifdef mat_DEBUG
    std::cout << "copy con: " << this << "\n";
#endif
    copyMatrix(m);
}

Matrix::Matrix(const double* const arr, const unsigned r, const unsigned c): r(r), c(c), ip(0)
{
#ifdef mat_DEBUG
    std::cout << "con double*: " << this << "\n";
#endif
    if (r && c) {
        const int s = sizeof(double*) * r + sizeof(double) * r * c;
        e = (double**)malloc(s);
        double *p = (double*)&e[r];
        for (unsigned i = 0; i < r; ++i, p += c) {
            e[i] = p;
            for (unsigned j = 0; j < c; ++j)
                e[i][j] = arr[i * c + j];
        }
    }
    else {
        this->r = this->c = 0;
        e = NULL;
    }
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
        e = NULL;
        return;
    }
    const int s = sizeof(double*) * r + sizeof(double) * r * c;
    double **new_e = (double**)malloc(s);
    double *p = (double*)&new_e[r];
    for (unsigned i = 0; i < r; ++i, p += c) {
        new_e[i] = p;
        if (i < this->r) {
            memcpy(new_e[i], e[i], sizeof(double) * ((c < this->c) ? c : this->c));
            if (c > this->c)
                memset(&new_e[i][this->c], 0, sizeof(double) * (c - this->c));
        }
        else
            memset(new_e[i], 0, sizeof(double) * c);
        
    }
    free(e);
    e = new_e;
    this->r = r;
    this->c = c;
}

void Matrix::erase(const unsigned start, const unsigned end, const int mType)
{
    assert(start <= end && ((mType == ROW && end <= r) || (mType == COL && end <= c)));
    if (start == end)
        return;
    if (start == 0 && ((mType == ROW && end == r) || (mType == COL && end == c))) {
        this->r = this->c = 0;
        if (e)
            free(e);
        e = NULL;
        return;
    }
    const unsigned new_r = (mType == ROW) ? r - (end - start) : r;
    const unsigned new_c = (mType == COL) ? c - (end - start) : c;
    double **new_e = (double**)malloc(sizeof(double*) * new_r + sizeof(double) * new_r * new_c);
    double *p = (double*)&new_e[new_r];
    for (unsigned i = 0; i < new_r; ++i, p += new_c) {
        new_e[i] = p;
        if (mType == ROW) {
            if (i < start)
                memcpy(new_e[i], e[i], sizeof(double) * new_c);
            else
                memcpy(new_e[i], e[i + end - start], sizeof(double) * new_c);
        }
        else {
            memcpy(new_e[i], e[i], sizeof(double) * start);
            memcpy(&new_e[i][start], &e[i][end], sizeof(double) * (c - end));
        }
    }
    free(e);
    e = new_e;
    r = new_r;
    c = new_c;
}

void Matrix::insert(const unsigned pos, const Matrix &m, const int mType)
{
    assert((mType == ROW && pos <= r) || (mType == COL && pos <= c));
    if (m.empty())
        return;
    if (e == NULL && pos == 0) {
        *this = Matrix(m);
        return;
    }
    unsigned ci, ri;
    m.dim(ri, ci);
    assert((mType == ROW && ci == c) || (mType == COL && ri == r));
    const unsigned new_r = (mType == ROW) ? r + ri : r;
    const unsigned new_c = (mType == COL) ? c + ci : c;
    double **new_e = (double**)malloc(sizeof(double*) * new_r + sizeof(double) * new_r * new_c);
    double *p = (double*)&new_e[new_r];
    for (unsigned i = 0; i < new_r; ++i, p += new_c) {
        new_e[i] = p;
        if (mType == ROW) {
            if (i < pos)
                memcpy(new_e[i], e[i], sizeof(double) * new_c);
            else if (pos <= i && i < pos + ri)
                memcpy(new_e[i], m.e[i - pos], sizeof(double) * new_c);
            else
                memcpy(new_e[i], e[i - ri], sizeof(double) * new_c);
        }
        else {
            memcpy(new_e[i], e[i], sizeof(double) * pos);
            memcpy(&new_e[i][pos], m.e[i], sizeof(double) * ci);
            memcpy(&new_e[i][pos + ci], &e[i][pos], sizeof(double) * (c - pos));
        }
    }
    free(e);
    e = new_e;
    r = new_r;
    c = new_c;
}

void Matrix::opt(const unsigned pv, const double d, const unsigned tg, const int mType, const unsigned start, const unsigned end)
{
    assert(mType == ROW || mType == COL);
    if (mType == ROW)
        rowOpt(pv, d, tg, start, end);
    else
        colOpt(pv, d, tg, start, end);
}

void Matrix::swap(const unsigned s1, const unsigned s2, const int mType, const unsigned start, const unsigned end)
{
    assert(mType == ROW || mType == COL);
    if (mType == ROW)
        rowSwap(s1, s2, start, end);
    else
        colSwap(s1, s2, start, end);
}

Matrix& Matrix::operator<<(const double d)
{
    e[ip / c][ip % c] = d;
    ++ip;
    return *this;
}

Matrix& Matrix::operator<<(const double *d)
{
    for (unsigned i = 0; ip < r * c; ++ip, ++i)
        e[ip / c][ip % c] = d[i];
    return *this;
}

Matrix& Matrix::operator=(const Matrix &m)
{
#ifdef mat_DEBUG
    std::cout << "operator =: " << this << "\n";
#endif
    if (e)
        free(e);
    r = m.r;
    c = m.c;
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

const Matrix& Matrix::operator+=(const Matrix &m)
{
#ifdef mat_DEBUG
    std::cout << "operator +=: " << this << "\n";
#endif
    assert(r == m.r && c == m.c);
    MATMATOPT(e, +=, m.e);
    return *this;
}

const Matrix& Matrix::operator+=(const double d)
{
#ifdef mat_DEBUG
    std::cout << "operator += double: " << this << "\n";
#endif
    MATSCAOPT(e, +=, d);
    return *this;
}

const Matrix& Matrix::operator-=(const Matrix &m)
{
#ifdef mat_DEBUG
    std::cout << "operator -=: " << this << "\n";
#endif
    assert(r == m.r && c == m.c);
    MATMATOPT(e, -=, m.e);
    return *this;
}

const Matrix& Matrix::operator-=(const double d)
{
#ifdef mat_DEBUG
    std::cout << "operator -= double: " << this << "\n";
#endif
    MATSCAOPT(e, -=, d);
    return *this;
}

const Matrix& Matrix::operator*=(const Matrix &m)
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
    return *this = mulM;
}

const Matrix& Matrix::operator*=(const double d)
{
#ifdef mat_DEBUG
    std::cout << "operator *= double: " << this << "\n";
#endif
    MATSCAOPT(e, *=, d);
    return *this;
}

const Matrix& Matrix::operator/=(const double d)
{
#ifdef mat_DEBUG
    std::cout << "operator /= double: " << this << "\n";
#endif
    MATSCAOPT(e, /=, d);
    return *this;
}

void Matrix::rowOpt(const unsigned pv, const double d, const unsigned tg, const unsigned startC, unsigned endC)
{
    assert(pv < r && tg < r && startC <= c);
    if (endC == UINT_MAX)
        endC = c;
    assert(endC <= c);
    for (unsigned j = startC; j < endC; ++j)
        e[tg][j] += e[pv][j] * d;
}

void Matrix::rowSwap(const unsigned r1, const unsigned r2, const unsigned startC, unsigned endC)
{
    assert(r1 < r && r2 < r && startC <= c);
    if (endC == UINT_MAX)
        endC = c;
    assert(endC <= c);
    if (r1 == r2)
        return;
    const unsigned sizeC = sizeof(double) * (endC - startC);
    double *bufRow = (double*)malloc(sizeC);
    memcpy(bufRow, &e[r1][startC], sizeC);
    memcpy(&e[r1][startC], &e[r2][startC], sizeC);
    memcpy(&e[r2][startC], bufRow, sizeC);
    free(bufRow);
}

void Matrix::colOpt(const unsigned pv, const double d, const unsigned tg, const unsigned startR, unsigned endR)
{
    assert(pv < c && tg < c && startR <= r);
    if (endR == UINT_MAX)
        endR = r;
    assert(endR <= r);
    for (unsigned i = 0; i < r; ++i)
        e[i][tg] += e[i][pv] * d;
}

void Matrix::colSwap(const unsigned c1, const unsigned c2, const unsigned startR, unsigned endR)
{
    assert(c1 < c && c2 < c && startR <= r);
    if (endR == UINT_MAX)
        endR = r;
    assert(endR <= r);
    if (c1 == c2)
        return;
    for (unsigned i = startR; i < endR; ++i) {
        const double d = e[i][c1];
        e[i][c1] = e[i][c2];
        e[i][c2] = d;
    }
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

Matrix Matrix::inv() const
{
    assert(r == c);
    if (r == 0)
        return Matrix();
    Matrix lm = Matrix(*this), rm = eye(r);
    for (unsigned i = 0; i < r; ++i) {
        unsigned check = i;
        for (; check < r; ++check) {
            if (fabs(lm.e[check][i]) >= mat_TOL)
                break;
        }
        if (check == r)
            return Matrix();
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
    return rm;
}

void Matrix::show() const
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
    unsigned r, c;
    rhs.dim(r, c);
    Matrix minusM(r, c, d);
    return minusM -= rhs;
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

Matrix operator/(const Matrix &lhs, const double d)
{
#ifdef mat_DEBUG
    std::cout << "operator /: " << &lhs << ", " << d << "\n";
#endif
    return Matrix(lhs) /= d;
}

Matrix operator/(const double d, const Matrix &rhs)
{
#ifdef mat_DEBUG
    std::cout << "operator /: " << d << ", " << &rhs << "\n";
#endif
    unsigned r, c;
    rhs.dim(r, c);
    Matrix divM(r, c, d);
    MATMATOPT(divM, /=, rhs);
    return divM;
}

Matrix dotDiv(const Matrix &lhs, const Matrix &rhs)
{
#ifdef mat_DEBUG
    std::cout << "ewd: " << &lhs << ", " << &rhs << "\n";
#endif
    Matrix m(lhs);
    unsigned r, c, rr, cr;
    lhs.dim(r, c);
    rhs.dim(rr, cr);
    assert(r == rr && c == cr);
    MATMATOPT(m, /=, rhs);
    return m;
}

Matrix eye(const unsigned s)
{
    Matrix m(s, s);
    for (unsigned i = 0; i < s; ++i)
        m[i][i] = 1.0;
    return m;
}

//Solve linear equations

Matrix conjgrad(const Matrix &a, const Matrix &b, Matrix x, unsigned maxIter)
{
    if (x.empty())
        x = Matrix(b.dim(ROW), 1);
    assert(a.dim(COL) == b.dim(ROW) && b.dim(ROW) == x.dim(ROW) && b.dim(COL) == 1 && x.dim(COL) == 1);
    Matrix r = b - a * x;
    Matrix p = r;
    double rsold = (r.t() * r)[0][0];
    for (unsigned i = 0; i < maxIter; ++i) {
        const Matrix ap = a * p;
        const double alpha = rsold / (p.t() * ap)[0][0];
        x += alpha * p;
        r -= alpha * ap;
        const double rsnew = (r.t() * r)[0][0];
        if (sqrt(rsnew) < mat_TOL)
            break;
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    }
    return x;
}

static double loss(const Matrix &a, const Matrix &b, const Matrix &x)
{
    //loss = 1/2 * sum((ax - b).^2)
    Matrix m = a * x - b;
    m = dot(m, m);
    double sum = 0.0;
    for (unsigned i = 0; i < m.dim(ROW); ++i)
        sum += m[i][0];
    return sum / 2.0;
}

static Matrix gradient(const Matrix &a, const Matrix &b, const Matrix &x)
{
    //gradient j = sum(aij * (ai * x - yi)
    const unsigned r = a.dim(ROW), c = a.dim(COL);
    const Matrix ax_b = a * x - b;
    Matrix grad(c, 1);
    for (unsigned j = 0; j < c; ++j)
        for (unsigned i = 0; i < r; ++i)
            grad[j][0] += a[i][j] * ax_b[i][0];
    return grad;
}

Matrix graddesc(const Matrix &a, const Matrix &b, Matrix x, unsigned maxIter)
{
    if (x.empty())
        x = Matrix(b.dim(ROW), 1);
    assert(a.dim(COL) == b.dim(ROW) && b.dim(ROW) == x.dim(ROW) && b.dim(COL) == 1 && x.dim(COL) == 1);
    double alpha = 1.0;
    double lossOld = loss(a, b, x);
    for (unsigned i = 0; i < maxIter; ++i) {
        const Matrix xNew = x - alpha * gradient(a, b, x);
        const double lossNew = loss(a, b, xNew);
        if (fabs(lossOld - lossNew) < mat_TOL * mat_TOL || lossOld < mat_TOL * mat_TOL) {
            //std::cout << i << std::endl;
            break;
        }
        /*if (i % 10000 == 0 || i < 100) {
            std::cout << i << ":" << lossNew << ";" << lossOld << ";" << alpha << std::endl;
        }*/
        if (lossOld > lossNew) {
            lossOld = lossNew;
            x = xNew;
        }
        else
            alpha *= 0.95;
    }
    return x;
}
