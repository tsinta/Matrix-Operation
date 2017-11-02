#include "Matrix.h"
#include <iostream>

using namespace MatOpt;
using namespace std;

int main()
{
    //show identity matrix
    Matrix Indntity = eye(3);
    cout << "Indentity\n";
    Indntity.show();
    
    //matrix multiplication
    double a[3][2] = {{1, 2}, {10, 11}, {3, 4}};
    double b[2][3] = {{1, 5, 3}, {6 ,4, 2}};
    Matrix Ma((double*)a, 3, 2), Mb((double*)b, 2, 3);
    Matrix Mc;
    cout << "a*b=\n";
    Mc = Ma * Mb;
    Mc.show();
    
    cout << "c*5=\n";
    Mc = Mc * 5;
    Mc.show();
    
    //Transpose
    cout << "Transpose c\n";
    Mc.t().show();
    
    Mc[2][2] = 1;
    //Determinant
    cout << "Determinant\n";
    cout << "det = " << det(Mc) << endl;
    
    //Inverse
    cout << "Inverse\n";
    Matrix inv_c = inv(Mc);
    inv_c.show();
    
    //c*inv_c
    cout << "c*inv_c\n";
    (Mc * inv_c).show();
    
    //insert & erase
    double d[2][2] = {{1, 2}, {3, 4}};
    double e[2] = {5, 6};
    Matrix Md((double*)d, 2, 2), Me(e, 1, 2);
    Md.insert(1, Me, ROW);
    cout << "Insert\n";
    Md.show();
    
    Md.erase(0, 1, COL);
    cout << "Erase\n";
    Md.show();
}
