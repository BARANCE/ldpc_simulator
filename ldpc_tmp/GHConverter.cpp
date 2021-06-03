/*
    線形符号の検査行列と生成行列を相互に変換するためのクラス(実装)
*/

#include "GHConverter.h"

// 検査行列から生成行列に変換(組織符号のみ)
Matrix<BinaryFiniteField> GHConverter::ToGeneratorMatrix(const Matrix<BinaryFiniteField>& H) {
    int hr = H.row();
    int hc = H.col();
    Matrix<BinaryFiniteField> P = H.getSubMatrix(0, hr, 0, hc - hr).transposed();
    Matrix<BinaryFiniteField> I(hc - hr, hc - hr);
    I.setIdentity();
    Matrix<BinaryFiniteField> G = P.unionLeft(I);
    //cout << "rowSize = " << G.row() << ", colSize = " << G.col() << "\n";
    //cout << G;

    return G;
}

// 生成行列から検査行列に変換(組織符号のみ)
Matrix<BinaryFiniteField> GHConverter::ToCheckMatrix(const Matrix<BinaryFiniteField>& G) {
    int gr = G.row();
    int gc = G.col();
    Matrix<BinaryFiniteField> P = G.getSubMatrix(0, gr, gr, gc).transposed();
    Matrix<BinaryFiniteField> I(gc - gr, gc - gr);
    I.setIdentity();
    Matrix<BinaryFiniteField> H = P.unionRight(I);
    return H;
}
