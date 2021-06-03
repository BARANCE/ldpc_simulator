/*
    ���`�����̌����s��Ɛ����s��𑊌݂ɕϊ����邽�߂̃N���X(����)
*/

#include "GHConverter.h"

// �����s�񂩂琶���s��ɕϊ�(�g�D�����̂�)
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

// �����s�񂩂猟���s��ɕϊ�(�g�D�����̂�)
Matrix<BinaryFiniteField> GHConverter::ToCheckMatrix(const Matrix<BinaryFiniteField>& G) {
    int gr = G.row();
    int gc = G.col();
    Matrix<BinaryFiniteField> P = G.getSubMatrix(0, gr, gr, gc).transposed();
    Matrix<BinaryFiniteField> I(gc - gr, gc - gr);
    I.setIdentity();
    Matrix<BinaryFiniteField> H = P.unionRight(I);
    return H;
}
