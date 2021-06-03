/*
    ���`�����̌����s��Ɛ����s��𑊌݂ɕϊ����邽�߂̃N���X
*/

#if !defined(___Class_GHConverter)
#define ___Class_GHConverter

#include "Matrix.h"
#include "BinaryFiniteField.h"

class GHConverter {
public :
    //==== �ϊ���̍s����i�[����\���� ====//
    typedef struct ___CONVMAT {
        Matrix<BinaryFiniteField> H;    // �����s��
        Matrix<BinaryFiniteField> G;    // �����s��
    } ConvMatrix;

    //===== �R���X�g���N�^ =====//
    GHConverter(){}
    GHConverter(const GHConverter& ghc){}

    //===== �ϊ� =====//
    // �����s�񂩂琶���s��ɕϊ�
    Matrix<BinaryFiniteField> ToGeneratorMatrix(const Matrix<BinaryFiniteField>& H);
    // �����s�񂩂猟���s��ɕϊ�
    Matrix<BinaryFiniteField> ToCheckMatrix(const Matrix<BinaryFiniteField>& G);
};

#endif
