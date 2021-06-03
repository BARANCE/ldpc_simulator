/*
    線形符号の検査行列と生成行列を相互に変換するためのクラス
*/

#if !defined(___Class_GHConverter)
#define ___Class_GHConverter

#include "Matrix.h"
#include "BinaryFiniteField.h"

class GHConverter {
public :
    //==== 変換後の行列を格納する構造体 ====//
    typedef struct ___CONVMAT {
        Matrix<BinaryFiniteField> H;    // 検査行列
        Matrix<BinaryFiniteField> G;    // 生成行列
    } ConvMatrix;

    //===== コンストラクタ =====//
    GHConverter(){}
    GHConverter(const GHConverter& ghc){}

    //===== 変換 =====//
    // 検査行列から生成行列に変換
    Matrix<BinaryFiniteField> ToGeneratorMatrix(const Matrix<BinaryFiniteField>& H);
    // 生成行列から検査行列に変換
    Matrix<BinaryFiniteField> ToCheckMatrix(const Matrix<BinaryFiniteField>& G);
};

#endif
