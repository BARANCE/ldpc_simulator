/*
    2元有限体を定義するクラス
    (有限体FiniteFieldクラスの派生)
*/
#if !defined (___Class_BinaryFiniteField)
#define ___Class_BinaryFiniteField

#include "FiniteField.h"

class BinaryFiniteField : public FiniteField {
public :
    //===== コンストラクタ =====//
    explicit BinaryFiniteField();
    BinaryFiniteField(const int& value);    // 初期化付き
    BinaryFiniteField(const FiniteField& ff);

    //==== 変換コンストラクタ =====//
    operator FiniteField() {
        FiniteField ff;
        ff.setModular(2);
        if (getValueExist()) ff.setValue(value);
        return ff;
    }
};

//===== 演算子の多重定義 =====//
bool operator==(const BinaryFiniteField& A, const BinaryFiniteField& B);
bool operator==(const vector<BinaryFiniteField>& A, const vector<BinaryFiniteField>& B);
bool operator!=(const BinaryFiniteField& A, const BinaryFiniteField& B);

#endif
