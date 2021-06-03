/*
    2���L���̂��`����N���X
    (�L����FiniteField�N���X�̔h��)
*/
#if !defined (___Class_BinaryFiniteField)
#define ___Class_BinaryFiniteField

#include "FiniteField.h"

class BinaryFiniteField : public FiniteField {
public :
    //===== �R���X�g���N�^ =====//
    explicit BinaryFiniteField();
    BinaryFiniteField(const int& value);    // �������t��
    BinaryFiniteField(const FiniteField& ff);

    //==== �ϊ��R���X�g���N�^ =====//
    operator FiniteField() {
        FiniteField ff;
        ff.setModular(2);
        if (getValueExist()) ff.setValue(value);
        return ff;
    }
};

//===== ���Z�q�̑��d��` =====//
bool operator==(const BinaryFiniteField& A, const BinaryFiniteField& B);
bool operator==(const vector<BinaryFiniteField>& A, const vector<BinaryFiniteField>& B);
bool operator!=(const BinaryFiniteField& A, const BinaryFiniteField& B);

#endif
