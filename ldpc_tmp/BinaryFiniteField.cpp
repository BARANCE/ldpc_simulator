/*
    2���L���̂��`����N���X(��������)
*/

#include "BinaryFiniteField.h"

//===== �R���X�g���N�^ =====//
BinaryFiniteField::BinaryFiniteField() : FiniteField(2, 0) {
    setModular(2);
}
BinaryFiniteField::BinaryFiniteField(const int& value) : FiniteField(2, value) {
}
BinaryFiniteField::BinaryFiniteField(const FiniteField& ff) : FiniteField() {
    setModular(2);
    if (ff.getValueExist()) setValue(ff.getValue());
}

//===== ���Z�q�̑��d��` =====//
bool operator==(const BinaryFiniteField& A, const BinaryFiniteField& B) {
    if (!A.getOperationPropriety()) return false;
    if (!B.getOperationPropriety()) return false;
    if (A.getValue() == B.getValue()) return true;
    else return false;
}
bool operator==(const vector<BinaryFiniteField>& A, const vector<BinaryFiniteField>& B) {
    int Asize = A.size();
    if (Asize != B.size()) return false;
    bool tmp = true;
    for(int i = 0; i < Asize; i++) {
        tmp = tmp && (A[i] == B[i]);
    }
    return tmp;
}
bool operator!=(const BinaryFiniteField& A, const BinaryFiniteField& B) {
    if (!A.getOperationPropriety()) return true;
    if (!B.getOperationPropriety()) return true;
    if (A.getValue() == B.getValue()) return false;
    else return true;
}
