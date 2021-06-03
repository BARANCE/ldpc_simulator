#include "LDPCEncoder2.h"

//==== �R���X�g���N�^ ====//
LDPCEncoder2::LDPCEncoder2(const SPMatrix* checkMatrix) : H(checkMatrix) {
}

//==== ������ ====//
LDPCEncoder2::LENOutput LDPCEncoder2::encode(LENInput inputParameter) {
    //==== �������E�O���� ====//
    const vector<BinaryFiniteField> message = inputParameter.message;

    //==== �s��̕ϊ� ====//
    Matrix<BinaryFiniteField> HN = (Matrix<BinaryFiniteField>)(*H);
    Matrix<BinaryFiniteField> newH = HN.getRENF_H();    // �����s��
    GHConverter ghc;
    Matrix<BinaryFiniteField> newG = ghc.ToGeneratorMatrix(newH);    // �����s��
    unsigned long hr = newH.row();
    unsigned long hc = newH.col();
    unsigned long k = newG.row();    // ����

    Matrix<BinaryFiniteField> message2(1, k);
    for(int i = 0; i < k; i++) {
        message2[0][i] = message[i];
    }
    Matrix<BinaryFiniteField> codeWord2 = message2 * newG;
    vector<BinaryFiniteField> codeWord(hc, 0);
    for(int i = 0; i < hc; i++) {
        codeWord[i] = codeWord2[0][i];
    }

    LENOutput retObj = {newH, newG, codeWord};
    return retObj;

}
