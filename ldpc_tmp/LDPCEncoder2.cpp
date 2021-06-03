#include "LDPCEncoder2.h"

//==== コンストラクタ ====//
LDPCEncoder2::LDPCEncoder2(const SPMatrix* checkMatrix) : H(checkMatrix) {
}

//==== 符号化 ====//
LDPCEncoder2::LENOutput LDPCEncoder2::encode(LENInput inputParameter) {
    //==== 初期化・前処理 ====//
    const vector<BinaryFiniteField> message = inputParameter.message;

    //==== 行列の変換 ====//
    Matrix<BinaryFiniteField> HN = (Matrix<BinaryFiniteField>)(*H);
    Matrix<BinaryFiniteField> newH = HN.getRENF_H();    // 検査行列
    GHConverter ghc;
    Matrix<BinaryFiniteField> newG = ghc.ToGeneratorMatrix(newH);    // 生成行列
    unsigned long hr = newH.row();
    unsigned long hc = newH.col();
    unsigned long k = newG.row();    // 次元

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
