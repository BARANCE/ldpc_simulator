#if !defined (___Class_LDPCEncoder2)
#define ___Class_LDPCEncoder2

#include "SPMatrix.h"
#include "GHConverter.h"

class LDPCEncoder2 {
public:
    //==== 入力データの構造 ====//
    typedef struct ___LENINPUT {
        vector<BinaryFiniteField> message;    // メッセージ
    } LENInput;

    //==== 出力データの構造 ====//
    typedef struct ___LENOUTPUT {
        Matrix<BinaryFiniteField> checkMatrix;    // 変換後の検査行列
        Matrix<BinaryFiniteField> generatorMatrix;    // 生成行列
        vector<BinaryFiniteField> sendingWord;    // 送信語
    } LENOutput;

private:
    const SPMatrix* H;    // 検査行列
public :
    //==== コンストラクタ ====//
    LDPCEncoder2(const SPMatrix* checkMatrix);

    //==== 符号化 ====//
    LENOutput encode(LENInput encodeParameter);

    //==== ランダムなメッセージを取得 ====//
    static vector<BinaryFiniteField> getRandomMessage(unsigned long size) {
        vector<BinaryFiniteField> retVec(size);
        for(int i = 0; i < size; i++) {
            double rnd = genrand_real2();
            if(rnd < 0.5) {
                retVec[i] = 0;
            } else {
                retVec[i] = 1;
            }
        }
        return retVec;
    }
};

#endif
