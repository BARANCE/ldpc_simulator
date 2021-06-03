#if !defined (___Class_LDPCEncoder2)
#define ___Class_LDPCEncoder2

#include "SPMatrix.h"
#include "GHConverter.h"

class LDPCEncoder2 {
public:
    //==== ���̓f�[�^�̍\�� ====//
    typedef struct ___LENINPUT {
        vector<BinaryFiniteField> message;    // ���b�Z�[�W
    } LENInput;

    //==== �o�̓f�[�^�̍\�� ====//
    typedef struct ___LENOUTPUT {
        Matrix<BinaryFiniteField> checkMatrix;    // �ϊ���̌����s��
        Matrix<BinaryFiniteField> generatorMatrix;    // �����s��
        vector<BinaryFiniteField> sendingWord;    // ���M��
    } LENOutput;

private:
    const SPMatrix* H;    // �����s��
public :
    //==== �R���X�g���N�^ ====//
    LDPCEncoder2(const SPMatrix* checkMatrix);

    //==== ������ ====//
    LENOutput encode(LENInput encodeParameter);

    //==== �����_���ȃ��b�Z�[�W���擾 ====//
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
