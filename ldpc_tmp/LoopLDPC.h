/*
    Loop��Ԍ���LDPC����
    �����̑т������ł���悤�ɂ�����Ԍ��������B

    ���̓p�����[�^�F
*/

#if !defined(___Class_LoopLDPC)
#define ___Class_LoopLDPC

#include "Matrix.h"
#include <exception>
#include <vector>

class LoopLDPC {
public:
    //==== ���̓f�[�^�̍\�� ====//
    typedef struct ___LBANDDATA {    // ��Ԍ����v���g�O���t1���̃f�[�^
        ___LBANDDATA() : l(3), r(6), L(7), conBandLeft(-1), conBandRight(-1), conPosLeft(0), conPosRight(0) {}
        int l;    // ��d��
        int r;    // �s�d��
        int L;    // ������
        int conBandLeft;    // �������ɐڑ�����т̔ԍ�(�ڑ����Ȃ��ꍇ��-1���w��)
        int conPosLeft;    // �������ɐڑ�����т̗�ʒu
        int conBandRight;    // �E�����ɐڑ�����т̔ԍ�(�ڑ����Ȃ��ꍇ��-1���w��)
        int conPosRight;    // �E�����ɐڑ�����т̗�ʒu
    } LBandData;

    //==== �G���[�N���X ====//
    class ParameterExeption {};

private:
    //==== �N���X�ϐ� ====//
    const vector<LBandData>& bandData;
    Matrix<BinaryFiniteField> mat;    // �s��{��
    vector<unsigned long> bsrp;    // �e�т̐擪�ʒu(�s)
    vector<unsigned long> bscp;    // �e�т̐퓬�ʒu(��)
public:
    //==== �R���X�g���N�^ ====//
    LoopLDPC(const vector<LBandData>& input);

    //==== ���J���\�b�h ====//
    // �p�����[�^�̓��e�ɏ]����LoopLDPC�����𐶐�����B
    Matrix<BinaryFiniteField> generate();

    // ���͂��ꂽ�т̏������J����B
    void printBandData();

private:
    //==== ����J���\�b�h ====//
    // �s��̃T�C�Y�����肷��
    unsigned long calculateRowSize(const vector<LBandData>& input);
    unsigned long calculateColSize(const vector<LBandData>& input);

    // �e�т̍s�񒆂̊J�n�ʒu(�s�Ɨ�)���擾���Absrp��bscp���X�V
    void calculateBandStartPos();

    // bandNum�Ԗڂ̑тɂ�����blockNum�Ԗڂ̃u���b�N�̐擪�ʒu(�s�Ɨ�)���擾
    unsigned long calculateBlockStartRowPos(int bandNum, int blockNum);
    unsigned long calculateBlockStartColPos(int bandNum, int blockNum);

    // bandNum�Ԗڂ̑тɂ�����blockNum�Ԗڂ̃u���b�N�𐶐�����
    void generateBlock(int bandNum, int blockNum);

    // bandNum�Ԗڂ̑т��s�񒆂ɔz�u����
    void generateBand(int bandNum);

    // �S�Ă̋�Ԍ����v���g�O���t�������s�񒆂ɔz�u����
    void generateBaseSCCProtograph();

    // ������s�̎Z�o
    unsigned long calculateBaseConRowLeft(int bandNum);
    unsigned long calculateBaseConRowRight(int bandNum);
    // �������̎Z�o
    unsigned long calculateBaseConColLeft(int bandNum);
    unsigned long calculateBaseConColRight(int bandNum);

    // �z�u�u���b�N�̒P�ʒ����Z�o
    unsigned long calculateConBlockLengthLeft(int bandNum);
    unsigned long calculateConBlockLengthRight(int bandNum);

    // step�X�e�b�v�ڂɂ�����c��z�u�r�b�g�����Z�o
    int calculateRemainSetBits(int bandNum, int step);

    // �z�u�u���b�N�����Z�o
    int calculateSetBlockNumLeft(int bandNum, int remainBits);
    int calculateSetBlockNumRight(int bandNum, int remainBits);

    // ������s�̍X�V
    unsigned long calculateBRowLeft(int bandNum, int step);
    unsigned long calculateBRowRight(int bandNum, int step);

    // �w��͈͂Ɍ����u���b�N��z�u
    void generateConnectBar(long bRow, long bColS, long bColE);

    // �����̔�������1��
    int connectWithStepLeft(int bandNum, int step);
    int connectWithStepRight(int bandNum, int step);

    // ���E�̌������̍쐬
    void generateConnectParts(int bandNum);

    // �������S�̂̍쐬
    void generateConnections();
};

#endif
