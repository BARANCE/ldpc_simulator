/*
    LDPC�����̕����V�~�����[�^
*/

#if !defined(___Class_DecodeSimulator2)
#define ___Class_DecodeSimulator2

#include "LDPCDecoder2.h"
#include <time.h>
//#include "mt19937ar-cok.cpp"

static int randFrag = 0;

class DecodeSimulator2 {
public:
    //==== �ʐM�H ====//
    typedef enum ___ECHANNEL {
        C_AWGN, C_BSC, C_BEC
    } EChannel;

    //==== �����@ ====//
    typedef enum ___EDMETHOD {
        C_SUM_PRODUCT, C_MIN_SUM, C_HARD_DECISION, C_MMIN_SUM, C_SHORTEN_MIN_SUM
    } EDMethod;

    //==== ���̓f�[�^�̍\�� ====//
    typedef struct ___DSINPUT {
        EDMethod decodingMethod;
        EChannel channel;
        string fileName;
        vector<BinaryFiniteField> cordWord;    // ���M��
        double SNR;    // SN��
        unsigned long maxLoop;    // �ő厎�s��
        unsigned long maxError;    // �ő����
        unsigned long convergenceFailNumber;    // �������s����(�ő唽����)
        double scale;    // min-sum�̃X�P�[�����O�萔(�s�v�ȏꍇ��0.0�ɂł�)
        double shortenRate;    // �V���[�g�����銄��(�s�v�ȏꍇ��0�ɂł�)
        int convCheckFrag;    // �����󋵂��`�F�b�N���邩�ǂ���
        int convCheckSize;    // �����󋵂��`�F�b�N����o�b�t�@�̃T�C�Y
    } DSInput;

    //==== ���ςƕ��U ====//
    typedef struct ___DSAVEVAR {
        double average;
        double variance;
    } DSAveVar;

    //==== �o�̓f�[�^�̍\�� ====//
    typedef struct ___DSOUTPUT {
        unsigned long loop;            // ���ۂ̎��s��
        double blockErrorRate;        // �u���b�N��藦
        double bitErrorRate;        // �r�b�g��藦
        double decodeSuccessRate;    // ����������
        double error1Rate;            // �������
        double error2Rate;            // �������s��
        double error3Rate;            // �������
        double averageConvergenceNumber;    // ���ώ�����
        unsigned long minConvergenceNumber;    // �ŏ�������
        unsigned long maxConvergenceNumber;    // �ő������
        vector<unsigned long> distributionOfNumberOfBitError;    // �������s���̌��r�b�g���̕��z
        double averageOfNumberOfBitError;    // �������s���̌��r�b�g���̕���
        //double varianceOfNumberOfBitError;    // �������s���̌��r�b�g���̕��U
        double averageOfGaussianNoise;    // �K�E�X�G���̕W�{����
        double varianceOfGaussianNoise;    // �K�E�X�G���̕W�{���U
        double averageOfBurstLength;    // �o�[�X�g���̕W�{����
        unsigned long execTime;    // ���s����

        //___DSOUTPUT() : loop(0), blockErrorRate(0.0), bitErrorRate(0.0), decodeSuccessRate(0.0), error1Rate(0.0), error2Rate(0.0), error3Rate(0.0), averageConvergenceNumber(0.0), varianceConvergenceNumber(0.0), minConvergenceNumber(0), maxConvergenceNumber(0), DistributionOfNumberOfBitError(0), averageOfNumberOfBitError(0.0), varianceOfNumberOfBitError(0.0) {}
    } DSOutput;

private:
    const SPMatrix* H;    // �����s��
    clock_t startClock, nowClock, endClock;    // �^�C�}�[
public:
    //==== �R���X�g���N�^ ====//
    DecodeSimulator2(const SPMatrix* checkMatrix);

    //==== �������� ====//
    DSOutput decodingSimulation(DSInput simulationParameter);

    //===== �G������ ====//
    double getGaussianNoise(double average, double variance) const;
    DSAveVar checkAverageAndVariance(const vector<double>* vec) const;

    //==== ������Ƒ��M�����ꂪ���������ǂ��� ====//
    bool checkDecodeWord(const vector<BinaryFiniteField>* codeWord, const vector<BinaryFiniteField>* recepWord) const;

    //==== �i���󋵂̕\�� ====//
    void printProgress(double progressRate) const;

    //==== �g���q����菜�� ====//
    string removeExtension(string str) const;

};

#endif
