/*
    ���x���W�@
*/

#if !defined (___Class_DensityEvolution)
#define ___Class_DensityEvolution

#include <math.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include "SPMatrix.h"
using namespace std;

class DensityEvolution{
public:
    //==== ���̓f�[�^�̍\�� ====//
    typedef struct ___DEINPUT {
        double fripProbability;    // ���]�m��
        long double borderProb;    // ���0�Ƃ݂Ȃ��m��
        unsigned long convergenceFailNumber;    // �������s����
        string optionalResultFileName;    // �����I�Ȏ������ʂ�ۑ�����t�@�C����
    } DEInput;
    typedef struct ___DESINPUT {
        long double zeroBorderProb;    // ���0�Ƃ݂Ȃ��m��
        unsigned long convergenceFailNumber;    // �������s����
        unsigned long binarySearchLoop;    // �񕪒T���̎��s��
        string checkMatrixFileName;    // �����s��̃t�@�C����
        string optionalResultFileName;    // �����I�Ȏ������ʂ�ۑ�����t�@�C����
    } DESInput;

    //==== �o�̓f�[�^�̍\�� ====//
    typedef struct ___DEOUTPUT {
        vector<long double> convProbList;    // �����I�����̊e�ϐ��m�[�h�̎����m����ێ�����x�N�g��
        bool failFrag;    // ��藦��0�Ɏ������Ȃ��ꍇ��true
        unsigned long convergenceNumber;    // ������
    } DEOutput;
    typedef struct ___DEWSOUTPUT {
        vector<long double> convProbList;    // �����I�����̊e�ϐ��m�[�h�̎����m����ێ�����x�N�g��
        bool failFrag;    // ��藦��0�Ɏ������Ȃ��ꍇ��true
        unsigned long convergenceNumber;    // ������
        Matrix<long double> valueChecker;    // �������Ƃ̕ϐ��m�[�h���b�Z�[�W�̌�藦�l(�����V���{�����M�m��)��ێ�
        vector<long double> maxValueChecker;    // �ő��藦
        vector<long double> aveValueChecker;    // ���ό�藦
    } DEWSOutput;
private:
    //==== BP�ŗp���郁�b�Z�[�W�o�b�t�@ ====//
    // BP�ŗp����ϐ��m�[�h���`�F�b�N�m�[�h�A�`�F�b�N�m�[�h���ϐ��m�[�h�̃��b�Z�[�W���Ǘ���������N���X
    class MessageBuffer {
        const SPMatrix* H;    // �����s��
    public:
        vector<vector<double> > variableBuf;    // �ϐ��m�[�h�̃o�b�t�@
        vector<vector<double> > checkBuf;    // �`�F�b�N�m�[�h�̃o�b�t�@

        class BufferOperationError{
            const MessageBuffer* ident;
        public:
            BufferOperationError(const MessageBuffer* a) : ident(a) {}
        };

        //==== �R���X�g���N�^ ====//
        MessageBuffer(const SPMatrix* checkMatrix) : H(checkMatrix), variableBuf(H->col()), checkBuf(H->row()) {
            int hr = H->row();
            int hc = H->col();
            const vector<vector<int> >* matCol = H->matCol;
            const vector<vector<int> >* matRow = H->matRow;
            for(int i = 0; i < hc; i++) {
                int size = (*matCol)[i].size();
                variableBuf[i].resize(size, 0.0);
            }
            for(int i = 0; i < hr; i++) {
                int size = (*matRow)[i].size();
                checkBuf[i].resize(size, 0.0);
            }

            // �f�o�O�p
            /*vector<int> testVec = (*matCol)[300];
            unsigned long sz = testVec.size();
            for(int i = 0; i < sz; i++) {
                int tv = testVec[i];
                cout << tv <<" ";
                variableToCheck(300, tv - 1, 1.0);
            }
            cout << "\n";
            cout << getCheckMessage(89, 300) << "\n";*/

        }

        //==== ���b�Z�[�W���M ====//
        int checkToVariable(unsigned long checkId, unsigned long variableId, double message) {    // checkId�Ԗڂ̃`�F�b�N�m�[�h����variableId�Ԗڂ̕ϐ��m�[�h�ւ̃��b�Z�[�W
            if(checkId < 0 || checkId >= (unsigned long)(H->row()) || variableId < 0 || variableId >= (unsigned long)(H->col())) throw BufferOperationError(this);
            long num = getCheckId(variableId, checkId);
            if(num == -1) {
                cout << "Error : " << checkId << "�Ԗڂ̃`�F�b�N�m�[�h�Ɨאڂ��Ă��Ȃ��ϐ��m�[�h" << variableId << "�ւ̃��b�Z�[�W�ł��B\n";
                throw BufferOperationError(this);
                return -1;    // ���M���s
            }
            variableBuf[variableId][num] = message;
            return 0;
        }
        int variableToCheck(unsigned long variableId, unsigned long checkId, double message) {    // variableId�Ԗڂ̕ϐ��m�[�h����checkId�Ԗڂ̃`�F�b�N�m�[�h�ւ̃��b�Z�[�W
            if(checkId < 0 || checkId >= (unsigned long)(H->row()) || variableId < 0 || variableId >= (unsigned long)(H->col())) throw BufferOperationError(this);
            long num = getVariableId(checkId, variableId);
            if(num == -1) {
                cout << "Error : " << variableId << "�Ԗڂ̕ϐ��m�[�h�Ɨאڂ��Ă��Ȃ��`�F�b�N�m�[�h" << checkId << "�ւ̃��b�Z�[�W�ł��B\n";
                throw BufferOperationError(this);
                return -1;    // ���M���s
            }
            checkBuf[checkId][num] = message;
            return 0;
        }

        //==== ���b�Z�[�W��M ====//
        double getCheckMessage(unsigned long checkId, unsigned long variableId) const {    // checkId�Ԗڂ̃`�F�b�N�m�[�h�ɂ����āAvariableId�Ԗڂ̕ϐ��m�[�h����̃��b�Z�[�W���擾
            long vid = getVariableId(checkId, variableId);
            if(vid == -1) {
                cout << "Error : " << checkId << "�Ԗڂ̃`�F�b�N�m�[�h�Ɨאڂ��Ă��Ȃ��ϐ��m�[�h" << variableId << "����̃��b�Z�[�W�͎�M�ł��܂���B\n";
                throw BufferOperationError(this);
                return 0.0;    // ��M���s
            }
            return checkBuf[checkId][vid];
        }
        double getVariableMessage(unsigned long variableId, unsigned long checkId) {    // varibaleId�Ԗڂ̕ϐ��m�[�h�ɂ����āAcheckId�Ԗڂ̃`�F�b�N�m�[�h����̃��b�Z�[�W���擾
            long vid = getCheckId(variableId, checkId);
            if(vid == -1) {
                cout << "Error : " << variableId << "�Ԗڂ̕ϐ��m�[�h�Ɨאڂ��Ă��Ȃ��`�F�b�N�m�[�h" << checkId << "����̃��b�Z�[�W�͎�M�ł��܂���B\n";
                throw BufferOperationError(this);
                return 0.0;    // ��M���s
            }
            return variableBuf[variableId][vid];
        }
        //==== �ڑ�����Ă���m�[�h�ɑΉ����郁�b�Z�[�W�o�b�t�@�x�N�^�̔ԍ����擾 ====//
        long getCheckId(unsigned long variableId, unsigned long checkId) const {    // variable�Ԗڂ̕ϐ��m�[�h�ɕt�^���ꂽcheckId�Ԗڂ̃`�F�b�N�m�[�h�̊i�[ID���擾
            const vector<vector<int> >* matCol = H->matCol;
            const vector<int>& vec = (*matCol)[variableId];
            int size = vec.size();
            int i = 0;
            for(i = 0; i < size; i++) {
                if(vec[i] == checkId + 1) {
                    break;
                }
            }
            if(i >= size) i = -1;
            return i;
        }
        long getVariableId(unsigned long checkId, unsigned long variableId) const {    // checkId�Ԗڂ̃`�F�b�N�m�[�h�ɕt�^���ꂽvariableId�Ԗڂ̕ϐ��m�[�h�̊i�[ID���擾
            const vector<vector<int> >* matRow = H->matRow;
            const vector<int>& vec = (*matRow)[checkId];
            int size = vec.size();
            int i = 0;
            for(i = 0; i < size; i++) {
                if(vec[i] == variableId + 1) {
                    break;
                }
            }
            if(i >= size) i = -1;
            return i;
        }
    };

    //==== �ϐ��� ====//
    const SPMatrix* H;    // �����s��
    double fp;    // ���]�m��
    MessageBuffer mb;    // ���b�Z�[�W�o�b�t�@
    clock_t startClock, nowClock, endClock;    // �^�C�}�[

public:
    //==== �R���X�g���N�^ ====//
    DensityEvolution(const SPMatrix* checkMatrix) : H(checkMatrix), mb(checkMatrix), fp(0.0) {}

    //==== ���s ====//
    long double execDESimulation(DESInput simulationParameter);    // ���͌����s��ɑ΂��A臒l�ƂȂ锽�]�m�������߂�V�~�����[�V���������s
    DEOutput execDE(DEInput parameter);    // �v���g�O���t��
    DEOutput execSubDE(DEInput parameter, vector<int> connectionArea, vector<int> blockEndPos);    // �v���g�O���t��
    DEWSOutput execDEwithStepCheck(DEInput parameter, string fileName, unsigned long defConvNum);
    long double borderProb(int valDim, int checkDim, int depth) const;
    bool deStep(int valDim, int checkDim, long double fripProbability) const;

    //==== ���x���W�@�̌v�Z�� ====//
    inline long double deeq(long double x, long double p, int checkDim, int valDim) const {
        long double retVal = p * pow((long double)1.0 - pow((long double)1.0 - x, checkDim - 1), valDim - 1);
        return retVal;
    }

    //==== �i���󋵂̕\�� ====//
    void printProgress(double progressRate) const;

    //==== ���ʂ�ۑ� ====//
    void saveDEData(string fileName, double epsilon, unsigned long loop, unsigned long girth) const;

};

#endif
