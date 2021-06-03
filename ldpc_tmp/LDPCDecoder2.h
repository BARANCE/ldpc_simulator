#if !defined (___Class_LDPCDecoder2)
#define ___Class_LDPCDecoder2

#include <vector>
//#include <math.h>
#include "SPMatrix.h"
#include "CombinationArray.h"
using namespace std;

static const double ERASURE = -100.0;    // �����V���{��
static const double ONE = 1.0;
static const double ZERO = 0.0;

class LDPCDecoder2 {
public:
    //==== �ʐM�H ====//
    typedef enum ___ECHANNEL {
        C_AWGN, C_BSC, C_BEC
    } EChannel;

    //==== ���̎�ނ�\�� ====//
    typedef enum ___EERRORTYPE {
        DECODE_SUCCESS, OVER_MAX_CONVERGENCE, BAD_CONVERGENCE
    } EErrorType;

    //==== ���̓f�[�^�̍\�� ====//
    typedef struct ___ESPINPUT {
        vector<double> receptionWord;
        unsigned long convergenceFailNumber;    // �������s����
    } ESPInput;

    //==== ���������̍\�� ====//
    typedef struct ___ECDLIST {
        int convCheckFrag;    // �����󋵂��`�F�b�N���邩�ǂ���
        int convCheckSize;    // �����󋵂��`�F�b�N���邽�߂̃o�b�t�@�̃T�C�Y
    } ECDList;

    //==== �o�̓f�[�^�̍\�� ====//
    typedef struct ___ESPOUTPUT {
        vector<BinaryFiniteField> decodeWord;    // ���蕄����
        bool failFrag;    // ����������������0,�������s������1
        unsigned long convergenceNumber;    // ������(�������s�̏ꍇ�͎������s�����ɓ������l)
        EErrorType errorType;
    } ESPOutput;

    // ���������p
    typedef struct ___ERSPOUTPUT {
        vector<double> decodeWord;    // ���蕄����
        bool failFrag;    // ����������������0,�������s������1
        unsigned long convergenceNumber;    // ������(�������s�̏ꍇ�͎������s�����ɓ������l)
        EErrorType errorType;
        vector<vector<double> > convCheckVec;    // �����󋵂��`�F�b�N����x�N�g��
    } ERSPOutput;
private:
    //==== ���b�Z�[�W�o�b�t�@ ====//
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
    EChannel channel;    // �ʐM�H
    MessageBuffer mb;    // ���b�Z�[�W�o�b�t�@
    double SNR;    // SN��(dB)
    double fp;    // ���]�m��
public:

    //==== �R���X�g���N�^ ====//
    LDPCDecoder2(const SPMatrix* checkMatrix);

    //==== �ʐM�H�̑I�� ====//
    void setChannel(EChannel channelName);

    //==== �ʐM�H�Ɋւ���p�����[�^�ݒ� ====//
    void setSNR(double rate);    // AWGN�ʐM�H��SN���ݒ肷��
    void setFlipProbability(double flipProbability);    // �񌳑Ώ̒ʐM�H�̔��]�m����ݒ肷��

    //==== ���� ====//
    ESPOutput hardDecision(ESPInput decodeParameter);    // �d����
    ESPOutput sumProductDecode(ESPInput decodeParameter);    // sum-product����
    ESPOutput minSumDecode(ESPInput decodeParameter, double scale);    // min-sum����
    ESPOutput modMinSumDecode(ESPInput decodeParameter, double scale, vector<int>* assumptionVector);    // �C��min-sum����
    ESPOutput iterativeMinSumDecode(ESPInput decodeParameter, double scale, const vector<int>* shortenIds);    // �Œ�l�����C���^min-sum����
    ESPOutput shortenMinSumDecode(ESPInput decodeParameter, double scale, double shortenRate);

    ERSPOutput sumProductDecodeOnBEC(ESPInput decodeParameter, ECDList researchParameter);


private:
    //==== ���b�Z�[�W�̈ꊇ�v�Z ====//
    double sumMessageFromVariableNode(unsigned long checkId) const;
    double sumMessageFromCheckNode(unsigned long variableId) const;

    //==== �A�[�N�n�C�p�{���b�N�^���W�F���g ====//
    double atanh(double x) const;

    double absolute(double x) const;

    //==== sign�֐� ====//
    double sign(double x) const;

    //==== gallager��f�֐� ====//
    double gallagerf(double x) const;

    //==== �p���e�B���� ====//
    bool parityCheck(const vector<BinaryFiniteField>* nHat) const;    // ��ʗp
    bool isZero(const vector<BinaryFiniteField>* nHat) const;    // ���x�d���p

    template <class Type> void printVector(const vector<Type>* vec) const {
        unsigned long size = vec->size();
        for(int i = 0; i < size; i++) {
            cout << (*vec)[i] << " ";
        }
        cout << "\n";
    }
};

#endif
