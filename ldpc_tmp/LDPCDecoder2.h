#if !defined (___Class_LDPCDecoder2)
#define ___Class_LDPCDecoder2

#include <vector>
//#include <math.h>
#include "SPMatrix.h"
#include "CombinationArray.h"
using namespace std;

static const double ERASURE = -100.0;    // 消失シンボル
static const double ONE = 1.0;
static const double ZERO = 0.0;

class LDPCDecoder2 {
public:
    //==== 通信路 ====//
    typedef enum ___ECHANNEL {
        C_AWGN, C_BSC, C_BEC
    } EChannel;

    //==== 誤りの種類を表す ====//
    typedef enum ___EERRORTYPE {
        DECODE_SUCCESS, OVER_MAX_CONVERGENCE, BAD_CONVERGENCE
    } EErrorType;

    //==== 入力データの構造 ====//
    typedef struct ___ESPINPUT {
        vector<double> receptionWord;
        unsigned long convergenceFailNumber;    // 収束失敗条件
    } ESPInput;

    //==== 検査条件の構造 ====//
    typedef struct ___ECDLIST {
        int convCheckFrag;    // 収束状況をチェックするかどうか
        int convCheckSize;    // 収束状況をチェックするためのバッファのサイズ
    } ECDList;

    //==== 出力データの構造 ====//
    typedef struct ___ESPOUTPUT {
        vector<BinaryFiniteField> decodeWord;    // 推定符号語
        bool failFrag;    // 訂正が成功したら0,収束失敗したら1
        unsigned long convergenceNumber;    // 収束回数(収束失敗の場合は収束失敗条件に等しい値)
        EErrorType errorType;
    } ESPOutput;

    // 消失訂正用
    typedef struct ___ERSPOUTPUT {
        vector<double> decodeWord;    // 推定符号語
        bool failFrag;    // 訂正が成功したら0,収束失敗したら1
        unsigned long convergenceNumber;    // 収束回数(収束失敗の場合は収束失敗条件に等しい値)
        EErrorType errorType;
        vector<vector<double> > convCheckVec;    // 収束状況をチェックするベクトル
    } ERSPOutput;
private:
    //==== メッセージバッファ ====//
    // BPで用いる変数ノード→チェックノード、チェックノード→変数ノードのメッセージを管理する内部クラス
    class MessageBuffer {
        const SPMatrix* H;    // 検査行列
    public:
        vector<vector<double> > variableBuf;    // 変数ノードのバッファ
        vector<vector<double> > checkBuf;    // チェックノードのバッファ

        class BufferOperationError{
            const MessageBuffer* ident;
        public:
            BufferOperationError(const MessageBuffer* a) : ident(a) {}
        };

        //==== コンストラクタ ====//
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

            // デバグ用
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

        //==== メッセージ送信 ====//
        int checkToVariable(unsigned long checkId, unsigned long variableId, double message) {    // checkId番目のチェックノードからvariableId番目の変数ノードへのメッセージ
            if(checkId < 0 || checkId >= (unsigned long)(H->row()) || variableId < 0 || variableId >= (unsigned long)(H->col())) throw BufferOperationError(this);
            long num = getCheckId(variableId, checkId);
            if(num == -1) {
                cout << "Error : " << checkId << "番目のチェックノードと隣接していない変数ノード" << variableId << "へのメッセージです。\n";
                throw BufferOperationError(this);
                return -1;    // 送信失敗
            }
            variableBuf[variableId][num] = message;
            return 0;
        }
        int variableToCheck(unsigned long variableId, unsigned long checkId, double message) {    // variableId番目の変数ノードからcheckId番目のチェックノードへのメッセージ
            if(checkId < 0 || checkId >= (unsigned long)(H->row()) || variableId < 0 || variableId >= (unsigned long)(H->col())) throw BufferOperationError(this);
            long num = getVariableId(checkId, variableId);
            if(num == -1) {
                cout << "Error : " << variableId << "番目の変数ノードと隣接していないチェックノード" << checkId << "へのメッセージです。\n";
                throw BufferOperationError(this);
                return -1;    // 送信失敗
            }
            checkBuf[checkId][num] = message;
            return 0;
        }

        //==== メッセージ受信 ====//
        double getCheckMessage(unsigned long checkId, unsigned long variableId) const {    // checkId番目のチェックノードにおいて、variableId番目の変数ノードからのメッセージを取得
            long vid = getVariableId(checkId, variableId);
            if(vid == -1) {
                cout << "Error : " << checkId << "番目のチェックノードと隣接していない変数ノード" << variableId << "からのメッセージは受信できません。\n";
                throw BufferOperationError(this);
                return 0.0;    // 受信失敗
            }
            return checkBuf[checkId][vid];
        }
        double getVariableMessage(unsigned long variableId, unsigned long checkId) {    // varibaleId番目の変数ノードにおいて、checkId番目のチェックノードからのメッセージを取得
            long vid = getCheckId(variableId, checkId);
            if(vid == -1) {
                cout << "Error : " << variableId << "番目の変数ノードと隣接していないチェックノード" << checkId << "からのメッセージは受信できません。\n";
                throw BufferOperationError(this);
                return 0.0;    // 受信失敗
            }
            return variableBuf[variableId][vid];
        }
        //==== 接続されているノードに対応するメッセージバッファベクタの番号を取得 ====//
        long getCheckId(unsigned long variableId, unsigned long checkId) const {    // variable番目の変数ノードに付与されたcheckId番目のチェックノードの格納IDを取得
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
        long getVariableId(unsigned long checkId, unsigned long variableId) const {    // checkId番目のチェックノードに付与されたvariableId番目の変数ノードの格納IDを取得
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

    //==== 変数類 ====//
    const SPMatrix* H;    // 検査行列
    EChannel channel;    // 通信路
    MessageBuffer mb;    // メッセージバッファ
    double SNR;    // SN比(dB)
    double fp;    // 反転確率
public:

    //==== コンストラクタ ====//
    LDPCDecoder2(const SPMatrix* checkMatrix);

    //==== 通信路の選択 ====//
    void setChannel(EChannel channelName);

    //==== 通信路に関するパラメータ設定 ====//
    void setSNR(double rate);    // AWGN通信路のSN比を設定する
    void setFlipProbability(double flipProbability);    // 二元対称通信路の反転確率を設定する

    //==== 復号 ====//
    ESPOutput hardDecision(ESPInput decodeParameter);    // 硬判定
    ESPOutput sumProductDecode(ESPInput decodeParameter);    // sum-product復号
    ESPOutput minSumDecode(ESPInput decodeParameter, double scale);    // min-sum復号
    ESPOutput modMinSumDecode(ESPInput decodeParameter, double scale, vector<int>* assumptionVector);    // 修正min-sum復号
    ESPOutput iterativeMinSumDecode(ESPInput decodeParameter, double scale, const vector<int>* shortenIds);    // 固定値反復修正型min-sum復号
    ESPOutput shortenMinSumDecode(ESPInput decodeParameter, double scale, double shortenRate);

    ERSPOutput sumProductDecodeOnBEC(ESPInput decodeParameter, ECDList researchParameter);


private:
    //==== メッセージの一括計算 ====//
    double sumMessageFromVariableNode(unsigned long checkId) const;
    double sumMessageFromCheckNode(unsigned long variableId) const;

    //==== アークハイパボリックタンジェント ====//
    double atanh(double x) const;

    double absolute(double x) const;

    //==== sign関数 ====//
    double sign(double x) const;

    //==== gallagerのf関数 ====//
    double gallagerf(double x) const;

    //==== パリティ検査 ====//
    bool parityCheck(const vector<BinaryFiniteField>* nHat) const;    // 一般用
    bool isZero(const vector<BinaryFiniteField>* nHat) const;    // 速度重視用

    template <class Type> void printVector(const vector<Type>* vec) const {
        unsigned long size = vec->size();
        for(int i = 0; i < size; i++) {
            cout << (*vec)[i] << " ";
        }
        cout << "\n";
    }
};

#endif
