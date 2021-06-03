/*
    密度発展法
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
    //==== 入力データの構造 ====//
    typedef struct ___DEINPUT {
        double fripProbability;    // 反転確率
        long double borderProb;    // 誤り0とみなす確率
        unsigned long convergenceFailNumber;    // 収束失敗条件
        string optionalResultFileName;    // 副次的な実験結果を保存するファイル名
    } DEInput;
    typedef struct ___DESINPUT {
        long double zeroBorderProb;    // 誤り0とみなす確率
        unsigned long convergenceFailNumber;    // 収束失敗条件
        unsigned long binarySearchLoop;    // 二分探索の実行回数
        string checkMatrixFileName;    // 検査行列のファイル名
        string optionalResultFileName;    // 副次的な実験結果を保存するファイル名
    } DESInput;

    //==== 出力データの構造 ====//
    typedef struct ___DEOUTPUT {
        vector<long double> convProbList;    // 反復終了時の各変数ノードの収束確率を保持するベクトル
        bool failFrag;    // 誤り率が0に収束しない場合はtrue
        unsigned long convergenceNumber;    // 収束回数
    } DEOutput;
    typedef struct ___DEWSOUTPUT {
        vector<long double> convProbList;    // 反復終了時の各変数ノードの収束確率を保持するベクトル
        bool failFrag;    // 誤り率が0に収束しない場合はtrue
        unsigned long convergenceNumber;    // 収束回数
        Matrix<long double> valueChecker;    // 反復ごとの変数ノードメッセージの誤り率値(消失シンボル送信確率)を保持
        vector<long double> maxValueChecker;    // 最大誤り率
        vector<long double> aveValueChecker;    // 平均誤り率
    } DEWSOutput;
private:
    //==== BPで用いるメッセージバッファ ====//
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
    double fp;    // 反転確率
    MessageBuffer mb;    // メッセージバッファ
    clock_t startClock, nowClock, endClock;    // タイマー

public:
    //==== コンストラクタ ====//
    DensityEvolution(const SPMatrix* checkMatrix) : H(checkMatrix), mb(checkMatrix), fp(0.0) {}

    //==== 実行 ====//
    long double execDESimulation(DESInput simulationParameter);    // 入力検査行列に対し、閾値となる反転確率を求めるシミュレーションを実行
    DEOutput execDE(DEInput parameter);    // プロトグラフ版
    DEOutput execSubDE(DEInput parameter, vector<int> connectionArea, vector<int> blockEndPos);    // プロトグラフ版
    DEWSOutput execDEwithStepCheck(DEInput parameter, string fileName, unsigned long defConvNum);
    long double borderProb(int valDim, int checkDim, int depth) const;
    bool deStep(int valDim, int checkDim, long double fripProbability) const;

    //==== 密度発展法の計算式 ====//
    inline long double deeq(long double x, long double p, int checkDim, int valDim) const {
        long double retVal = p * pow((long double)1.0 - pow((long double)1.0 - x, checkDim - 1), valDim - 1);
        return retVal;
    }

    //==== 進捗状況の表示 ====//
    void printProgress(double progressRate) const;

    //==== 結果を保存 ====//
    void saveDEData(string fileName, double epsilon, unsigned long loop, unsigned long girth) const;

};

#endif
