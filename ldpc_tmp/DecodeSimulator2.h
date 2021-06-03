/*
    LDPC符号の復号シミュレータ
*/

#if !defined(___Class_DecodeSimulator2)
#define ___Class_DecodeSimulator2

#include "LDPCDecoder2.h"
#include <time.h>
//#include "mt19937ar-cok.cpp"

static int randFrag = 0;

class DecodeSimulator2 {
public:
    //==== 通信路 ====//
    typedef enum ___ECHANNEL {
        C_AWGN, C_BSC, C_BEC
    } EChannel;

    //==== 復号法 ====//
    typedef enum ___EDMETHOD {
        C_SUM_PRODUCT, C_MIN_SUM, C_HARD_DECISION, C_MMIN_SUM, C_SHORTEN_MIN_SUM
    } EDMethod;

    //==== 入力データの構造 ====//
    typedef struct ___DSINPUT {
        EDMethod decodingMethod;
        EChannel channel;
        string fileName;
        vector<BinaryFiniteField> cordWord;    // 送信語
        double SNR;    // SN比
        unsigned long maxLoop;    // 最大試行回数
        unsigned long maxError;    // 最大誤り回数
        unsigned long convergenceFailNumber;    // 収束失敗条件(最大反復回数)
        double scale;    // min-sumのスケーリング定数(不要な場合は0.0にでも)
        double shortenRate;    // ショートンする割合(不要な場合は0にでも)
        int convCheckFrag;    // 収束状況をチェックするかどうか
        int convCheckSize;    // 収束状況をチェックするバッファのサイズ
    } DSInput;

    //==== 平均と分散 ====//
    typedef struct ___DSAVEVAR {
        double average;
        double variance;
    } DSAveVar;

    //==== 出力データの構造 ====//
    typedef struct ___DSOUTPUT {
        unsigned long loop;            // 実際の試行回数
        double blockErrorRate;        // ブロック誤り率
        double bitErrorRate;        // ビット誤り率
        double decodeSuccessRate;    // 訂正成功率
        double error1Rate;            // 誤訂正率
        double error2Rate;            // 収束失敗率
        double error3Rate;            // 誤収束率
        double averageConvergenceNumber;    // 平均収束回数
        unsigned long minConvergenceNumber;    // 最小収束回数
        unsigned long maxConvergenceNumber;    // 最大収束回数
        vector<unsigned long> distributionOfNumberOfBitError;    // 訂正失敗時の誤りビット数の分布
        double averageOfNumberOfBitError;    // 訂正失敗時の誤りビット数の平均
        //double varianceOfNumberOfBitError;    // 訂正失敗時の誤りビット数の分散
        double averageOfGaussianNoise;    // ガウス雑音の標本平均
        double varianceOfGaussianNoise;    // ガウス雑音の標本分散
        double averageOfBurstLength;    // バースト長の標本平均
        unsigned long execTime;    // 実行時間

        //___DSOUTPUT() : loop(0), blockErrorRate(0.0), bitErrorRate(0.0), decodeSuccessRate(0.0), error1Rate(0.0), error2Rate(0.0), error3Rate(0.0), averageConvergenceNumber(0.0), varianceConvergenceNumber(0.0), minConvergenceNumber(0), maxConvergenceNumber(0), DistributionOfNumberOfBitError(0), averageOfNumberOfBitError(0.0), varianceOfNumberOfBitError(0.0) {}
    } DSOutput;

private:
    const SPMatrix* H;    // 検査行列
    clock_t startClock, nowClock, endClock;    // タイマー
public:
    //==== コンストラクタ ====//
    DecodeSimulator2(const SPMatrix* checkMatrix);

    //==== 復号実験 ====//
    DSOutput decodingSimulation(DSInput simulationParameter);

    //===== 雑音生成 ====//
    double getGaussianNoise(double average, double variance) const;
    DSAveVar checkAverageAndVariance(const vector<double>* vec) const;

    //==== 復号語と送信符号語が等しいかどうか ====//
    bool checkDecodeWord(const vector<BinaryFiniteField>* codeWord, const vector<BinaryFiniteField>* recepWord) const;

    //==== 進捗状況の表示 ====//
    void printProgress(double progressRate) const;

    //==== 拡張子を取り除く ====//
    string removeExtension(string str) const;

};

#endif
