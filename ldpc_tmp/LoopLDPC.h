/*
    Loop空間結合LDPC符号
    複数の帯を結合できるようにした空間結合符号。

    入力パラメータ：
*/

#if !defined(___Class_LoopLDPC)
#define ___Class_LoopLDPC

#include "Matrix.h"
#include <exception>
#include <vector>

class LoopLDPC {
public:
    //==== 入力データの構造 ====//
    typedef struct ___LBANDDATA {    // 空間結合プロトグラフ1つ分のデータ
        ___LBANDDATA() : l(3), r(6), L(7), conBandLeft(-1), conBandRight(-1), conPosLeft(0), conPosRight(0) {}
        int l;    // 列重み
        int r;    // 行重み
        int L;    // 結合数
        int conBandLeft;    // 左方向に接続する帯の番号(接続しない場合は-1を指定)
        int conPosLeft;    // 左方向に接続する帯の列位置
        int conBandRight;    // 右方向に接続する帯の番号(接続しない場合は-1を指定)
        int conPosRight;    // 右方向に接続する帯の列位置
    } LBandData;

    //==== エラークラス ====//
    class ParameterExeption {};

private:
    //==== クラス変数 ====//
    const vector<LBandData>& bandData;
    Matrix<BinaryFiniteField> mat;    // 行列本体
    vector<unsigned long> bsrp;    // 各帯の先頭位置(行)
    vector<unsigned long> bscp;    // 各帯の戦闘位置(列)
public:
    //==== コンストラクタ ====//
    LoopLDPC(const vector<LBandData>& input);

    //==== 公開メソッド ====//
    // パラメータの内容に従ってLoopLDPC符号を生成する。
    Matrix<BinaryFiniteField> generate();

    // 入力された帯の情報を公開する。
    void printBandData();

private:
    //==== 非公開メソッド ====//
    // 行列のサイズを決定する
    unsigned long calculateRowSize(const vector<LBandData>& input);
    unsigned long calculateColSize(const vector<LBandData>& input);

    // 各帯の行列中の開始位置(行と列)を取得し、bsrpとbscpを更新
    void calculateBandStartPos();

    // bandNum番目の帯におけるblockNum番目のブロックの先頭位置(行と列)を取得
    unsigned long calculateBlockStartRowPos(int bandNum, int blockNum);
    unsigned long calculateBlockStartColPos(int bandNum, int blockNum);

    // bandNum番目の帯におけるblockNum番目のブロックを生成する
    void generateBlock(int bandNum, int blockNum);

    // bandNum番目の帯を行列中に配置する
    void generateBand(int bandNum);

    // 全ての空間結合プロトグラフを検査行列中に配置する
    void generateBaseSCCProtograph();

    // 結合基準行の算出
    unsigned long calculateBaseConRowLeft(int bandNum);
    unsigned long calculateBaseConRowRight(int bandNum);
    // 結合基準列の算出
    unsigned long calculateBaseConColLeft(int bandNum);
    unsigned long calculateBaseConColRight(int bandNum);

    // 配置ブロックの単位長を算出
    unsigned long calculateConBlockLengthLeft(int bandNum);
    unsigned long calculateConBlockLengthRight(int bandNum);

    // stepステップ目における残り配置ビット数を算出
    int calculateRemainSetBits(int bandNum, int step);

    // 配置ブロック数を算出
    int calculateSetBlockNumLeft(int bandNum, int remainBits);
    int calculateSetBlockNumRight(int bandNum, int remainBits);

    // 結合基準行の更新
    unsigned long calculateBRowLeft(int bandNum, int step);
    unsigned long calculateBRowRight(int bandNum, int step);

    // 指定範囲に結合ブロックを配置
    void generateConnectBar(long bRow, long bColS, long bColE);

    // 結合の反復処理1回分
    int connectWithStepLeft(int bandNum, int step);
    int connectWithStepRight(int bandNum, int step);

    // 左右の結合部の作成
    void generateConnectParts(int bandNum);

    // 結合部全体の作成
    void generateConnections();
};

#endif
