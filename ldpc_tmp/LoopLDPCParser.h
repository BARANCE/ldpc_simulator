/*
    LoopLDPC符号のパラメータを設定ファイルから読み取るためのパーサ
*/

#if !defined(___Class_LoopLDPCParser)
#define ___Class_LoopLDPCParser

#include "LoopLDPC.h"
#include "HashOperator.h"

class LoopLDPCParser {
    vector<LoopLDPC::LBandData> bandData;
public:
    //==== コンストラクタ ====//
    LoopLDPCParser(const HashOperator& hash);

    //==== 公開メソッド ====//
    // Hashから設定を読み出し、vectorデータに変換する
    vector<LoopLDPC::LBandData> translate();
};

#endif
