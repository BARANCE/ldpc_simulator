#include "LoopLDPCParser.h"

LoopLDPCParser::LoopLDPCParser(const HashOperator& hash) : bandData(0) {
    bool check = hash.tryKey("copyNum");
    if(check) {
        int size = hash.getValueByKey<int>("copyNum");    // 帯の個数
        for(int i = 0; i < size; i++) {
            LoopLDPC::LBandData nowRead;
            stringstream ssl, ssr, ssL, ssbl, ssbr, sscl, sscr;
            ssl << "l" << i;
            ssr << "r" << i;
            ssL << "L" << i;
            ssbl << "ssbl" << i;
            ssbr << "ssbr" << i;
            sscl << "sscl" << i;
            sscr << "sscr" << i;
            if(hash.tryKey(ssl.str())) nowRead.l = hash.getValueByKey<int>(ssl.str());
            if(hash.tryKey(ssr.str())) nowRead.r = hash.getValueByKey<int>(ssr.str());
            if(hash.tryKey(ssL.str())) nowRead.L = hash.getValueByKey<int>(ssL.str());
            if(hash.tryKey(ssbl.str())) nowRead.conBandLeft = hash.getValueByKey<int>(ssbl.str());
            if(hash.tryKey(ssbr.str())) nowRead.conBandRight = hash.getValueByKey<int>(ssbr.str());
            if(hash.tryKey(sscl.str())) nowRead.conPosLeft = hash.getValueByKey<int>(sscl.str());
            if(hash.tryKey(sscr.str())) nowRead.conPosRight = hash.getValueByKey<int>(sscr.str());

            bandData.push_back(nowRead);
        }
    }
}

// Hashから設定を読み出し、vectorデータに変換する
vector<LoopLDPC::LBandData> LoopLDPCParser::translate() {
    return bandData;
}
