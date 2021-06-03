#include "LDPCDecoder2.h"
#include <math.h>

//==== コンストラクタ ====//
LDPCDecoder2::LDPCDecoder2(const SPMatrix* checkMatrix) : H(checkMatrix), mb(checkMatrix), channel(C_AWGN), SNR(0.0), fp(0.0) {
    //ESPInput espi = {vector<double>(4000, -1.0), 100};
    //sumProductDecode(espi);
}

//==== 通信路の選択 ====//
void LDPCDecoder2::setChannel(EChannel channelName) {
    channel = channelName;
}

//==== 通信路に関するパラメータの設定 ====//
void LDPCDecoder2::setSNR(double rate) {
    SNR = rate;
}
void LDPCDecoder2::setFlipProbability(double flipProbability) {
    fp = flipProbability;
}

//==== 復号 ====//
// 硬判定
LDPCDecoder2::ESPOutput LDPCDecoder2::hardDecision(ESPInput decodeParameter) {
    //==== 変数宣言・前処理 ====//
    unsigned long hc = H->col();
    const vector<double> receptionWord = decodeParameter.receptionWord;    // 受信語
    vector<BinaryFiniteField> nHat(hc);    // 推定語

    for(int i = 0; i < hc; i++) {
        if(receptionWord[i] > 0.0) {
            nHat[i] = 0;
        } else if(receptionWord[i] < 0.0) {
            nHat[i] = 1;
        } else {    // タイブレーク
            double rnd = genrand_real2();
            if(rnd < 0.5) nHat[i] = 0;
            else nHat[i] = 1;
        }
    }

    //==== パリティ検査 ====//
    bool failFrag = true;    // 復号が失敗したらtrue
    if(parityCheck(&nHat)) {
        failFrag = false;
    }

    //==== 実験データの整形 ====//
    ESPOutput retObj = {vector<BinaryFiniteField>(hc), false, 0, DECODE_SUCCESS};    // 返却構造体
    for(int i = 0; i < hc; i++) {
        retObj.decodeWord[i] = nHat[i];
    }
    retObj.convergenceNumber = 0;
    if(failFrag) {
        retObj.failFrag = true;
        retObj.errorType = OVER_MAX_CONVERGENCE;    // 訂正失敗
    } else {
        retObj.failFrag = false;
        retObj.errorType = DECODE_SUCCESS;    // 訂正成功(復号成功ではない)
    }

    return retObj;
}

// sum-product復号
LDPCDecoder2::ESPOutput LDPCDecoder2::sumProductDecode(ESPInput decodeParameter) {
    //==== 変数宣言・前処理 ====//
    unsigned long hr = H->row();    // チェックノード数
    unsigned long hc = H->col();    // 変数ノード数
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;

    unsigned long loop = 0;    // 反復回数
    const unsigned long maxLoop = decodeParameter.convergenceFailNumber;    // 最大反復回数(収束失敗条件)
    unsigned long conv = 0;    // 収束回数
    const vector<double> receptionWord = decodeParameter.receptionWord;    // 受信語

    const double codeRate = (double)(hc - hr) / hc;    // 符号化率
    const double variance = 0.5 * (1.0 / pow(10.0, SNR / 10.0)) * (double)hc / (double)(hc - hr);    // SN比から求めた分散値
    vector<double> logLikeRate(hc);    // 変数ノードごとの対数尤度比

    vector<double> renewChecker(hc, -100.0);    // 更新が進んでいるかを調査するベクトル
    bool badConvergenceFrag = false;

    vector<BinaryFiniteField> nHat(hc);    // 一時推定語(パリティ検査を通過したら返却構造体に格納)

    ESPOutput retObj = {vector<BinaryFiniteField>(hc), false, 0, DECODE_SUCCESS};    // 返却構造体

    //==== 全変数ノードの対数尤度比を計算 ====//
    for(int i = 0; i < hc; i++) {
        logLikeRate[i] = 2.0 * receptionWord[i] / variance;
    }

    //==== 全てのチェックノードから隣接する変数ノードに0を送る ====//
    for(int i = 0; i < hr; i++) {
        const vector<int> variableList = (*matRow)[i];
        int size = variableList.size();
        for(int j = 0; j < size; j++) {
            int variableId = variableList[j] - 1;
            mb.checkToVariable(i, variableId, 0.0);
        }
    }

    for(loop = 0; loop < maxLoop; loop++) {
        //cout << "loop = " << loop << "\n";

        //cout << "====変数ノード処理====\n";

        //==== 変数ノード処理 & 一時推定ビットの決定 ====//
        unsigned long renewCounter = 0;
        for(int i = 0; i < hc; i++) {
            double llr = logLikeRate[i];    // 対数尤度比
            //double vSum = sumMessageFromCheckNode(i);    // メッセージの和
            //double sendTmp = llr + vSum;

            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            double vSum = llr;
            for(int j = 0; j < size; j++) {
                vSum += mb.variableBuf[i][j];
            }

            for(int j = 0; j < size; j++) {    // チェックノードへのメッセージを算出する
                unsigned long checkId = (*checkList)[j] - 1;
                double writeValue = vSum - mb.getVariableMessage(i, checkId);
                //cout << "x" << i << "→ c" << checkId << " : " << writeValue << "\n";
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << sendTmp << " ";

            //==== 推定ビットの決定 ====//
            if(vSum < 0.0) {
                nHat[i] = 1;
            } else {
                nHat[i] = 0;
            }

            //==== 更新が進んでいるか ====//
            if(renewChecker[i] == vSum) {
                renewCounter++;
            }
            renewChecker[i] = vSum;
        }

        //==== パリティ検査 ====//
        /*for(int i = 0; i < hc; i++) {
            cout << nHat[i] << " ";
        }
        cout << "\n";*/
        if(parityCheck(&nHat)) {
        //if(isZero(&nHat)) {    // 全ゼロ符号語にしか通用しない
            //cout << "復号成功!!\n";
            break;
        }
        /*if(renewCounter >= hr) {    // 誤収束
            badConvergenceFrag = true;
            break;
        }*/
        // 収束失敗はメインループのループ回数で判定

        //cout << "====チェックノード処理====\n";

        //==== チェックノード処理 ====//
        for(int i = 0; i < hr; i++) {
            const vector<int>* variableList = &((*matRow)[i]);
            unsigned long size = variableList->size();

            double cSum = 0.0;
            double cProd = 1.0;
            for(int j = 0; j < size; j++) {
                double tmp1 = mb.checkBuf[i][j];
                cProd *= sign(tmp1);
                cSum += gallagerf(absolute(tmp1));
            }

            for(int j = 0; j < size; j++) {    // 変数ノードへのメッセージを算出する
                unsigned long variableId = (*variableList)[j] - 1;    // 送信先変数ノードのID

                //cout << "weiteValue = " << writeValue << ", writeValue2 = " << writeValue2 << "\n";
                double tmp2 = mb.getCheckMessage(i, variableId);    // 送信先からのデータ
                double signValue = sign(tmp2);
                if(tmp2 != 0.0) {
                    signValue = cProd / signValue;
                } else {
                    signValue = cProd;
                }

                double fValue = gallagerf(absolute(tmp2));
                fValue = gallagerf(cSum - fValue);

                double writeValue = signValue * fValue;

                //cout << "c" << i << "→ x" << variableId << " : " << writeValue << "\n";
                mb.checkToVariable(i, variableId, writeValue);
            }
        }



        //break;

    }
    if(loop < maxLoop) conv = loop + 1;    // 収束回数
    else conv = loop;

    //==== 実験データの整形 ====//
    for(int i = 0; i < hc; i++) {
        retObj.decodeWord[i] = nHat[i];
    }
    retObj.convergenceNumber = conv;
    if(conv >= maxLoop) {
        retObj.failFrag = true;
        retObj.errorType = OVER_MAX_CONVERGENCE;    // 収束失敗
    } else {
        if(badConvergenceFrag) {
            retObj.failFrag = true;
            retObj.errorType = BAD_CONVERGENCE;    // 誤収束
        } else {
            retObj.failFrag = false;
            retObj.errorType = DECODE_SUCCESS;    // 訂正成功(復号成功ではない)
        }
    }

    //printVector(&retObj.decodeWord);

    return retObj;
}

// min-sum復号
LDPCDecoder2::ESPOutput LDPCDecoder2::minSumDecode(ESPInput decodeParameter, double scale) {
    //==== 変数宣言・前処理 ====//
    unsigned long hr = H->row();    // チェックノード数
    unsigned long hc = H->col();    // 変数ノード数
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;

    unsigned long loop = 0;    // 反復回数
    const unsigned long maxLoop = decodeParameter.convergenceFailNumber;    // 最大反復回数(収束失敗条件)
    unsigned long conv = 0;    // 収束回数
    const vector<double> receptionWord = decodeParameter.receptionWord;    // 受信語

    const double codeRate = (double)(hc - hr) / hc;    // 符号化率
    const double variance = 0.5 * (1.0 / pow(10.0, SNR / 10.0)) * (double)hc / (double)(hc - hr);    // SN比から求めた分散値
    vector<double> logLikeRate(hc);    // 変数ノードごとの対数尤度比

    vector<double> renewChecker(hc, -100.0);    // 更新が進んでいるかを調査するベクトル
    bool badConvergenceFrag = false;

    vector<BinaryFiniteField> nHat(hc);    // 一時推定語(パリティ検査を通過したら返却構造体に格納)

    ESPOutput retObj = {vector<BinaryFiniteField>(hc), false, 0, DECODE_SUCCESS};    // 返却構造体

    //----テスト用ここから----
    //Matrix<double> hatCheck(maxLoop, hc);    // 推定ビットの変遷を記録する行列
    //----テスト用ここまで----

    //==== 全変数ノードの対数尤度比を計算 ====//
    for(int i = 0; i < hc; i++) {
        logLikeRate[i] = 2.0 * receptionWord[i] / variance;
    }

    //==== 全てのチェックノードから隣接する変数ノードに0を送る ====//
    for(int i = 0; i < hr; i++) {
        const vector<int> variableList = (*matRow)[i];
        int size = variableList.size();
        for(int j = 0; j < size; j++) {
            int variableId = variableList[j] - 1;
            mb.checkToVariable(i, variableId, 0.0);
        }
    }

    for(loop = 0; loop < maxLoop; loop++) {
        //cout << "loop = " << loop << "\n";

        //cout << "====変数ノード処理====\n";

        //==== 変数ノード処理 & 一時推定ビットの決定 ====//
        unsigned long renewCounter = 0;
        for(int i = 0; i < hc; i++) {
            double llr = logLikeRate[i];    // 対数尤度比
            //double vSum = sumMessageFromCheckNode(i);    // メッセージの和
            //double sendTmp = llr + vSum;

            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            double vSum = llr;
            for(int j = 0; j < size; j++) {
                vSum += mb.variableBuf[i][j];
            }

            for(int j = 0; j < size; j++) {    // チェックノードへのメッセージを算出する
                unsigned long checkId = (*checkList)[j] - 1;
                double writeValue = vSum - mb.getVariableMessage(i, checkId);
                //cout << "x" << i << "→ c" << checkId << " : " << writeValue << "\n";
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << sendTmp << " ";

            //---- テスト用ここから ----//
            //hatCheck[loop][i] = vSum;
            //---- テスト用ここまで ----//
            if(vSum != vSum) {
                cout << "NaNが発生しています。\n";
            }

            //==== 推定ビットの決定 ====//
            if(vSum < 0.0) {
                nHat[i] = 1;
            } else if(vSum > 0.0) {
                nHat[i] = 0;
            } else {
                double rnd = genrand_real2();
                if(rnd < 0.5) nHat[i] = 0;
                else nHat[i] = 1;
            }

            //==== 更新が進んでいるか ====//
            if(renewChecker[i] == vSum) {
                renewCounter++;
            }
            renewChecker[i] = vSum;
        }

        //==== パリティ検査 ====//
        /*for(int i = 0; i < hc; i++) {
            cout << nHat[i] << " ";
        }
        cout << "\n";*/
        if(parityCheck(&nHat)) {
        //if(isZero(&nHat)) {    // 全ゼロ符号語にしか通用しない
            //cout << "復号成功!!\n";
            break;
        }
        /*if(renewCounter >= hr) {    // 誤収束
            badConvergenceFrag = true;
            break;
        }*/
        // 収束失敗はメインループのループ回数で判定

        //cout << "====チェックノード処理====\n";

        //==== チェックノード処理 ====//
        for(int i = 0; i < hr; i++) {
            const vector<int>* variableList = &((*matRow)[i]);
            unsigned long size = variableList->size();

            //double cMin = absolute(mb.checkBuf[i][0]);
            double cProd;
            if(size != 0) {
                cProd = sign(mb.checkBuf[i][0]);
            } else {
                cProd = 0.0;
            }
            for(int j = 1; j < size; j++) {
                double tmp1 = mb.checkBuf[i][j];
                cProd *= sign(tmp1);
                //if(absolute(tmp1) < cMin) {
                //    cMin = absolute(tmp1);
                //}
            }

            for(int j = 0; j < size; j++) {    // 変数ノードへのメッセージを算出する
                unsigned long variableId = (*variableList)[j] - 1;    // 送信先変数ノードのID

                //cout << "weiteValue = " << writeValue << ", writeValue2 = " << writeValue2 << "\n";
                double tmp2 = mb.getCheckMessage(i, variableId);    // 送信先からのデータ
                double signValue = sign(tmp2);
                if(tmp2 != 0.0) {
                    signValue = cProd / signValue;
                } else {
                    signValue = cProd;
                }

                //double fValue = gallagerf(absolute(tmp2));
                //fValue = gallagerf(cSum - fValue);
                double cMin = 0.0;
                if(size > 0 && (*variableList)[0] - 1 == variableId) {
                    cMin = absolute(mb.checkBuf[i][1]);
                } else if(size > 0) {
                    cMin = absolute(mb.checkBuf[i][0]);
                }
                for(int k = 1; k < size; k++) {
                    if((*variableList)[k] - 1 != variableId) {    // 送る相手を除く
                        if(absolute(mb.checkBuf[i][k]) < cMin) {
                            cMin = absolute(mb.checkBuf[i][k]);
                        }
                    }
                }

                double writeValue = scale * signValue * cMin;

                //cout << "c" << i << "→ x" << variableId << " : " << writeValue << "\n";
                mb.checkToVariable(i, variableId, writeValue);
            }
        }



        //break;

    }
    if(loop < maxLoop) conv = loop + 1;    // 収束回数
    else conv = loop;

    //==== 実験データの整形 ====//
    for(int i = 0; i < hc; i++) {
        retObj.decodeWord[i] = nHat[i];
    }
    retObj.convergenceNumber = conv;
    if(conv >= maxLoop) {
        retObj.failFrag = true;
        retObj.errorType = OVER_MAX_CONVERGENCE;    // 収束失敗
    } else {
        if(badConvergenceFrag) {
            retObj.failFrag = true;
            retObj.errorType = BAD_CONVERGENCE;    // 誤収束
        } else {
            retObj.failFrag = false;
            retObj.errorType = DECODE_SUCCESS;    // 訂正成功(復号成功ではない)
        }
    }

    //---- テスト用ここから ----//
    //hatCheck.exportMatrix("hatCheck_RPEGSCCm400_n4000_wc4_b40.txt");
    //---- テスト用ここまで ----//


    //printVector(&retObj.decodeWord);

    return retObj;
}

// 消失訂正sum-product復号法
LDPCDecoder2::ERSPOutput LDPCDecoder2::sumProductDecodeOnBEC(ESPInput decodeParameter, ECDList researchParameter) {
    //==== 変数宣言・前処理 ====//
    unsigned long hr = H->row();    // チェックノード数
    unsigned long hc = H->col();    // 変数ノード数
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;

    unsigned long loop = 0;    // 反復回数
    const unsigned long maxLoop = decodeParameter.convergenceFailNumber;    // 最大反復回数(収束失敗条件)
    unsigned long conv = 0;    // 収束回数
    const vector<double> receptionWord = decodeParameter.receptionWord;    // 受信語

    bool badConvergenceFrag = false;

    vector<double> nHat(hc);    // 一時推定語(パリティ検査を通過したら返却構造体に格納)

    int convCheckFrag = researchParameter.convCheckFrag;    // 収束チェックを実行するかどうか
    int convCheckSize = researchParameter.convCheckSize;    // 収束チェックを行うためのバッファのサイズ
    vector<vector<double> > convCheck(convCheckSize, vector<double>(hc, 0.0));    // 収束状況をチェックする2次元ベクトル

    ERSPOutput retObj = {vector<double>(hc), false, 0, DECODE_SUCCESS, vector<vector<double> >(0)};    // 返却構造体

    //==== 初期推定値は受信語 ====//
    for(int i = 0; i < hc; i++) {
        nHat[i] = receptionWord[i];
        if(convCheckFrag) {
            convCheck[0][i] = receptionWord[i];
        }
    }

    //==== 全てのチェックノードから隣接する変数ノードにeを送る ====//
    for(int i = 0; i < hr; i++) {
        const vector<int> variableList = (*matRow)[i];
        int size = variableList.size();
        for(int j = 0; j < size; j++) {
            int variableId = variableList[j] - 1;
            mb.checkToVariable(i, variableId, ERASURE);
        }
    }

    for(loop = 0; loop < maxLoop; loop++) {

        //==== 変数ノード処理 & 一時推定ビットの決定 ====//
        unsigned long renewCounter = 0;
        for(int i = 0; i < hc; i++) {
            const vector<int>* checkList = &((*matCol)[i]);    // i番目の変数ノードに接続しているチェックノード番号のリスト
            unsigned long size = checkList->size();    // i番目の列の変数ノードの次数(列重み)

            if(nHat[i] == ONE || nHat[i] == ZERO) {    // 変数ノードの推定値が0または1の場合
                // 隣接する全てのチェックノードに推定値を送信する
                for(int j = 0; j < size; j++) {
                    unsigned long checkId = (*checkList)[j] - 1;
                    double writeValue = nHat[i];
                    mb.variableToCheck(i, checkId, writeValue);
                }

            } else if(nHat[i] == ERASURE) {    // 変数ノードの推定値がeの場合
                // 送られてきたチェックノード→変数ノードのメッセージに値があるかをカウントする
                int valueCounter = 0;    // 0または1のメッセージの数をカウントする
                long messagePointer = -1;    // 0または1のメッセージが送られてきたメッセージの送り相手(無ければ-1)
                double value = ERASURE;    // 推定値

                for(int j = 0; j < size; j++) {
                    unsigned long checkId = (*checkList)[j] - 1;
                    double message = mb.getVariableMessage(i, checkId);
                    if(message == ZERO || message == ONE) {
                        valueCounter++;
                        messagePointer = checkId;
                        value = message;
                        if(valueCounter >= 2) break;    // これ以上やる必要なし
                    } else if(message == ERASURE) {
                        // 何もしない
                    } else {
                        cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]メッセージ値に、0・1・eの3種類以外の値が使用されています。\n";
                    }
                }

                // カウンターの値によって条件分岐
                if(valueCounter == 0) {    // 全てERASURE
                    // 送信メッセージは全てERASURE
                    for(int j = 0; j < size; j++) {
                        unsigned long checkId = (*checkList)[j] - 1;
                        double writeValue = ERASURE;
                        mb.variableToCheck(i, checkId, writeValue);
                    }
                } else if(valueCounter >= 2) {
                    // 送信メッセージは推定値
                    for(int j = 0; j < size; j++) {
                        unsigned long checkId = (*checkList)[j] - 1;
                        double writeValue = value;
                        mb.variableToCheck(i, checkId, writeValue);
                    }
                } else if(valueCounter == 1) {
                    // 0,1の値を送ってきたチェックノードにはERASUREを、それ以外のERASUREを送ってきたチェックノードには推定値を送る
                    if(messagePointer == -1) cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]消失シンボル以外の値を送ってきたチェックノードの値が記録されていません。\n";
                    for(int j = 0; j < size; j++) {
                        unsigned long checkId = (*checkList)[j] - 1;
                        double writeValue = value;
                        if(checkId == messagePointer) {
                            writeValue = ERASURE;
                        }
                        mb.variableToCheck(i, checkId, writeValue);
                    }
                } else {
                    cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]カウンターの値が異常です。\n";
                    throw this;
                }

                //==== 推定値の更新 ====//
                if(value != ERASURE) nHat[i] = value;
            } else {
                cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]推定値に、0・1・eの3種類以外の値が使用されています。\n";
                throw this;
            }
        }

        //==== パリティ検査 ====//
        int parityCounter = 0;
        for(int i = 0; i <hc; i++) {
            if(convCheckFrag && loop < convCheckSize - 1) {
                convCheck[loop + 1][i] = nHat[i];
                if(convCheck[loop][i] != ERASURE && convCheck[loop + 1][i] == ERASURE) {
                    cout << "[loop : " << loop << ", i = " << i << "] おかしいぜ\n";
                }
            }
            if(nHat[i] == ERASURE) {
                parityCounter++;
                //break;
            }
        }
        if(parityCounter == 0) break;    // 復号完了

        // デバッグ用
        /*if(loop == 0) {
        for(int i = 0; i < hc; i++) {
            if(nHat[i] == ZERO) cout << "0";
            else if(nHat[i] == ONE) cout << "1";
            else if(nHat[i] == ERASURE) cout << "e";
            else cout << "?";
        }
        cout << "\n";
        }*/

        //==== チェックノード処理 ====//
        for(int i = 0; i < hr; i++) {
            const vector<int>* variableList = &((*matRow)[i]);
            unsigned long size = variableList->size();

            // 送られてきた変数ノード→チェックノードのメッセージに消失シンボルが何個存在するかをカウントする
            int erasureCounter = 0;
            long messagePointer = -1;    // erasureを送ってきた変数ノードの位置(無ければ-1)
            double cSum = ZERO;    // パリティチェックの合計値(メッセージ算出に使用)

            for(int j = 0; j < size; j++) {
                unsigned long variableId = (*variableList)[j] - 1;
                double message = mb.getCheckMessage(i, variableId);
                if(message == ZERO || message == ONE) {
                    if(message == ONE) {
                        if(cSum == ZERO) cSum = ONE;
                        else if(cSum == ONE) cSum = ZERO;
                        else cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]メッセージ合計値に、0・1の2種類以外の値が使用されています。\n";
                    }
                } else if(message == ERASURE) {
                    erasureCounter++;
                    messagePointer = variableId;
                } else {
                    cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]メッセージ値に、0・1・eの3種類以外の値が使用されています。\n";
                }
            }

            // カウンターの値によって条件分岐
            if(erasureCounter == 0) {    // 全て0か1
                // 送信メッセージは、送信相手からのメッセージを除いた受信メッセージの和
                for(int j = 0; j < size; j++) {
                    unsigned long variableId = (*variableList)[j] - 1;
                    double writeValue = ZERO;
                    if(mb.getCheckMessage(i, variableId) == ONE) {    // writeValue = cSum + mb.getCheckMessage(i, variableId)をやろうとしている部分
                        if(cSum == ONE) writeValue = ZERO;
                        else if(cSum == ZERO) writeValue = ONE;
                        else cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]メッセージ合計値に、0・1の2種類以外の値が使用されています。\n";
                    } else if(mb.getCheckMessage(i, variableId) == ZERO) {
                        if(cSum == ONE) writeValue = ONE;
                        else if(cSum == ZERO) writeValue = ZERO;
                        else cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]メッセージ合計値に、0・1の2種類以外の値が使用されています。\n";
                    } else {
                        cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]ここではメッセージがeであることを想定していません。\n";
                    }
                    mb.checkToVariable(i, variableId, writeValue);
                }
            } else if(erasureCounter >= 2) {
                // 送信メッセージは全てERASURE
                for(int j = 0; j < size; j++) {
                    unsigned long variableId = (*variableList)[j] - 1;
                    double writeValue = ERASURE;
                    mb.checkToVariable(i, variableId, writeValue);
                }
            } else if(erasureCounter == 1) {
                // ERASUREを送ってきた変数ノードには受信メッセージの和を、そうでない変数ノードにはERASUREを送る
                for(int j = 0; j < size; j++) {
                    unsigned long variableId = (*variableList)[j] - 1;
                    double writeValue = ERASURE;
                    if(variableId == messagePointer) {
                        writeValue = cSum;
                    }
                    mb.checkToVariable(i, variableId, writeValue);
                }
            } else {
                cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]カウンターの値が異常です。\n";
            }

        }


    }

    if(loop < maxLoop) conv = loop + 1;    // 収束回数
    else conv = loop;

    //==== 実験データの整形 ====//
    for(int i = 0; i < hc; i++) {
        retObj.decodeWord[i] = nHat[i];
    }
    retObj.convergenceNumber = conv;
    if(conv >= maxLoop) {
        retObj.failFrag = true;
        retObj.errorType = OVER_MAX_CONVERGENCE;    // 収束失敗
    } else {
        if(badConvergenceFrag) {
            retObj.failFrag = true;
            retObj.errorType = BAD_CONVERGENCE;    // 誤収束
        } else {
            retObj.failFrag = false;
            retObj.errorType = DECODE_SUCCESS;    // 訂正成功(復号成功ではない)
        }
    }
    retObj.convCheckVec = convCheck;

    return retObj;
}

// 修正min-sum復号
// assumptionVector : ショートンした箇所の値を格納するベクトル
// assumptionVector[i] > 0 : i番目の変数ノードは0を仮定
// assumptionVector[i] < 0 : i番目の変数ノードは1を仮定
// assumptionVector[i] = 0 : ショートンしない変数ノード
LDPCDecoder2::ESPOutput LDPCDecoder2::modMinSumDecode(ESPInput decodeParameter, double scale, vector<int>* assumptionVector) {
    //==== 変数宣言・前処理 ====//
    unsigned long hr = H->row();    // チェックノード数
    unsigned long hc = H->col();    // 変数ノード数
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;

    unsigned long loop = 0;    // 反復回数
    const unsigned long maxLoop = decodeParameter.convergenceFailNumber;    // 最大反復回数(収束失敗条件)
    unsigned long conv = 0;    // 収束回数
    const vector<double> receptionWord = decodeParameter.receptionWord;    // 受信語

    const double codeRate = (double)(hc - hr) / hc;    // 符号化率
    const double variance = 0.5 * (1.0 / pow(10.0, SNR / 10.0)) * (double)hc / (double)(hc - hr);    // SN比から求めた分散値
    vector<double> logLikeRate(hc);    // 変数ノードごとの対数尤度比

    vector<double> renewChecker(hc, -100.0);    // 更新が進んでいるかを調査するベクトル
    bool badConvergenceFrag = false;

    vector<BinaryFiniteField> nHat(hc);    // 一時推定語(パリティ検査を通過したら返却構造体に格納)

    ESPOutput retObj = {vector<BinaryFiniteField>(hc), false, 0, DECODE_SUCCESS};    // 返却構造体

    //----テスト用ここから----
    //Matrix<double> hatCheck(maxLoop, hc);    // 推定ビットの変遷を記録する行列
    //----テスト用ここまで----

    //==== 全変数ノードの対数尤度比を計算 ====//
    const double A_SCALE = 10000000;
    for(int i = 0; i < hc; i++) {
        if((*assumptionVector)[i] > 0) {
            logLikeRate[i] = A_SCALE;
        } else if ((*assumptionVector)[i] < 0) {
            logLikeRate[i] = -1.0 * A_SCALE;
        } else {
            logLikeRate[i] = 2.0 * receptionWord[i] / variance;
        }
    }
    //==== 全てのチェックノードから隣接する変数ノードに0を送る ====//
    for(int i = 0; i < hr; i++) {
        const vector<int> variableList = (*matRow)[i];
        int size = variableList.size();
        for(int j = 0; j < size; j++) {
            int variableId = variableList[j] - 1;
            mb.checkToVariable(i, variableId, 0.0);
        }
    }

    for(loop = 0; loop < maxLoop; loop++) {
        //cout << "loop = " << loop << "\n";

        //cout << "====変数ノード処理====\n";

        //==== 変数ノード処理 & 一時推定ビットの決定 ====//
        unsigned long renewCounter = 0;
        for(int i = 0; i < hc; i++) {
            double llr = logLikeRate[i];    // 対数尤度比
            //double vSum = sumMessageFromCheckNode(i);    // メッセージの和
            //double sendTmp = llr + vSum;

            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            double vSum = llr;
            for(int j = 0; j < size; j++) {
                vSum += mb.variableBuf[i][j];
            }

            for(int j = 0; j < size; j++) {    // チェックノードへのメッセージを算出する
                unsigned long checkId = (*checkList)[j] - 1;
                double writeValue = vSum - mb.getVariableMessage(i, checkId);
                //cout << "x" << i << "→ c" << checkId << " : " << writeValue << "\n";
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << sendTmp << " ";

            //---- テスト用ここから ----//
            //hatCheck[loop][i] = vSum;
            //---- テスト用ここまで ----//

            //==== 推定ビットの決定 ====//
            if(vSum < 0.0) {
                nHat[i] = 1;
            } else if(vSum > 0.0) {
                nHat[i] = 0;
            } else {
                double rnd = genrand_real2();
                if(rnd < 0.5) nHat[i] = 0;
                else nHat[i] = 1;
            }

            //==== 更新が進んでいるか ====//
            if(renewChecker[i] == vSum) {
                renewCounter++;
            }
            renewChecker[i] = vSum;
        }

        //==== パリティ検査 ====//
        /*for(int i = 0; i < hc; i++) {
            cout << nHat[i] << " ";
        }
        cout << "\n";*/
        if(parityCheck(&nHat)) {
        //if(isZero(&nHat)) {    // 全ゼロ符号語にしか通用しない
            //cout << "復号成功!!\n";
            break;
        }
        /*if(renewCounter >= hr) {    // 誤収束
            badConvergenceFrag = true;
            break;
        }*/
        // 収束失敗はメインループのループ回数で判定

        //cout << "====チェックノード処理====\n";

        //==== チェックノード処理 ====//
        for(int i = 0; i < hr; i++) {
            const vector<int>* variableList = &((*matRow)[i]);
            unsigned long size = variableList->size();

            //double cMin = absolute(mb.checkBuf[i][0]);
            double cProd = sign(mb.checkBuf[i][0]);
            for(int j = 1; j < size; j++) {
                double tmp1 = mb.checkBuf[i][j];
                cProd *= sign(tmp1);
                //if(absolute(tmp1) < cMin) {
                //    cMin = absolute(tmp1);
                //}
            }

            for(int j = 0; j < size; j++) {    // 変数ノードへのメッセージを算出する
                unsigned long variableId = (*variableList)[j] - 1;    // 送信先変数ノードのID

                //cout << "weiteValue = " << writeValue << ", writeValue2 = " << writeValue2 << "\n";
                double tmp2 = mb.getCheckMessage(i, variableId);    // 送信先からのデータ
                double signValue = sign(tmp2);
                if(tmp2 != 0.0) {
                    signValue = cProd / signValue;
                } else {
                    signValue = cProd;
                }

                //double fValue = gallagerf(absolute(tmp2));
                //fValue = gallagerf(cSum - fValue);
                double cMin = 0.0;
                if(size > 0 && (*variableList)[0] - 1 == variableId) {
                    cMin = absolute(mb.checkBuf[i][1]);
                } else if(size > 0) {
                    cMin = absolute(mb.checkBuf[i][0]);
                }
                for(int k = 1; k < size; k++) {
                    if((*variableList)[k] - 1 != variableId) {    // 送る相手を除く
                        if(absolute(mb.checkBuf[i][k]) < cMin) {
                            cMin = absolute(mb.checkBuf[i][k]);
                        }
                    }
                }

                double writeValue = scale * signValue * cMin;

                //cout << "c" << i << "→ x" << variableId << " : " << writeValue << "\n";
                mb.checkToVariable(i, variableId, writeValue);
            }
        }



        //break;

    }
    if(loop < maxLoop) conv = loop + 1;    // 収束回数
    else conv = loop;

    //==== 実験データの整形 ====//
    for(int i = 0; i < hc; i++) {
        retObj.decodeWord[i] = nHat[i];
    }
    retObj.convergenceNumber = conv;
    if(conv >= maxLoop) {
        retObj.failFrag = true;
        retObj.errorType = OVER_MAX_CONVERGENCE;    // 収束失敗
    } else {
        if(badConvergenceFrag) {
            retObj.failFrag = true;
            retObj.errorType = BAD_CONVERGENCE;    // 誤収束
        } else {
            retObj.failFrag = false;
            retObj.errorType = DECODE_SUCCESS;    // 訂正成功(復号成功ではない)
        }
    }

    //---- テスト用ここから ----//
    //hatCheck.exportMatrix("hoge.txt");
    //---- テスト用ここまで ----//

    //printVector(&retObj.decodeWord);

    return retObj;
}

LDPCDecoder2::ESPOutput LDPCDecoder2::iterativeMinSumDecode(ESPInput decodeParameter, double scale, const vector<int>* shortenIds) {
    unsigned long hr = H->row();
    unsigned long hc = H->col();
    vector<int> assumptionVector(hc);
    for(int i = 0; i < hc; i++) {
        if(shortenIds > 0) {
            assumptionVector[i] = 1;
        } else {
            assumptionVector[i] = 0;
        }
    }
    return modMinSumDecode(decodeParameter, scale, &assumptionVector);
}

LDPCDecoder2::ESPOutput LDPCDecoder2::shortenMinSumDecode(ESPInput decodeParameter, double scale, double shortenRate) {
    if(shortenRate < 0) shortenRate = 0.0;
    else if(shortenRate > 1) shortenRate = 1.0;
    unsigned long hr = H->row();
    unsigned long hc = H->col();
    ESPOutput result;
    vector<int> shortenIds(0);
    for(int i = 0; i < (unsigned long)(hc * shortenRate) / 2; i++) {
        shortenIds.push_back(i);
    }
    for(int i = hc - ((unsigned long)(hc * shortenRate) / 2); i < hc; i++) {
        shortenIds.push_back(i);
    }

    CombinationArray ca(&shortenIds);
    bool firstFrag = true;
    while(1) {
        vector<int> gv = ca.next();
        if(gv.empty() && firstFrag == false) {break;}
        else firstFrag = false;
        vector<int> svec(hc, 0);
        for(int i = 0; i < (unsigned long)(hc * shortenRate) / 2; i++) {
            svec[i] = 1;
        }
        for(int i = hc - ((unsigned long)(hc * shortenRate) / 2); i < hc; i++) {
            svec[i] = 1;
        }
        for(int i = 0; i < gv.size(); i++) {
            svec[gv[i]] = -1;
        }
        /*cout << "STAGE : " <<ca.getCount() << "\n";

        for(int i = 0; i < gv.size(); i++) {
            cout << gv[i] << " ";
        }
        cout << "\n";*/

        result = modMinSumDecode(decodeParameter, scale, &svec);
        if(result.failFrag == false) break;
    }

    return result;
}

//==== メッセージの一括計算 ====//
double LDPCDecoder2::sumMessageFromVariableNode(unsigned long checkId) const {
    const vector<double>* vec = &mb.checkBuf[checkId];
    unsigned long size = vec->size();
    double retVal = 0.0;
    for(int i = 0; i < size; i++) {
        retVal += (*vec)[i];
    }
    return retVal;
}
double LDPCDecoder2::sumMessageFromCheckNode(unsigned long variableId) const {
    const vector<double>* vec = &mb.variableBuf[variableId];
    unsigned long size = vec->size();
    double retVal = 0.0;
    for(int i = 0; i < size; i++) {
        retVal += (*vec)[i];
    }
    return retVal;
}

//==== アークハイパボリックタンジェント ====//
double LDPCDecoder2::atanh(double x) const {
    double retVal = 0.5;
    retVal *= log((1 + x) / (1 - x));
    return retVal;
}

//==== abs関数 ====//
double LDPCDecoder2::absolute(double x) const {
    if(x >= 0) return x;
    else return -x;
}

//==== sign関数 ====//
double LDPCDecoder2::sign(double x) const {
    if(x > 0) return 1.0;
    else if(x < 0) return -1.0;
    else return 0.0;
}

//==== gallagerのf関数 ====//
double LDPCDecoder2::gallagerf(double x) const {
    //if(x > 0.00001) {
        return log( (exp(x) + 1.0) / (exp(x) - 1.0) );
    //} else {
    //    return 12.21;
    //}
}

//==== パリティ検査 ====//
bool LDPCDecoder2::parityCheck(const vector<BinaryFiniteField>* nHat) const {    // 一般用。nHatが符号語であればtrueを返す
    unsigned long hr = H->row();
    unsigned long hc = H->col();
    const vector<vector<int> >* matRow = H->matRow;
    for(int i = 0; i < hr; i++) {    // それぞれの行に対して実行
        BinaryFiniteField bff = 0;
        const vector<int>* variableList = &((*matRow)[i]);
        unsigned long size = variableList->size();
        for(int j = 0; j < size; j++) {
            bff += (*nHat)[(*variableList)[j] - 1];
        }
        if((int)bff != 0) return false;
    }
    return true;
}
bool LDPCDecoder2::isZero(const vector<BinaryFiniteField>* nHat) const {    // 速度重視用。誤訂正を検出できない上、全ゼロ符号語にしか通用しないので注意
    unsigned long size = H->col();
    for(int i = 0; i < size; i++) {
        if((*nHat)[i] != 0) return false;
    }
    return true;
}
