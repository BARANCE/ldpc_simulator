#include "DecodeSimulator2.h"
#include <math.h>

//==== コンストラクタ ====//
DecodeSimulator2::DecodeSimulator2(const SPMatrix* checkMatrix) : H(checkMatrix) {
    randFrag = 0;
}

//==== 復号実験 ====//
DecodeSimulator2::DSOutput DecodeSimulator2::decodingSimulation(DSInput simulationParameter) {
    startClock = clock();    // 時間計測開始

    //==== 変数宣言・前処理 ====//
    unsigned long hr = H->row();
    unsigned long hc = H->col();
    vector<BinaryFiniteField> codeWord = simulationParameter.cordWord;    // 符号語
    vector<double> transWord(hc);    // 送信語

    //==== 通信路の設定 ====//
    EChannel channel = simulationParameter.channel;
    double SNR = simulationParameter.SNR;    // SN比[dB] or 反転確率・消失確率
    const double codeRate = (double)(hc - hr) / (double)hc;    // 符号化率
    double variance = 0.5 * (1.0 / pow(10.0, SNR / 10.0)) * (double)hc / (double)(hc - hr);    // SN比から求めた雑音の分散
    double averageOfNoise = 0.0;    // 雑音の標本平均
    double varianceOfNoise = 0.0;    // 雑音の標本分散

    unsigned long maxLoop = simulationParameter.maxLoop;    // 最大試行回数
    unsigned long loop = 0;    // 試行回数
    double scale = simulationParameter.scale;

    unsigned long convergenceFailNumber = simulationParameter.convergenceFailNumber;

    //==== 検証用変数 ====//
    unsigned long maxError = simulationParameter.maxError;    // 最大ブロック誤り回数
    unsigned long numberOfSuccess = 0;    // 復号が成功した回数
    unsigned long numberOfBlockError = 0;    // ブロック誤りが発生した回数
    unsigned long numberOfBitError = 0;    // ビット誤りが発生した回数
    unsigned long error1 = 0;    // 誤訂正誤りが発生した回数
    unsigned long error2 = 0;    // 収束失敗による誤りが発生した回数
    unsigned long error3 = 0;    // 誤収束による誤りが発生した回数
    unsigned long maxConvergenceNumber = 0;    // 最大収束回数
    unsigned long minConvergenceNumber = convergenceFailNumber;    // 最小収束回数
    unsigned long averageConvergenceNumber = 0;    // 平均収束回数
    //unsigned long averageOfBitError;    // ビット誤りの平均個数
    unsigned long minOfBitError = hc;    // ビット誤りの最小個数
    unsigned long maxOfBitError = 0;    // ビット誤りの最大個数
    vector<unsigned long> distributionOfNumberOfBitError(hc, 0);    // ビット誤り個数の分布
    unsigned long averageOfBurstLength = 0;    // バースト長の平均
    unsigned long errorFirst = hc;
    unsigned long errorEnd = 0;
    DSOutput retObj;    // 返却構造体
    int convCheckFrag = simulationParameter.convCheckFrag;    // 収束チェックをするかどうか
    int convCheckSize = simulationParameter.convCheckSize;
    vector<vector<unsigned long> > convCheckCounter(convCheckSize, vector<unsigned long>(hc, 0));    // 反復回数ごとの誤った回数を記録
    vector<vector<double> > convCheckProb(convCheckSize, vector<double>(hc, 0.0));    // 反復回数ごとの誤り率を記録
    LDPCDecoder2::ECDList researchParameter = {convCheckFrag, convCheckSize};

    //==== 特殊な復号実験用のパラメータ ====//
    double shortenRate = simulationParameter.shortenRate;

    //==== 復号器の設定 ====//
    EDMethod decodingMethod = simulationParameter.decodingMethod;    // 復号法
    LDPCDecoder2 dec(H);    // 復号器
    //dec.setChannel(LDPCDecoder2::EChannel::C_AWGN);    // AWGN通信路
    dec.setSNR(SNR);    // SN比をセット
    LDPCDecoder2::ESPInput inputParameter = {vector<double>(hc, 0.0), convergenceFailNumber};    // 復号器への入力データ(送信語, 収束失敗条件)
    // ※送信語は雑音を付加して生成するため、ここでは初期化のみ行う。
    LDPCDecoder2::ESPOutput outputParameter;    // 復号器からの出力データ
    LDPCDecoder2::ERSPOutput outputBECParameter;

    //==== ショートン復号専用設定 ====//
    vector<int> shortenIds(hc);    // ビットをショートンする箇所
    if(decodingMethod == C_MMIN_SUM) {
        //const unsigned long shortenRate = 400;
        for(unsigned long i = 0; i < hc; i++) {
            /*if(i < (hc / shortenRate) || i > hc - (hc / shortenRate) ) {
                shortenIds[i] = 1;
                codeWord[i] = 0;
            } else {
                shortenIds[i] = 0;
            }*/
            unsigned long blockNumb = (unsigned long)(hc * shortenRate);
            unsigned long blockId = (unsigned long)((i * 1.0) / (double)blockNumb);
            if(i == (int)(blockNumb * blockId)) {
                shortenIds[i] = 1;
                codeWord[i] = 0;
            } else {
                shortenIds[i] = 0;
            }
        }

    }

    //==== バイナリ-バイポーラ変換(BPSKモデル) ====//
    for(unsigned long i = 0; i < hc; i++) {
        if((int)codeWord[i] == 0) transWord[i] = 1.0;
        else transWord[i] = -1.0;
    }

    for(loop = 0; loop < maxLoop; loop++) {

        //==== 雑音付加(雑音の統計処理を含む) ====//
        vector<double> noiseList(hc);
        double noiseAve = 0.0;
        double noiseVar = 0.0;
        if(channel == C_AWGN) {
            for(unsigned long i = 0; i < hc; i++) {
                inputParameter.receptionWord[i] = transWord[i];
                noiseList[i] = getGaussianNoise(0, variance);
                noiseAve += noiseList[i];
                inputParameter.receptionWord[i] += noiseList[i];
            }
            noiseAve /= hc;
            for(unsigned long i = 0; i < hc; i++) {
                noiseVar += pow(noiseAve - noiseList[i], 2.0);
            }
            noiseVar /= hc;
            averageOfNoise += noiseAve;
            varianceOfNoise += noiseVar;
        } else if(channel == C_BEC) {
            for(unsigned long i = 0; i < hc; i++) {
                if((int)codeWord[i] == 0) inputParameter.receptionWord[i] = ZERO;
                else if((int)codeWord[i] == 1) inputParameter.receptionWord[i] = ONE;
                double rand = genrand_real1(); //一様実乱数[0,1] (32ビット精度)
                if(SNR > rand) inputParameter.receptionWord[i] = ERASURE;
            }
        }

        //if(channel == EChannel::C_AWGN) cout << "AWGN\n";
        //else if(channel == EChannel::C_BEC) cout << "BEC\n";

        //==== 復号 ====//
        if(channel == C_AWGN) {
            if(decodingMethod == C_SUM_PRODUCT) {
                outputParameter = dec.sumProductDecode(inputParameter);
            } else if(decodingMethod == C_MIN_SUM) {
                outputParameter = dec.minSumDecode(inputParameter, scale);
            } else if(decodingMethod == C_HARD_DECISION) {
                outputParameter = dec.hardDecision(inputParameter);
            } else if(decodingMethod == C_MMIN_SUM) {
                outputParameter = dec.iterativeMinSumDecode(inputParameter, scale, NULL);
            } else if(decodingMethod == C_SHORTEN_MIN_SUM) {
                outputParameter = dec.shortenMinSumDecode(inputParameter, scale, shortenRate);
            }
        } else if(channel == C_BEC) {
            outputBECParameter = dec.sumProductDecodeOnBEC(inputParameter, researchParameter);
            outputParameter.convergenceNumber = outputBECParameter.convergenceNumber;
            outputParameter.errorType = outputBECParameter.errorType;
            outputParameter.failFrag = outputBECParameter.failFrag;
        }

        //==== 検証 ====//
        vector<BinaryFiniteField>* decodeWord = &outputParameter.decodeWord;
        unsigned long convn = outputParameter.convergenceNumber;
        if(channel == C_AWGN) {
        if(checkDecodeWord(&codeWord, decodeWord)) {    // 復号成功
            numberOfSuccess++;
            if(maxConvergenceNumber < convn) maxConvergenceNumber = convn;    // 最大収束回数を更新
            if(minConvergenceNumber > convn) minConvergenceNumber = convn;    // 最小収束回数を更新
            averageConvergenceNumber += convn;    // 平均収束回数を更新

            if(loop % 10000 == 0) {
                //==== 経過時間 ====//
                nowClock = clock();
                unsigned long span = nowClock - startClock;
                unsigned int tmpExecTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
                unsigned int tmpExecTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
                unsigned int tmpExecTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;

                //==== 表示 ====//
                double progressRate = (double)numberOfBlockError / maxError;
                printProgress(progressRate);
                cout << "検査行列：" << simulationParameter.fileName << "(" << hr << " × " << hc << ", R = " << codeRate << ")\n";
                cout << "経過時間：" << tmpExecTimeHour << "時間" << tmpExecTimeMin << "分" << tmpExecTimeSec << "秒\n";
                cout << "SN比：" << SNR << " (分散：" << variance << ")\n";
                cout << "雑音標本平均：" << averageOfNoise / (loop + 1.0) << ", 標本分散：" << varianceOfNoise / (loop + 1.0) << "\n";
                cout << "試行回数：" << loop + 1 << " (最大：" << maxLoop << ")\n";
                cout << "ブロック誤り率：" << numberOfBlockError / (loop + 1.0) << "\n";
                cout << "ビット誤り率：" << numberOfBitError / ((loop + 1.0) * (double)hc) << "\n";
                cout << "エラー別誤り率：" << (double)error1 / (loop + 1.0) << " / " << (double)error2 / (loop + 1.0) << " / " << (double)error3 / (loop + 1.0) << "\n";
                cout << "平均誤りビット数：" << (double)numberOfBitError / numberOfBlockError << "(最小：" << minOfBitError << " / 最大：" << maxOfBitError <<")\n";
                cout << "平均誤りバースト長：" << (double)averageOfBurstLength / numberOfBlockError << "(先頭：" << errorFirst << " / 最終：" << errorEnd << ")\n";
                cout << "平均収束回数：" << (double)averageConvergenceNumber / numberOfSuccess << "(最小：" << minConvergenceNumber << " / 最大：" << maxConvergenceNumber << ")\n";
            }
        } else {    // 復号失敗
            //==== 復号成功率 ====//
            double tmpSuccess = (double)numberOfSuccess / (loop + 1.0);    // 暫定復号成功率

            //==== ブロック誤り率の測定 ====//
            numberOfBlockError++;
            double tmpBlockError = (double)numberOfBlockError / (loop + 1.0);    // 暫定ブロック誤り率
            double progressRate = (double)numberOfBlockError / maxError;    // 進行状況

            //==== ビット誤り数の測定 ====//
            unsigned long bitCount = 0;    // 誤りの個数
            unsigned long firstBitError = hc;    // 最初の誤りの位置
            unsigned long finalBitError = 0;    // 最後の誤りの位置
            unsigned long burstLength = 0;    // 誤りの長さ
            bool firstFrag = true;    // 最初の誤りであるかどうか
            for(unsigned long i = 0; i < hc; i++) {
                if(codeWord[i] != outputParameter.decodeWord[i]) {
                    finalBitError = i + 1;
                    bitCount++;
                    if(firstFrag) {
                        firstFrag = false;
                        firstBitError = i;
                    }
                    if(i < errorFirst) errorFirst = i;
                    if(i > errorEnd) errorEnd = i;
                }
            }
            burstLength = finalBitError - firstBitError;
            averageOfBurstLength += burstLength;

            numberOfBitError += bitCount;
            if(bitCount < minOfBitError) minOfBitError = bitCount;
            if(bitCount > maxOfBitError) maxOfBitError = bitCount;

            distributionOfNumberOfBitError[bitCount - 1]++;    // ビット誤り個数の分布(0要素目に誤り1個の個数が入る)

            double tmpBitError = (double)numberOfBitError / ((loop + 1.0) * (double)hc);    // 暫定ビット誤り率
            double tmpAveNumBitError = (double)numberOfBitError / numberOfBlockError;    // 誤りビットの平均個数
            double tmpBurstLength = (double)averageOfBurstLength / numberOfBlockError;    // バースト長の平均

            //==== 誤りの種類の測定 ====//
            if(outputParameter.failFrag) {
                if(outputParameter.errorType == LDPCDecoder2::OVER_MAX_CONVERGENCE) {    // 収束失敗
                    error2++;
                } else if(outputParameter.errorType == LDPCDecoder2::BAD_CONVERGENCE) {    // 誤収束
                    error3++;
                } else {
                    cout << "Warning : 規定されていない誤りが発生しています。\n";
                }
            } else {    // 誤訂正
                error1++;
            }

            double tmpError1 = (double)error1 / (loop + 1.0);
            double tmpError2 = (double)error2 / (loop + 1.0);
            double tmpError3 = (double)error3 / (loop + 1.0);

            //==== 収束回数 ====//
            double tmpAveConvergence = (double)averageConvergenceNumber / numberOfSuccess;

            //==== 実行時間 ====//
            nowClock = clock();
            unsigned long span = nowClock - startClock;
            unsigned int tmpExecTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
            unsigned int tmpExecTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
            unsigned int tmpExecTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;

            //==== 表示 ====//
            printProgress(progressRate);
            cout << "検査行列：" << simulationParameter.fileName << "(" << hr << " × " << hc << ", R = " << codeRate << ")\n";
            cout << "経過時間：" << tmpExecTimeHour << "時間" << tmpExecTimeMin << "分" << tmpExecTimeSec << "秒\n";
            cout << "SN比：" << SNR << " (分散：" << variance << ")\n";
            cout << "雑音標本平均：" << averageOfNoise / (loop + 1.0) << ", 標本分散：" << varianceOfNoise / (loop + 1.0) << "\n";
            cout << "試行回数：" << loop + 1 << " (最大：" << maxLoop << ")\n";
            cout << "ブロック誤り率：" << tmpBlockError << "\n";
            cout << "ビット誤り率：" << tmpBitError << "\n";
            cout << "エラー別誤り率：" << tmpError1 << " / " << tmpError2 << " / " << tmpError3 << "\n";
            cout << "平均誤りビット数：" << tmpAveNumBitError << "(最小：" << minOfBitError << " / 最大：" << maxOfBitError <<")\n";
            cout << "平均誤りバースト長：" << tmpBurstLength << "(先頭：" << errorFirst << " / 最終：" << errorEnd << ")\n";
            cout << "平均収束回数：" << tmpAveConvergence << "(最小：" << minConvergenceNumber << " / 最大：" << maxConvergenceNumber << ")\n";

            //==== 終了条件 ====//
            if(numberOfBlockError >= maxError) break;
        }
        } else if(channel == C_BEC) {
            vector<double> decodeWord = outputBECParameter.decodeWord;    // 復号語
            int errorCounter = 0;    // 消失シンボルの数
            for(unsigned long i = 0; i < hc; i++) {
                if(decodeWord[i] == ERASURE) errorCounter++;
            }
            vector<vector<double> > convCheckVec = outputBECParameter.convCheckVec;    // 収束状況100回まで記録
            if(convCheckFrag) {
                for(int i = 0; i < convCheckSize; i++) {
                    for(unsigned long j = 0; j < hc; j++) {
                        if(convCheckVec[i][j] == ERASURE) {
                            convCheckCounter[i][j]++;
                        }
                    }
                }
            }


            if(errorCounter == 0) {    // 復号成功
                numberOfSuccess++;
                if(maxConvergenceNumber < convn) maxConvergenceNumber = convn;    // 最大収束回数を更新
                if(minConvergenceNumber > convn) minConvergenceNumber = convn;    // 最小収束回数を更新
                averageConvergenceNumber += convn;    // 平均収束回数を更新

                if(loop % 10000 == 0) {
                    //==== 経過時間 ====//
                    nowClock = clock();
                    unsigned long span = nowClock - startClock;
                    unsigned int tmpExecTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
                    unsigned int tmpExecTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
                    unsigned int tmpExecTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;

                    //==== 表示 ====//
                    double progressRate = (double)numberOfBlockError / maxError;
                    printProgress(progressRate);
                    cout << "検査行列：" << simulationParameter.fileName << "(" << hr << " × " << hc << ", R = " << codeRate << ")\n";
                    cout << "経過時間：" << tmpExecTimeHour << "時間" << tmpExecTimeMin << "分" << tmpExecTimeSec << "秒\n";
                    cout << "消失確率：" << SNR << " (分散：" << variance << ")\n";
                    cout << "雑音標本平均：" << averageOfNoise / (loop + 1.0) << ", 標本分散：" << varianceOfNoise / (loop + 1.0) << "\n";
                    cout << "試行回数：" << loop + 1 << " (最大：" << maxLoop << ")\n";
                    cout << "ブロック誤り率：" << numberOfBlockError / (loop + 1.0) << "\n";
                    cout << "ビット誤り率：" << numberOfBitError / ((loop + 1.0) * (double)hc) << "\n";
                    cout << "エラー別誤り率：" << (double)error1 / (loop + 1.0) << " / " << (double)error2 / (loop + 1.0) << " / " << (double)error3 / (loop + 1.0) << "\n";
                    cout << "平均誤りビット数：" << (double)numberOfBitError / numberOfBlockError << "(最小：" << minOfBitError << " / 最大：" << maxOfBitError <<")\n";
                    cout << "平均誤りバースト長：" << (double)averageOfBurstLength / numberOfBlockError << "(先頭：" << errorFirst << " / 最終：" << errorEnd << ")\n";
                    cout << "平均収束回数：" << (double)averageConvergenceNumber / numberOfSuccess << "(最小：" << minConvergenceNumber << " / 最大：" << maxConvergenceNumber << ")\n";
                }
            } else {    // 復号失敗
                //==== 復号成功率 ====//
                double tmpSuccess = (double)numberOfSuccess / (loop + 1.0);    // 暫定復号成功率

                //==== ブロック誤り率の測定 ====//
                numberOfBlockError++;
                double tmpBlockError = (double)numberOfBlockError / (loop + 1.0);    // 暫定ブロック誤り率
                double progressRate = (double)numberOfBlockError / maxError;    // 進行状況

                //==== ビット誤り率 ====//
                unsigned long bitCount = 0;    // 誤りの個数
                unsigned long firstBitError = hc;    // 最初の誤りの位置
                unsigned long finalBitError = 0;    // 最後の誤りの位置
                unsigned long burstLength = 0;    // 誤りの長さ
                bool firstFrag = true;    // 最初の誤りであるかどうか
                for(unsigned long i = 0; i < hc; i++) {
                    if(decodeWord[i] == ERASURE) {
                        finalBitError = i + 1;
                        bitCount++;
                        if(firstFrag) {
                            firstFrag = false;
                            firstBitError = i;
                        }
                        if(i < errorFirst) errorFirst = i;
                        if(i > errorEnd) errorEnd = i;
                    }
                }
                burstLength = finalBitError - firstBitError;
                averageOfBurstLength += burstLength;

                numberOfBitError += bitCount;
                if(bitCount < minOfBitError) minOfBitError = bitCount;
                if(bitCount > maxOfBitError) maxOfBitError = bitCount;

                distributionOfNumberOfBitError[bitCount - 1]++;    // ビット誤り個数の分布(0要素目に誤り1個の個数が入る)

                double tmpBitError = (double)numberOfBitError / ((loop + 1.0) * (double)hc);    // 暫定ビット誤り率
                double tmpAveNumBitError = (double)numberOfBitError / numberOfBlockError;    // 誤りビットの平均個数
                double tmpBurstLength = (double)averageOfBurstLength / numberOfBlockError;    // バースト長の平均

                //==== 誤りの種類の測定 ====//
                if(outputParameter.failFrag) {
                    if(outputParameter.errorType == LDPCDecoder2::OVER_MAX_CONVERGENCE) {    // 収束失敗
                        error2++;
                    } else if(outputParameter.errorType == LDPCDecoder2::BAD_CONVERGENCE) {    // 誤収束
                        error3++;
                    } else {
                        cout << "Warning : 規定されていない誤りが発生しています。\n";
                    }
                } else {    // 誤訂正
                    error1++;
                }

                double tmpError1 = (double)error1 / (loop + 1.0);
                double tmpError2 = (double)error2 / (loop + 1.0);
                double tmpError3 = (double)error3 / (loop + 1.0);

                //==== 収束回数 ====//
                double tmpAveConvergence = (double)averageConvergenceNumber / numberOfSuccess;

                //==== 実行時間 ====//
                nowClock = clock();
                unsigned long span = nowClock - startClock;
                unsigned int tmpExecTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
                unsigned int tmpExecTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
                unsigned int tmpExecTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;

                //==== 表示 ====//
                printProgress(progressRate);
                cout << "検査行列：" << simulationParameter.fileName << "(" << hr << " × " << hc << ", R = " << codeRate << ")\n";
                cout << "経過時間：" << tmpExecTimeHour << "時間" << tmpExecTimeMin << "分" << tmpExecTimeSec << "秒\n";
                cout << "消失確率：" << SNR << " (分散：" << variance << ")\n";
                cout << "雑音標本平均：" << averageOfNoise / (loop + 1.0) << ", 標本分散：" << varianceOfNoise / (loop + 1.0) << "\n";
                cout << "試行回数：" << loop + 1 << " (最大：" << maxLoop << ")\n";
                cout << "ブロック誤り率：" << tmpBlockError << "\n";
                cout << "ビット誤り率：" << tmpBitError << "\n";
                cout << "エラー別誤り率：" << tmpError1 << " / " << tmpError2 << " / " << tmpError3 << "\n";
                cout << "平均誤りビット数：" << tmpAveNumBitError << "(最小：" << minOfBitError << " / 最大：" << maxOfBitError <<")\n";
                cout << "平均誤りバースト長：" << tmpBurstLength << "(先頭：" << errorFirst << " / 最終：" << errorEnd << ")\n";
                cout << "平均収束回数：" << tmpAveConvergence << "(最小：" << minConvergenceNumber << " / 最大：" << maxConvergenceNumber << ")\n";

                //==== 終了条件 ====//
                if(numberOfBlockError >= maxError) break;
            }
        }

        //break;
    }

    // 収束状況のチェック
    if(convCheckFrag) {
        for(int i = 0; i < convCheckSize; i++) {
            for(unsigned long j = 0; j < hc; j++) {
                convCheckProb[i][j] = convCheckCounter[i][j] / (double)(loop + 1.0);
            }
        }
    }

    //==== 実験結果の整形 ====//
    retObj.loop = loop + 1;    // 実際の試行回数
    retObj.blockErrorRate = (double)numberOfBlockError / (loop + 1.0);    // ブロック誤り率
    retObj.bitErrorRate = (double)numberOfBitError / ((loop + 1.0) * (double)hc);    // ビット誤り率
    retObj.decodeSuccessRate = (double)numberOfSuccess / (loop + 1.0);    // 復号成功率
    retObj.error1Rate = (double)error1 / (loop + 1.0);    // 誤訂正率
    retObj.error2Rate = (double)error2 / (loop + 1.0);    // 収束失敗率
    retObj.error3Rate = (double)error3 / (loop + 1.0);    // 誤収束率
    retObj.averageConvergenceNumber = (double)averageConvergenceNumber / numberOfSuccess;    // 平均収束回数
    retObj.minConvergenceNumber = minConvergenceNumber;    // 最小収束回数(復号が一度も成功していない場合は符号長に等しい)
    retObj.maxConvergenceNumber = maxConvergenceNumber;    // 最大収束回数(復号が一度も成功していない場合は0)
    retObj.distributionOfNumberOfBitError = distributionOfNumberOfBitError;    // 誤りビット数の分布
    retObj.averageOfNumberOfBitError = (double)numberOfBitError / numberOfBlockError;    // 平均誤りビット数
    retObj.averageOfGaussianNoise = averageOfNoise / (loop + 1.0);    // 雑音の平均
    retObj.varianceOfGaussianNoise = varianceOfNoise / (loop + 1.0);    // 雑音の分散
    retObj.averageOfBurstLength = (double)averageOfBurstLength / numberOfBlockError;    // 平均誤りバースト長
    retObj.execTime = startClock - clock();    // 実行時間

    //==== 結果の出力 ====//
    stringstream resultss;
    resultss << "result_" << simulationParameter.fileName;
    string resultFileName = resultss.str();    // 結果のファイル名
    ofstream ofs(resultFileName.c_str(), fstream::app);
        ofs << SNR << " ";
        ofs << retObj.blockErrorRate << " ";
        ofs << retObj.bitErrorRate << " ";
        ofs << retObj.loop << " ";
        ofs << retObj.decodeSuccessRate << " ";
        ofs << retObj.error1Rate << " ";
        ofs << retObj.error2Rate << " ";
        ofs << retObj.error3Rate << " ";
        ofs << retObj.averageConvergenceNumber << " ";
        ofs << retObj.minConvergenceNumber << " ";
        ofs << retObj.maxConvergenceNumber << " ";
        ofs << retObj.averageOfNumberOfBitError << " ";
        ofs << retObj.averageOfGaussianNoise << " ";
        ofs << retObj.varianceOfGaussianNoise << " ";
        ofs << retObj.averageOfBurstLength << " ";
        ofs << retObj.execTime << "\n";
    ofs.close();

    stringstream resultss2;    // 誤り分布用
    resultss2 << "dis_" << simulationParameter.fileName;
    string disFileName = resultss2.str();
    ofstream ofs2(disFileName.c_str(), fstream::app);
        ofs2 << SNR << " ";
        for(unsigned long i = 0; i < hc; i++) {
            ofs2 << retObj.distributionOfNumberOfBitError[i];
            if (i < hc - 1) ofs2 << " ";
        }
        ofs2 << "\n";
    ofs2.close();

    if(convCheckFrag) {
        stringstream resultss3;    // 収束チェック用
        resultss3 << "conv_SNR" << SNR << "_" << removeExtension(simulationParameter.fileName) << ".txt";
        string convFileName = resultss3.str();
        ofstream ofs3(convFileName.c_str(), fstream::out);    // 書き込み用
        ofs3.precision(10);
        //ofs3.setf(ios::fixed);
        ofs3 << "#node ";
        for(int j = 0; j < convCheckSize; j++) {
            ofs3 << j;
            if(j < convCheckSize - 1) ofs3 << " ";
        }
        ofs3 << "\n";
        for(unsigned long i = 0; i < hc; i++) {
            ofs3 << i << " ";
            for(int j = 0; j < convCheckSize; j++) {
                ofs3 << convCheckProb[j][i];
                if(j < convCheckSize - 1) ofs3 << " ";
            }
            ofs3 << "\n";
        }
        ofs3.close();

        stringstream resultss4;    // 収束チェック用2
        resultss4 << "conv_EP" << SNR << "_" << removeExtension(simulationParameter.fileName) << ".plt";
        string plotName = resultss4.str();
        ofstream ofs4(plotName.c_str(), fstream::out);    // gnuplot用
        ofs4 << "#!/I_download/gnuplot -persist\n";
        ofs4 << "set terminal postscript enhanced eps color solid\n";
        ofs4 << "set output \"conv_SNR" << SNR << "_" << removeExtension(simulationParameter.fileName) << ".eps\"\n";
        ofs4 << "set xlabel \"variable node\"\n";
        ofs4 << "set ylabel \"erasure message probability\"\n";
        ofs4 << "set xrange [0:" << hc - 1 << "]\n";
        ofs4 << "set yrange [1e-7:1]\n";
        ofs4 << "set logscale y\n";
        ofs4 << "set format y \"10^{%L}\"\n";
        ofs4 << "set size 1.0,0.7\n";
        ofs4 << "set nokey\n";
        ofs4 << "set grid xtics ytics mxtics mytics\n";
        ofs4 << "set style line 1 lt 1 lc 1 lw 3\n";
        ofs4 << "set style line 2 lt 1 lc rgb \"#32cd32\" lw 3\n";
        ofs4 << "set style line 3 lt 1 lc 3 lw 3\n";
        ofs4 << "set style line 4 lt 1 lc 4 lw 3\n";
        ofs4 << "set style line 5 lt 1 lc rgb \"#20b2aa\" lw 3\n";
        ofs4 << "set style line 6 lt 1 lc rgb \"#b8860b\" lw 3\n";
        ofs4 << "set style line 7 lt 1 lc 7 lw 3\n";
        ofs4 << "set style line 8 lt 1 lc 8 lw 3\n";
        ofs4 << "set style line 9 lt 1 lc 9 lw 3\n";
        ofs4 << "plot ";
        for(int i = 0; i < convCheckSize; i++) {
            ofs4 << "\"" << convFileName << "\" using 1:" << i + 2 << " w steps ls " << (i % 9) + 1 << " ti \"iteration" << i << "\"";
            if(i < convCheckSize - 1) ofs4 << ", \\";
            ofs4 << "\n";
        }
        ofs4 << "#\tEOF";
        ofs4.close();

        stringstream resultss6;
        resultss6 << "\"gnuplot " << plotName << "\"";
        string systemName = resultss6.str();
        system(systemName.c_str());

        stringstream resultss7;
        resultss7 << "\"DEL " << plotName << "\"";
        string delName = resultss7.str();
        system(delName.c_str());
    }

    return retObj;
}

//==== 雑音生成 ====//
double DecodeSimulator2::getGaussianNoise(double average, double variance) const {
    double retVal;
    double r1 = 1.0 - genrand_real2();    // (0, 1]一様乱数
    double r2 = 1.0 - genrand_real2();
    if(randFrag == 0) {
        retVal = sin(2.0 * 3.14159265 * r2);
    } else {
        retVal = cos(2.0 * 3.14159265 * r2);
    }
    retVal *= sqrt(-2.0 * log(r1)) * sqrt(variance);
    retVal += average;

    return retVal;
}
DecodeSimulator2::DSAveVar DecodeSimulator2::checkAverageAndVariance(const vector<double>* vec) const {
    unsigned long size = vec->size();
    double average = 0.0;
    double variance = 0.0;
    for(unsigned long i = 0; i < size; i++) {
        average += (*vec)[i];
    }
    average /= size;
    for(unsigned long i = 0; i < size; i++) {
        variance += pow(average - (*vec)[i], 2);
    }
    variance /= size;
    DSAveVar retObj = {average, variance};
    return retObj;
}

//==== 復号語と送信符号語が等しいかどうか ====//
bool DecodeSimulator2::checkDecodeWord(const vector<BinaryFiniteField>* codeWord, const vector<BinaryFiniteField>* recepWord) const {
    unsigned long size = codeWord->size();
    for(unsigned long i = 0; i < size; i++) {
        if((*codeWord)[i] != (*recepWord)[i]) return false;
    }
    return true;
}

//==== 進捗状況の表示 ====//
void DecodeSimulator2::printProgress(double progressRate) const {
    int maxLength = 30;
    int length = (int)(maxLength * progressRate);
    system("cls");
    for(int i = 0; i < length; i++) {
        cout << "■";
    }
    for(int i = length; i < maxLength; i++) {
        cout << "□";
    }
    cout << " (" << progressRate * 100 << "%完了)\n";
}

//==== strから拡張子を取り除く ====//
string DecodeSimulator2::removeExtension(string str) const {
    int pos = str.find_last_of(".");
    if(pos != -1) {
        return str.substr(0, pos);
    } else {
        return str;
    }
}
