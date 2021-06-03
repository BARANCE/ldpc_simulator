#include "DensityEvolution.h"

long double DensityEvolution::execDESimulation(DESInput simulationParameter) {
    startClock = clock();    // 時間計測開始

    //==== 変数宣言・前処理 ====//
    string fileName = simulationParameter.checkMatrixFileName;
    unsigned long hr = H->row();
    unsigned long hc = H->col();
    long double codeRate = (long double)(hc - hr) / hc;

    long double leftSuccess = 0.0;    // 左側成功確率
    long double rightFailure = 1.0;    // 右側失敗確率

    long double zeroBorderProb = simulationParameter.zeroBorderProb;    // 0とみなす確率
    unsigned long convergenceFailNumber = simulationParameter.convergenceFailNumber;    // 収束失敗回数
    unsigned long maxLoop = simulationParameter.binarySearchLoop;    // 二分探索の回数の上限
    DEInput dei = {0.0, zeroBorderProb, convergenceFailNumber};
    DEOutput deo;

    long double borderFripProbability;    // 閾値

    printProgress(0.0);
    cout.precision(10);
    for(unsigned long i = 0; i < maxLoop; i++) {
        //==== 経過時間 ====//
        nowClock = clock();
        unsigned long span = nowClock - startClock;
        unsigned int tmpExecTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
        unsigned int tmpExecTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
        unsigned int tmpExecTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;

        cout << "プロトグラフ：" << fileName << "(" << hr << " × " << hc << ", R = " << codeRate << ")\n";
        cout << "経過時間：" << tmpExecTimeHour << "時間" << tmpExecTimeMin << "分" << tmpExecTimeSec << "秒\n";
        cout << "試行回数：" << i + 1 << " (最大：" << maxLoop << ")\n";
        long double centerProb = (leftSuccess + rightFailure) / 2.0;
        cout << "左側：" << leftSuccess << "\n";
        cout << "右側：" << rightFailure << "\n";
        cout << "探索中：" << centerProb << "\n";

        dei.fripProbability = centerProb;
        deo = execDE(dei);
        if(deo.failFrag) {
            rightFailure = centerProb;
        } else {
            leftSuccess = centerProb;
        }

        printProgress((i + 1.0) / maxLoop);
    }
    //==== 実行時間 ====//
    endClock = clock();
    unsigned long span = endClock - startClock;
    unsigned int execTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
    unsigned int execTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
    unsigned int execTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;
    cout << "検査行列：" << fileName << "(" << hr << " × " << hc << ", R = " << codeRate << ")\n";
    cout << "経過時間：" << execTimeHour << "時間" << execTimeMin << "分" << execTimeSec << "秒\n";
    cout << "試行回数：" << maxLoop << " (最大：" << maxLoop << ")\n";
    cout << "左側：" << leftSuccess << "\n";
    cout << "右側：" << rightFailure << "\n";
    cout << "探索中：" << (leftSuccess + rightFailure) / 2 << "\n";


    borderFripProbability = leftSuccess;
    //cout << "閾値：" << borderFripProbability << "\n";

    return borderFripProbability;
}

DensityEvolution::DEOutput DensityEvolution::execDE(DEInput parameter) {
    //==== 変数宣言・前処理 ====//
    unsigned long hr = H->row();    // チェックノード数
    unsigned long hc = H->col();    // 変数ノード数
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;
    double p = parameter.fripProbability;    // 反転確率

    unsigned long loop = 0;    // 反復回数
    const unsigned long maxLoop = parameter.convergenceFailNumber;    // 最大反復回数(収束失敗条件)[debug]
    unsigned long conv = 0; // 収束回数

    long double border = parameter.borderProb;    // 0とみなす確率
    vector<BinaryFiniteField> borderChecker(hc);    // 変数ノードメッセージが境界値以下かどうかを格納
    vector<long double> probList(hc);    // 確変数ノードからの収束確率を保持するベクトル(返却用)

    DEOutput retObj = {vector<long double>(hc), false, 0};    // 返却構造体

    //==== 全ての変数ノードから隣接するチェックノードへpを送る ====//
    for(unsigned long i = 0; i < hc; i++) {
        const vector<int> checkList = (*matCol)[i];    // i番目の変数ノードと接続しているチェックノードの番号のリスト
        int size = checkList.size();
        for(int j = 0; j < size; j++) {
            int checkId = checkList[j] - 1;    // alist形式では、1から数え始めるので-1して添字をベクトルに合わせる
            mb.variableToCheck(i, checkId, p);    // i番目の変数ノードから、checkId番目のチェックノードに確率pを送信
        }
    }

    for(loop = 0; loop < maxLoop; loop++) {
        //==== チェックノード処理 ====//
        for(unsigned long i = 0; i < hr; i++) {
            const vector<int>* variableList = &((*matRow)[i]);
            unsigned long size = variableList->size();

            long double cProb = 1;
            for(unsigned long j = 0; j < size; j++) {
                cProb *= 1 - mb.checkBuf[i][j];
            }

            // 変数ノードへのメッセージを算出する
            for(unsigned long j = 0; j < size; j++) {
                unsigned long variableId = (*variableList)[j] - 1;
                long double writeValue = 1 - cProb;
                if(1 - mb.getCheckMessage(i, variableId) != 0.0) {
                    writeValue = 1 - (cProb / (1 - mb.getCheckMessage(i, variableId)));
                }
                mb.checkToVariable(i, variableId, writeValue);
            }
        }

        //==== 変数ノード処理 & 収束判定 ====//
        for(unsigned long i = 0; i < hc; i++) {
            //==== 変数ノード処理 ====//
            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            long double vProb = p;
            for(unsigned long j = 0; j < size; j++) {
                vProb *= mb.variableBuf[i][j];
            }

            // チェックノードへのメッセージを算出する
            for(unsigned long j = 0; j < size; j++) {
                unsigned long checkId = (*checkList)[j] - 1;
                long double writeValue = vProb;
                if(mb.getVariableMessage(i, checkId) != 0.0) {
                    writeValue = vProb / mb.getVariableMessage(i, checkId);
                }
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << "i = " << i << ", vProb = " << vProb << "\n";
            //==== 収束判定(境界値以下ならば1) ====//
            if(vProb < border) {
                borderChecker[i] = 1;
            } else {
                borderChecker[i] = 0;
            }
            probList[i] = vProb;
            //cout << (((int)borderChecker[i] == 1) ? "○" : "×") << "vProb[" << i << "] = " << vProb << "\n";
        }

        //==== 収束判定 ====//
        // borderCheckerの値が全て1ならば収束
        bool borderFrag = true;
        for(unsigned long i = 0; i < hc; i++) {
            if((int)borderChecker[i] == 0) {
                borderFrag = false;
                break;
            }
        }
        if(borderFrag) {
            break;
        }
    }

    //==== 収束回数の判定 ====//
    conv = loop + 1;

    //==== 実験データの整形 ====//
    for(unsigned long i = 0; i < hc; i++) {
        retObj.convProbList[i] = probList[i];
    }
    if(loop < maxLoop) {
        retObj.failFrag = false;
        //cout << "収束成功!!\n";
    } else {
        retObj.failFrag = true;
        //cout << "収束失敗・・・\n";
    }
    retObj.convergenceNumber = conv;
    return retObj;
}

DensityEvolution::DEWSOutput DensityEvolution::execDEwithStepCheck(DEInput parameter, string fileName, unsigned long defConvNum) {
    //==== 変数宣言・前処理 ====//
    unsigned long hr = H->row();    // チェックノード数
    unsigned long hc = H->col();    // 変数ノード数
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;
    double p = parameter.fripProbability;    // 反転確率

    unsigned long loop = 0;    // 反復回数
    const unsigned long maxLoop = parameter.convergenceFailNumber;    // 最大反復回数(収束失敗条件)[debug]
    unsigned long conv = 0; // 収束回数

    long double border = parameter.borderProb;    // 0とみなす確率
    vector<int> borderChecker(hc);    // 変数ノードメッセージが境界値以下かどうかを格納
    vector<long double> probList(hc);    // 確変数ノードからの収束確率を保持するベクトル(返却用)

    unsigned long maxWriteSize = maxLoop / 100 + 100;
    unsigned long writeCount = 0;    // データを取った個数
    DEWSOutput retObj = {vector<long double>(hc), false, 0, Matrix<long double>(maxWriteSize, hc), vector<long double>(maxWriteSize), vector<long double>(maxWriteSize)};    // 返却構造体

    vector<unsigned long> valueNumChecker(maxWriteSize);    // どの反復を記録したかを保持
    Matrix<long double> valueChecker(maxWriteSize, hc);
    vector<long double> maxValueChecker(maxWriteSize);    // 誤り率の最大値を保持
    vector<long double> minValueChecker(maxWriteSize);    // 誤り率の最小値を保持
    vector<long double> aveValueChecker(maxWriteSize);    // 平均誤り率を保持

    for(unsigned long i = 0; i < maxWriteSize; i++) maxValueChecker[i] = -1.0;
    for(unsigned long i = 0; i < maxWriteSize; i++) minValueChecker[i] = -1.0;
    for(unsigned long i = 0; i < maxWriteSize; i++) aveValueChecker[i] = 0.0;

    //==== 全ての変数ノードから隣接するチェックノードへpを送る ====//
    for(unsigned long i = 0; i < hc; i++) {
        const vector<int> checkList = (*matCol)[i];
        int size = checkList.size();
        for(int j = 0; j < size; j++) {
            int checkId = checkList[j] - 1;
            mb.variableToCheck(i, checkId, p);
        }
    }

    for(loop = 0; loop < maxLoop; loop++) {
        bool writeFrag = false;
        if(loop % 100 == 0) writeFrag = true;
        if(loop >= defConvNum - 1) writeFrag = true;
        if(loop < 10) writeFrag = true;
        if(loop > defConvNum - 50) writeFrag = true;

        //==== チェックノード処理 ====//
        for(unsigned long i = 0; i < hr; i++) {
            const vector<int>* variableList = &((*matRow)[i]);
            unsigned long size = variableList->size();

            long double cProb = 1;
            for(unsigned long j = 0; j < size; j++) {
                cProb *= 1 - mb.checkBuf[i][j];
            }

            // 変数ノードへのメッセージを算出する
            for(unsigned long j = 0; j < size; j++) {
                unsigned long variableId = (*variableList)[j] - 1;
                long double writeValue = 1 - cProb;
                if(1 - mb.getCheckMessage(i, variableId) != 0.0) {
                    writeValue = 1 - (cProb / (1 - mb.getCheckMessage(i, variableId)));
                }
                mb.checkToVariable(i, variableId, writeValue);
            }
        }

        //==== 変数ノード処理 & 収束判定 ====//
        for(unsigned long i = 0; i < hc; i++) {
            //==== 変数ノード処理 ====//
            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            long double vProb = p;
            for(unsigned long j = 0; j < size; j++) {
                vProb *= mb.variableBuf[i][j];
            }

            // チェックノードへのメッセージを算出する
            for(unsigned long j = 0; j < size; j++) {
                unsigned long checkId = (*checkList)[j] - 1;
                long double writeValue = vProb;
                if(mb.getVariableMessage(i, checkId) != 0.0) {
                    writeValue = vProb / mb.getVariableMessage(i, checkId);
                }
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << "i = " << i << ", vProb = " << vProb << "\n";
            //==== 収束判定(境界値以下ならば1) ====//
            if(vProb < border) {
                borderChecker[i] = 1;
            } else {
                borderChecker[i] = 0;
            }
            probList[i] = vProb;
            //cout << (((int)borderChecker[i] == 1) ? "○" : "×") << "vProb[" << i << "] = " << vProb << "\n";
            if(writeFrag) {
                valueChecker[writeCount][i] = vProb;
                if(maxValueChecker[writeCount] == -1.0) maxValueChecker[writeCount] = vProb;
                else if(maxValueChecker[writeCount] < vProb) maxValueChecker[writeCount] = vProb;
                if(minValueChecker[writeCount] == -1.0) minValueChecker[writeCount] = vProb;
                else if(minValueChecker[writeCount] > vProb) minValueChecker[writeCount] = vProb;
                aveValueChecker[writeCount] += vProb;
            }
        }
        if(writeFrag) {
            aveValueChecker[writeCount] /= hc;
            valueNumChecker[writeCount] = loop + 1;
            writeCount++;
        }

        //==== 収束判定 ====//
        // borderCheckerの値が全て1ならば収束
        bool borderFrag = true;
        for(unsigned long i = 0; i < hc; i++) {
            if((int)borderChecker[i] == 0) {
                borderFrag = false;
                break;
            }
        }
        if(borderFrag) {
            break;
        }
    }

    //==== 収束回数の判定 ====//
    conv = loop + 1;

    //==== 実験データの整形 ====//
    for(unsigned long i = 0; i < hc; i++) {
        retObj.convProbList[i] = probList[i];
    }
    if(loop < maxLoop) {
        retObj.failFrag = false;
        cout << "収束成功!!\n";
    } else {
        retObj.failFrag = true;
        cout << "収束失敗・・・\n";
    }
    retObj.convergenceNumber = conv;
    retObj.valueChecker = valueChecker;

    cout << "収束回数：" << conv << "回\n";
    cout << "ファイル書き出し開始\n";
    while(aveValueChecker[writeCount - 1] == 0.0) writeCount--;
    ofstream ofs("deDataFile1.txt");
    ofs << "#loop-> 0 ";
    for(unsigned long j = 0; j < writeCount; j++) {
        ofs << valueNumChecker[j] << " ";
    }
    ofs << "\n";
    for(unsigned long i = 0; i < valueChecker.col(); i++) {
        ofs << i << " " << p << " ";
        for(unsigned long j = 0; j < writeCount; j++) {
            ofs << valueChecker[j][i] << " ";
        }
        ofs << "\n";
    }
    ofs.close();
    ofstream ofs2("deDataFile2.txt");
    ofs2 << "0 " << p << " " << p << " " << p << "\n";
    for(unsigned long i = 0; i < writeCount; i++) {
        ofs2 << valueNumChecker[i] << " ";
        ofs2 << maxValueChecker[i] << " ";
        ofs2 << minValueChecker[i] << " ";
        ofs2 << aveValueChecker[i] << "\n";
    }
    ofs2.close();
    ofstream ofs3("dePlotFile.plt");
    ofs3 << "#!/I_download/gnuplot -persist\n";
    ofs3 << "set terminal postscript enhanced eps color solid\n";
    ofs3 << "set output \"deGraph.eps\"\n";
    ofs3 << "set size 0.7,0.7\n";
    ofs3 << "set key left bottom\n";
    ofs3 << "set grid xtics ytics mxtics mytics\n";
    ofs3 << "set xlabel \"Number of iteration\"\n";
    ofs3 << "set ylabel \"Probability sent erasure simbol\"\n";
    ofs3 << "set logscale y\n";
    ofs3 << "set format y \"10^{%L}\"\n";
    ofs3 << "set xrange [0:" << valueNumChecker[writeCount - 1] + valueNumChecker[writeCount - 1] / 50 << "]\n";
    ofs3 << "set yrange [" << aveValueChecker[writeCount - 1] / 10.0 << ":" << p * 10.0 << "]\n";
    ofs3 << "set style line 1 lt 1 lc 1 lw 3\n";
    ofs3 << "set style line 2 lt 1 lc 2 lw 3\n";
    ofs3 << "set style line 3 lt 1 lc 3 lw 3\n";
    ofs3 << "plot \"deDataFile2.txt\" using 1:2 w l ls 1 ti \"maxErasureRate\", \\\n";
    ofs3 << "\"deDataFile2.txt\" using 1:3 w l ls 2 ti \"minErasureRate\", \\\n";
    ofs3 << "\"deDataFile2.txt\" using 1:4 w l ls 3 ti \"aveErasureRate\"\n";
    ofs3 << "#\tEOF\n";
    ofs3.close();
    system("\"gnuplot dePlotFile.plt\"");
    cout << "ファイル書き出し完了\n";

    return retObj;
}

// ループ空間結合符号の形状を利用した複数の帯を接合した空間結合符号について、分離するまで回すDE
// connectionArea : 接合部の位置を表すベクトル
// blockEndPos : 各帯の終端部を表すベクトル
DensityEvolution::DEOutput DensityEvolution::execSubDE(DEInput parameter, vector<int> connectionArea, vector<int> blockEndPos) {
    //==== 変数宣言・前処理 ====//
    unsigned long hr = H->row();    // チェックノード数
    unsigned long hc = H->col();    // 変数ノード数
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;
    double p = parameter.fripProbability;    // 反転確率

    unsigned long loop = 0;    // 反復回数
    const unsigned long maxLoop = parameter.convergenceFailNumber;    // 最大反復回数(収束失敗条件)[debug]
    unsigned long conv = 0; // 収束回数

    long double border = parameter.borderProb;    // 0とみなす確率
    vector<BinaryFiniteField> borderChecker(hc);    // 変数ノードメッセージが境界値以下かどうかを格納
    vector<long double> probList(hc);    // 確変数ノードからの収束確率を保持するベクトル(返却用)

    unsigned long separateNum = 0;    // 分離した反復回数

    DEOutput retObj = {vector<long double>(hc), false, 0};    // 返却構造体

    //==== 全ての変数ノードから隣接するチェックノードへpを送る ====//
    for(unsigned long i = 0; i < hc; i++) {
        const vector<int> checkList = (*matCol)[i];    // i番目の変数ノードと接続しているチェックノードの番号のリスト
        int size = checkList.size();
        for(int j = 0; j < size; j++) {
            int checkId = checkList[j] - 1;    // alist形式では、1から数え始めるので-1して添字をベクトルに合わせる
            mb.variableToCheck(i, checkId, p);    // i番目の変数ノードから、checkId番目のチェックノードに確率pを送信
        }
    }

    for(loop = 0; loop < maxLoop; loop++) {
        //==== チェックノード処理 ====//
        for(unsigned long i = 0; i < hr; i++) {
            const vector<int>* variableList = &((*matRow)[i]);
            unsigned long size = variableList->size();

            long double cProb = 1;
            for(unsigned long j = 0; j < size; j++) {
                cProb *= 1 - mb.checkBuf[i][j];
            }

            // 変数ノードへのメッセージを算出する
            for(unsigned long j = 0; j < size; j++) {
                unsigned long variableId = (*variableList)[j] - 1;
                long double writeValue = 1 - cProb;
                if(1 - mb.getCheckMessage(i, variableId) != 0.0) {
                    writeValue = 1 - (cProb / (1 - mb.getCheckMessage(i, variableId)));
                }
                mb.checkToVariable(i, variableId, writeValue);
            }
        }

        //==== 変数ノード処理 & 収束判定 ====//
        for(unsigned long i = 0; i < hc; i++) {
            //==== 変数ノード処理 ====//
            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            long double vProb = p;
            for(unsigned long j = 0; j < size; j++) {
                vProb *= mb.variableBuf[i][j];
            }

            // チェックノードへのメッセージを算出する
            for(unsigned long j = 0; j < size; j++) {
                unsigned long checkId = (*checkList)[j] - 1;
                long double writeValue = vProb;
                if(mb.getVariableMessage(i, checkId) != 0.0) {
                    writeValue = vProb / mb.getVariableMessage(i, checkId);
                }
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << "i = " << i << ", vProb = " << vProb << "\n";
            //==== 収束判定(境界値以下ならば1) ====//
            if(vProb < border) {
                borderChecker[i] = 1;
            } else {
                borderChecker[i] = 0;
            }
            probList[i] = vProb;
            //cout << (((int)borderChecker[i] == 1) ? "○" : "×") << "vProb[" << i << "] = " << vProb << "\n";
        }

        //==== 収束判定 ====//
        // connectionAreaの要素に対応する変数ノードが全て収束しているかを判別
        bool separateFrag = true;
        int caSize = connectionArea.size();
        for(int i = 0; i < caSize; i++) {
            if((int)borderChecker[connectionArea[i]] == 0) {
                separateFrag = false;
                break;
            }
        }
        if(separateFrag) {
            separateNum = loop;    // 分離した反復回数
            cout << loop << "回目で分離!\n";
            break;
        }

        // borderCheckerの値が全て1ならば収束
        bool borderFrag = true;
        for(unsigned long i = 0; i < hc; i++) {
            if((int)borderChecker[i] == 0) {
                borderFrag = false;
                break;
            }
        }
        if(borderFrag) {
            break;
        }
    }

    //==== 収束回数の判定 ====//
    conv = loop + 1;

    //==== 実験データの整形 ====//
    for(unsigned long i = 0; i < hc; i++) {
        retObj.convProbList[i] = probList[i];
    }
    if(loop < maxLoop) {
        retObj.failFrag = false;
        //cout << "収束成功!!\n";
    } else {
        retObj.failFrag = true;
        //cout << "収束失敗・・・\n";
    }
    retObj.convergenceNumber = conv;
    return retObj;
}

// 深さ(depth) : 反転確率の境界値の精度(1.0 * 10^(-1 * depth)の精度で求める)
long double DensityEvolution::borderProb(int valDim, int checkDim, int depth) const {
    long double border = 0.0;    // 境界値の候補
    for(int i = 1; i <= depth; i++) {
        int j;    // j * pow(0.1, i)の値をborderに足す。
        for(j = 1; j <= 10; j++) {
            long double p = border + j * pow((long double)0.1, i);
            bool passChecker = deStep(valDim, checkDim, p);
            if(passChecker == false) {
                j = j - 1;
                break;
            }
        }
        border += j * pow((long double)0.1, i);
    }
    return border;
}

bool DensityEvolution::deStep(int valDim, int checkDim, long double fripProbability) const {
    bool convergence = false;    // 0に収束したらtrue
    long double p = fripProbability;
    unsigned long maxLoop = 100000;
    unsigned long loop = 1;
    long double border = pow((long double)0.1, 15);    // 0だと判定する境界値
    for(loop = 1; loop <= maxLoop; loop++) {
        p = deeq(p, fripProbability, checkDim, valDim);
        if(p < border) {
            //cout << "0に収束しました(fp = " << fripProbability <<", p = " << p << ", loop = " << loop << ")\n";
            convergence = true;
            break;
        }
    }
    if(loop >= maxLoop) {
        //cout << "0に収束しません(fp = " << fripProbability <<", p = " << p << ", loop = " << loop << ")\n";
    }

    if(convergence) return true;
    else return false;
}

//==== 進捗状況の表示 ====//
void DensityEvolution::printProgress(double progressRate) const {
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

//==== 実験結果を保存 ====//
void DensityEvolution::saveDEData(string fileName, double epsilon, unsigned long loop, unsigned long girth) const {
    ofstream ofs("deResult.txt", fstream::app);
    if(ofs.fail()) cout<< "ファイルを開けません。\n";
    else {
        ofs.precision(10);
        ofs << fileName << " " << epsilon << " " << loop << " " << girth << "\n";
    }
    ofs.close();
}
