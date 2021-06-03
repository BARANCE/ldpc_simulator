#include "LoopLDPC.h"

//==== コンストラクタ ====//
LoopLDPC::LoopLDPC(const vector<LBandData>& input) : bandData(input), mat(calculateRowSize(input), calculateColSize(input)), bsrp(input.size()), bscp(input.size()) {
    calculateBandStartPos();
    generateBaseSCCProtograph();
    generateConnections();
    //printBandData();
    //cout << "m = " << mat.row() << ", n = " << mat.col() << "\n";
}

//==== 公開メソッド ====//
// パラメータの内容に従ってLoopLDPC符号を生成する。
Matrix<BinaryFiniteField> LoopLDPC::generate() {
    return mat;
}

// 入力された帯の情報を表示する。
void LoopLDPC::printBandData() {
    int size = bandData.size();
    for(int i = 0; i < size; i++) {
        cout << "[" << i << "] ";
        cout << "l = " << bandData[i].l << ", ";
        cout << "r = " << bandData[i].r << ", ";
        cout << "L = " << bandData[i].L << "\n";
        cout << "\tcbl = " << bandData[i].conBandLeft << ", ";
        cout << "cbr = " << bandData[i].conBandRight << ", ";
        cout << "cpl = " << bandData[i].conPosLeft << ", ";
        cout << "cpr = " << bandData[i].conPosRight << "\n";
    }
}

//==== 非公開メソッド ====//
// 行列のサイズを決定する
unsigned long LoopLDPC::calculateRowSize(const vector<LBandData>& input) {

    int blockNum = input.size();
    unsigned long retVal = 0;
    for(int i = 0; i < blockNum; i++) {
        retVal += input[i].l + input[i].L;
        if(input[i].l <= 0 || input[i].L <= 0) throw ParameterExeption();
    }
    retVal -= blockNum;
    return retVal;
}
unsigned long LoopLDPC::calculateColSize(const vector<LBandData>& input) {
    int blockNum = input.size();
    unsigned long retVal = 0;
    for(int i = 0; i < blockNum; i++) {
        retVal += (input[i].r * input[i].L) / input[i].l;
        if(input[i].r <= 0 || input[i].L <= 0 || input[i].l <= 0) throw ParameterExeption();
    }
    return retVal;
}

// 各帯の行列中の開始位置(行と列)を取得し、bsrpとbscpを更新
void LoopLDPC::calculateBandStartPos() {
    int blockNum = bandData.size();
    unsigned long m = 0;
    unsigned long n = 0;
    for(int i = 0; i < blockNum; i++) {
        bsrp[i] = m;
        bscp[i] = n;
        m += bandData[i].l + bandData[i].L - 1;
        n += bandData[i].r * bandData[i].L / bandData[i].l;
        //cout << "bsrp[" << i << "] = " << bsrp[i] << ", bscp[" << i << "] = " << bscp[i] << "\n";
    }
}

// bandNum番目の帯におけるblockNum番目のブロックの先頭位置(行と列)を取得
unsigned long LoopLDPC::calculateBlockStartRowPos(int bandNum, int blockNum){
    return bsrp[bandNum] + blockNum;
}
unsigned long LoopLDPC::calculateBlockStartColPos(int bandNum, int blockNum) {
    return bscp[bandNum] + bandData[bandNum].r * blockNum / bandData[bandNum].l;
}

// bandNum番目の帯におけるblockNum番目のブロックを生成する
void LoopLDPC::generateBlock(int bandNum, int blockNum) {
    unsigned long startRow = calculateBlockStartRowPos(bandNum, blockNum);
    unsigned long endRow = startRow + bandData[bandNum].l;
    unsigned long startCol = calculateBlockStartColPos(bandNum, blockNum);
    unsigned long endCol = startCol + bandData[bandNum].r / bandData[bandNum].l;
    for(unsigned long i = startRow; i < endRow; i++) {
        for(unsigned long j = startCol; j < endCol; j++) {
            mat[i][j] = 1;
        }
    }
}

// bandNum番目の帯を行列中に配置する
void LoopLDPC::generateBand(int bandNum){
    int L = bandData[bandNum].L;
    for(int i = 0; i < L; i++) {
        generateBlock(bandNum, i);
    }
}

// 全ての空間結合プロトグラフを検査行列中に配置する
void LoopLDPC::generateBaseSCCProtograph() {
    int bandNum = bandData.size();
    for(int i = 0; i < bandNum; i++) {
        generateBand(i);
    }
}

// 結合基準行の算出
unsigned long LoopLDPC::calculateBaseConRowLeft(int bandNum) {
    return bsrp[bandNum];
}
unsigned long LoopLDPC::calculateBaseConRowRight(int bandNum) {
    return bsrp[bandNum] + bandData[bandNum].l + bandData[bandNum].L - 2;
}
// 結合基準列の算出
unsigned long LoopLDPC::calculateBaseConColLeft(int bandNum) {
    int conBandNum = bandData[bandNum].conBandLeft;    // 接続先の帯の番号
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    int conBlockPos = bandData[bandNum].conPosLeft;    // 接続先の帯の何番目のブロックに接続するのか
    if(conBlockPos < 0 || conBlockPos >= bandData[conBandNum].L) throw ParameterExeption();
    //cout << "[left] bandData[" << bandNum << "].conBandLeft = " << conBandNum << ", bandData[" << bandNum << "].conPosLeft = " << conBlockPos << "\n";
    //cout << "[left] bCol = " << bscp[conBandNum] + (bandData[conBandNum].r * conBlockPos) / bandData[conBandNum].l << "\n";
    return bscp[conBandNum] + (bandData[conBandNum].r * conBlockPos) / bandData[conBandNum].l;
}
unsigned long LoopLDPC::calculateBaseConColRight(int bandNum){
    int conBandNum = bandData[bandNum].conBandRight;    // 接続先の帯の番号
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    int conBlockPos = bandData[bandNum].conPosRight;    // 接続先の帯の何番目のブロックに接続するのか
    if(conBlockPos < 0 || conBlockPos >= bandData[conBandNum].L) throw ParameterExeption();
    //cout << "[right] bandData[" << bandNum << "].conBandLeft = " << conBandNum << ", bandData[" << bandNum << "].conPosLeft = " << conBlockPos << "\n";
    return bscp[conBandNum] + (bandData[conBandNum].r * conBlockPos) / bandData[conBandNum].l;
}

// 配置ブロックの単位長を算出
unsigned long LoopLDPC::calculateConBlockLengthLeft(int bandNum) {
    int conBandNum = bandData[bandNum].conBandLeft;    // 接続先の帯の番号
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    return bandData[conBandNum].r / bandData[conBandNum].l;
}
unsigned long LoopLDPC::calculateConBlockLengthRight(int bandNum) {
    int conBandNum = bandData[bandNum].conBandRight;    // 接続先の帯の番号
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    return bandData[conBandNum].r / bandData[conBandNum].l;
}

// stepステップ目における残り配置ビット数を算出
int LoopLDPC::calculateRemainSetBits(int bandNum, int step) {
    return bandData[bandNum].r - ((bandData[bandNum].r / bandData[bandNum].l) * (step + 1));
}

// 配置ブロック数を算出
int LoopLDPC::calculateSetBlockNumLeft(int bandNum, int remainBits) {
    int conBandNum = bandData[bandNum].conBandLeft;    // 接続先の帯の番号
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    return remainBits / (bandData[conBandNum].r / bandData[conBandNum].l);
}
int LoopLDPC::calculateSetBlockNumRight(int bandNum, int remainBits) {
    int conBandNum = bandData[bandNum].conBandRight;    // 接続先の帯の番号
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    return remainBits / (bandData[conBandNum].r / bandData[conBandNum].l);
}

// 結合基準行の更新
unsigned long LoopLDPC::calculateBRowLeft(int bandNum, int step) {
    return calculateBaseConRowLeft(bandNum) + step;
}
unsigned long LoopLDPC::calculateBRowRight(int bandNum, int step) {
    return calculateBaseConRowRight(bandNum) - step;
}

// 指定範囲に結合ブロックを配置
void LoopLDPC::generateConnectBar(long bRow, long bColS, long bColE) {
    //cout << "bRow = " << bRow << ", bColS = " << bColS << ", bColE = " << bColE << "\n";
    if(bRow < 0 || bRow >= mat.row()) throw ParameterExeption();
    //if(bColS < 0 || bColS >= mat.col()) throw ParameterExeption();
    //if(bColE <= 0 || bColE > mat.col()) throw ParameterExeption();
    if(bColS > bColE) throw ParameterExeption();
    for(int i = bColS; i < bColE; i++) {
        mat[bRow][i] = 1;
    }
}

// 結合の反復処理1回分
int LoopLDPC::connectWithStepLeft(int bandNum, int step) {
    int bits = calculateRemainSetBits(bandNum, step);    // 残り配置ビット数を計算
    if(bits <= 0) return -1;

    int conBandNum = bandData[bandNum].conBandLeft;    // 接続先の帯の番号
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    unsigned long bRow = calculateBRowLeft(bandNum, step);    // 基準行の更新
    unsigned long bCol = calculateBaseConColLeft(bandNum);    // 基準列
    //cout << "基準列：" << bCol << "\n";

    int blocks = calculateSetBlockNumLeft(bandNum, bits);    // 配置ブロック数を計算
    if(blocks % 2 == 1) {
        long startColPos = bCol;
        long endColPos = bCol + (bandData[conBandNum].r / bandData[conBandNum].l);
        //cout << "[left/center]startColPos = " << startColPos << ", endColPos = " << endColPos << "\n";
        generateConnectBar(bRow, startColPos, endColPos);
    }

    int putSteps = blocks / 2;
    for(int i = 0; i < putSteps; i++) {
        long startColPos = bCol - (bandData[conBandNum].r / bandData[conBandNum].l) * (i + 1);
        long endColPos = bCol - (bandData[conBandNum].r / bandData[conBandNum].l) * i;
        //cout << "[left/left]startColPos = " << startColPos << ", endColPos = " << endColPos << "\n";
        if(startColPos >= (long)bscp[conBandNum]) {
            //cout << "[left/left]bscp[" << conBandNum << "] = " << bscp[conBandNum] << "\n";
            generateConnectBar(bRow, startColPos, endColPos);
        }

        startColPos = bCol + (bandData[conBandNum].r / bandData[conBandNum].l) * (i + 1);
        endColPos = bCol + (bandData[conBandNum].r / bandData[conBandNum].l) * (i + 2);
        if(endColPos <= (long)bscp[conBandNum] + (bandData[conBandNum].r * bandData[conBandNum].L / bandData[conBandNum].l)) {
            //cout << "[left/right]startColPos = " << startColPos << ", endColPos = " << endColPos << "\n";
            generateConnectBar(bRow, startColPos, endColPos);
        }
    }
    return 0;
}
int LoopLDPC::connectWithStepRight(int bandNum, int step) {
    int bits = calculateRemainSetBits(bandNum, step);    // 残り配置ビット数を計算
    if(bits <= 0) return -1;

    int conBandNum = bandData[bandNum].conBandRight;    // 接続先の帯の番号
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    unsigned long bRow = calculateBRowRight(bandNum, step);    // 基準行の更新
    unsigned long bCol = calculateBaseConColRight(bandNum);    // 基準列

    int blocks = calculateSetBlockNumRight(bandNum, bits);    // 配置ブロック数を計算
    if(blocks % 2 == 1) {
        long startColPos = bCol;
        long endColPos = bCol + (bandData[conBandNum].r / bandData[conBandNum].l);
        //cout << "[right/center]startColPos = " << startColPos << ", endColPos = " << endColPos << "\n";
        generateConnectBar(bRow, startColPos, endColPos);
    }

    int putSteps = blocks / 2;
    for(int i = 0; i < putSteps; i++) {
        long startColPos = bCol - (bandData[conBandNum].r / bandData[conBandNum].l) * (i + 1);
        long endColPos = bCol - (bandData[conBandNum].r / bandData[conBandNum].l) * i;
        //cout << "[rightt/left]startColPos = " << startColPos << ", endColPos = " << endColPos << "\n";
        if(startColPos >= (long)bscp[conBandNum]) {
            generateConnectBar(bRow, startColPos, endColPos);
        }

        startColPos = bCol + (bandData[conBandNum].r / bandData[conBandNum].l) * (i + 1);
        endColPos = bCol + (bandData[conBandNum].r / bandData[conBandNum].l) * (i + 2);
        //cout << "[right/right]startColPos = " << startColPos << ", endColPos = " << endColPos << "\n";
        if(endColPos <= (long)bscp[conBandNum] + (bandData[conBandNum].r * bandData[conBandNum].L / bandData[conBandNum].l)) {
            generateConnectBar(bRow, startColPos, endColPos);
        }
    }
    return 0;
}

// 左右の結合部の作成
void LoopLDPC::generateConnectParts(int bandNum) {
    if(bandData[bandNum].conBandLeft >= 0) {
        int step = 0;
        while(1) {
            int check = connectWithStepLeft(bandNum, step);
            if(check == -1) break;
            step++;
        }
    }
    if(bandData[bandNum].conBandRight >= 0) {
        int step = 0;
        while(1) {
            int check = connectWithStepRight(bandNum, step);
            if(check == -1) break;
            step++;
        }
    }
}

// 結合部全体の作成
void LoopLDPC::generateConnections() {
    int bandNum = bandData.size();
    for(int i = 0; i < bandNum; i++) {
        generateConnectParts(i);
    }
}
