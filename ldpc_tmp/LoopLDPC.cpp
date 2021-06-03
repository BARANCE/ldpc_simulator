#include "LoopLDPC.h"

//==== �R���X�g���N�^ ====//
LoopLDPC::LoopLDPC(const vector<LBandData>& input) : bandData(input), mat(calculateRowSize(input), calculateColSize(input)), bsrp(input.size()), bscp(input.size()) {
    calculateBandStartPos();
    generateBaseSCCProtograph();
    generateConnections();
    //printBandData();
    //cout << "m = " << mat.row() << ", n = " << mat.col() << "\n";
}

//==== ���J���\�b�h ====//
// �p�����[�^�̓��e�ɏ]����LoopLDPC�����𐶐�����B
Matrix<BinaryFiniteField> LoopLDPC::generate() {
    return mat;
}

// ���͂��ꂽ�т̏���\������B
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

//==== ����J���\�b�h ====//
// �s��̃T�C�Y�����肷��
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

// �e�т̍s�񒆂̊J�n�ʒu(�s�Ɨ�)���擾���Absrp��bscp���X�V
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

// bandNum�Ԗڂ̑тɂ�����blockNum�Ԗڂ̃u���b�N�̐擪�ʒu(�s�Ɨ�)���擾
unsigned long LoopLDPC::calculateBlockStartRowPos(int bandNum, int blockNum){
    return bsrp[bandNum] + blockNum;
}
unsigned long LoopLDPC::calculateBlockStartColPos(int bandNum, int blockNum) {
    return bscp[bandNum] + bandData[bandNum].r * blockNum / bandData[bandNum].l;
}

// bandNum�Ԗڂ̑тɂ�����blockNum�Ԗڂ̃u���b�N�𐶐�����
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

// bandNum�Ԗڂ̑т��s�񒆂ɔz�u����
void LoopLDPC::generateBand(int bandNum){
    int L = bandData[bandNum].L;
    for(int i = 0; i < L; i++) {
        generateBlock(bandNum, i);
    }
}

// �S�Ă̋�Ԍ����v���g�O���t�������s�񒆂ɔz�u����
void LoopLDPC::generateBaseSCCProtograph() {
    int bandNum = bandData.size();
    for(int i = 0; i < bandNum; i++) {
        generateBand(i);
    }
}

// ������s�̎Z�o
unsigned long LoopLDPC::calculateBaseConRowLeft(int bandNum) {
    return bsrp[bandNum];
}
unsigned long LoopLDPC::calculateBaseConRowRight(int bandNum) {
    return bsrp[bandNum] + bandData[bandNum].l + bandData[bandNum].L - 2;
}
// �������̎Z�o
unsigned long LoopLDPC::calculateBaseConColLeft(int bandNum) {
    int conBandNum = bandData[bandNum].conBandLeft;    // �ڑ���̑т̔ԍ�
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    int conBlockPos = bandData[bandNum].conPosLeft;    // �ڑ���̑т̉��Ԗڂ̃u���b�N�ɐڑ�����̂�
    if(conBlockPos < 0 || conBlockPos >= bandData[conBandNum].L) throw ParameterExeption();
    //cout << "[left] bandData[" << bandNum << "].conBandLeft = " << conBandNum << ", bandData[" << bandNum << "].conPosLeft = " << conBlockPos << "\n";
    //cout << "[left] bCol = " << bscp[conBandNum] + (bandData[conBandNum].r * conBlockPos) / bandData[conBandNum].l << "\n";
    return bscp[conBandNum] + (bandData[conBandNum].r * conBlockPos) / bandData[conBandNum].l;
}
unsigned long LoopLDPC::calculateBaseConColRight(int bandNum){
    int conBandNum = bandData[bandNum].conBandRight;    // �ڑ���̑т̔ԍ�
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    int conBlockPos = bandData[bandNum].conPosRight;    // �ڑ���̑т̉��Ԗڂ̃u���b�N�ɐڑ�����̂�
    if(conBlockPos < 0 || conBlockPos >= bandData[conBandNum].L) throw ParameterExeption();
    //cout << "[right] bandData[" << bandNum << "].conBandLeft = " << conBandNum << ", bandData[" << bandNum << "].conPosLeft = " << conBlockPos << "\n";
    return bscp[conBandNum] + (bandData[conBandNum].r * conBlockPos) / bandData[conBandNum].l;
}

// �z�u�u���b�N�̒P�ʒ����Z�o
unsigned long LoopLDPC::calculateConBlockLengthLeft(int bandNum) {
    int conBandNum = bandData[bandNum].conBandLeft;    // �ڑ���̑т̔ԍ�
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    return bandData[conBandNum].r / bandData[conBandNum].l;
}
unsigned long LoopLDPC::calculateConBlockLengthRight(int bandNum) {
    int conBandNum = bandData[bandNum].conBandRight;    // �ڑ���̑т̔ԍ�
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    return bandData[conBandNum].r / bandData[conBandNum].l;
}

// step�X�e�b�v�ڂɂ�����c��z�u�r�b�g�����Z�o
int LoopLDPC::calculateRemainSetBits(int bandNum, int step) {
    return bandData[bandNum].r - ((bandData[bandNum].r / bandData[bandNum].l) * (step + 1));
}

// �z�u�u���b�N�����Z�o
int LoopLDPC::calculateSetBlockNumLeft(int bandNum, int remainBits) {
    int conBandNum = bandData[bandNum].conBandLeft;    // �ڑ���̑т̔ԍ�
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    return remainBits / (bandData[conBandNum].r / bandData[conBandNum].l);
}
int LoopLDPC::calculateSetBlockNumRight(int bandNum, int remainBits) {
    int conBandNum = bandData[bandNum].conBandRight;    // �ڑ���̑т̔ԍ�
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    return remainBits / (bandData[conBandNum].r / bandData[conBandNum].l);
}

// ������s�̍X�V
unsigned long LoopLDPC::calculateBRowLeft(int bandNum, int step) {
    return calculateBaseConRowLeft(bandNum) + step;
}
unsigned long LoopLDPC::calculateBRowRight(int bandNum, int step) {
    return calculateBaseConRowRight(bandNum) - step;
}

// �w��͈͂Ɍ����u���b�N��z�u
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

// �����̔�������1��
int LoopLDPC::connectWithStepLeft(int bandNum, int step) {
    int bits = calculateRemainSetBits(bandNum, step);    // �c��z�u�r�b�g�����v�Z
    if(bits <= 0) return -1;

    int conBandNum = bandData[bandNum].conBandLeft;    // �ڑ���̑т̔ԍ�
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    unsigned long bRow = calculateBRowLeft(bandNum, step);    // ��s�̍X�V
    unsigned long bCol = calculateBaseConColLeft(bandNum);    // ���
    //cout << "���F" << bCol << "\n";

    int blocks = calculateSetBlockNumLeft(bandNum, bits);    // �z�u�u���b�N�����v�Z
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
    int bits = calculateRemainSetBits(bandNum, step);    // �c��z�u�r�b�g�����v�Z
    if(bits <= 0) return -1;

    int conBandNum = bandData[bandNum].conBandRight;    // �ڑ���̑т̔ԍ�
    if(conBandNum < 0 || conBandNum >= bandData.size()) throw ParameterExeption();
    unsigned long bRow = calculateBRowRight(bandNum, step);    // ��s�̍X�V
    unsigned long bCol = calculateBaseConColRight(bandNum);    // ���

    int blocks = calculateSetBlockNumRight(bandNum, bits);    // �z�u�u���b�N�����v�Z
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

// ���E�̌������̍쐬
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

// �������S�̂̍쐬
void LoopLDPC::generateConnections() {
    int bandNum = bandData.size();
    for(int i = 0; i < bandNum; i++) {
        generateConnectParts(i);
    }
}
