#include "DensityEvolution.h"

long double DensityEvolution::execDESimulation(DESInput simulationParameter) {
    startClock = clock();    // ���Ԍv���J�n

    //==== �ϐ��錾�E�O���� ====//
    string fileName = simulationParameter.checkMatrixFileName;
    unsigned long hr = H->row();
    unsigned long hc = H->col();
    long double codeRate = (long double)(hc - hr) / hc;

    long double leftSuccess = 0.0;    // ���������m��
    long double rightFailure = 1.0;    // �E�����s�m��

    long double zeroBorderProb = simulationParameter.zeroBorderProb;    // 0�Ƃ݂Ȃ��m��
    unsigned long convergenceFailNumber = simulationParameter.convergenceFailNumber;    // �������s��
    unsigned long maxLoop = simulationParameter.binarySearchLoop;    // �񕪒T���̉񐔂̏��
    DEInput dei = {0.0, zeroBorderProb, convergenceFailNumber};
    DEOutput deo;

    long double borderFripProbability;    // 臒l

    printProgress(0.0);
    cout.precision(10);
    for(unsigned long i = 0; i < maxLoop; i++) {
        //==== �o�ߎ��� ====//
        nowClock = clock();
        unsigned long span = nowClock - startClock;
        unsigned int tmpExecTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
        unsigned int tmpExecTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
        unsigned int tmpExecTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;

        cout << "�v���g�O���t�F" << fileName << "(" << hr << " �~ " << hc << ", R = " << codeRate << ")\n";
        cout << "�o�ߎ��ԁF" << tmpExecTimeHour << "����" << tmpExecTimeMin << "��" << tmpExecTimeSec << "�b\n";
        cout << "���s�񐔁F" << i + 1 << " (�ő�F" << maxLoop << ")\n";
        long double centerProb = (leftSuccess + rightFailure) / 2.0;
        cout << "�����F" << leftSuccess << "\n";
        cout << "�E���F" << rightFailure << "\n";
        cout << "�T�����F" << centerProb << "\n";

        dei.fripProbability = centerProb;
        deo = execDE(dei);
        if(deo.failFrag) {
            rightFailure = centerProb;
        } else {
            leftSuccess = centerProb;
        }

        printProgress((i + 1.0) / maxLoop);
    }
    //==== ���s���� ====//
    endClock = clock();
    unsigned long span = endClock - startClock;
    unsigned int execTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
    unsigned int execTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
    unsigned int execTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;
    cout << "�����s��F" << fileName << "(" << hr << " �~ " << hc << ", R = " << codeRate << ")\n";
    cout << "�o�ߎ��ԁF" << execTimeHour << "����" << execTimeMin << "��" << execTimeSec << "�b\n";
    cout << "���s�񐔁F" << maxLoop << " (�ő�F" << maxLoop << ")\n";
    cout << "�����F" << leftSuccess << "\n";
    cout << "�E���F" << rightFailure << "\n";
    cout << "�T�����F" << (leftSuccess + rightFailure) / 2 << "\n";


    borderFripProbability = leftSuccess;
    //cout << "臒l�F" << borderFripProbability << "\n";

    return borderFripProbability;
}

DensityEvolution::DEOutput DensityEvolution::execDE(DEInput parameter) {
    //==== �ϐ��錾�E�O���� ====//
    unsigned long hr = H->row();    // �`�F�b�N�m�[�h��
    unsigned long hc = H->col();    // �ϐ��m�[�h��
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;
    double p = parameter.fripProbability;    // ���]�m��

    unsigned long loop = 0;    // ������
    const unsigned long maxLoop = parameter.convergenceFailNumber;    // �ő唽����(�������s����)[debug]
    unsigned long conv = 0; // ������

    long double border = parameter.borderProb;    // 0�Ƃ݂Ȃ��m��
    vector<BinaryFiniteField> borderChecker(hc);    // �ϐ��m�[�h���b�Z�[�W�����E�l�ȉ����ǂ������i�[
    vector<long double> probList(hc);    // �m�ϐ��m�[�h����̎����m����ێ�����x�N�g��(�ԋp�p)

    DEOutput retObj = {vector<long double>(hc), false, 0};    // �ԋp�\����

    //==== �S�Ă̕ϐ��m�[�h����אڂ���`�F�b�N�m�[�h��p�𑗂� ====//
    for(unsigned long i = 0; i < hc; i++) {
        const vector<int> checkList = (*matCol)[i];    // i�Ԗڂ̕ϐ��m�[�h�Ɛڑ����Ă���`�F�b�N�m�[�h�̔ԍ��̃��X�g
        int size = checkList.size();
        for(int j = 0; j < size; j++) {
            int checkId = checkList[j] - 1;    // alist�`���ł́A1���琔���n�߂�̂�-1���ēY�����x�N�g���ɍ��킹��
            mb.variableToCheck(i, checkId, p);    // i�Ԗڂ̕ϐ��m�[�h����AcheckId�Ԗڂ̃`�F�b�N�m�[�h�Ɋm��p�𑗐M
        }
    }

    for(loop = 0; loop < maxLoop; loop++) {
        //==== �`�F�b�N�m�[�h���� ====//
        for(unsigned long i = 0; i < hr; i++) {
            const vector<int>* variableList = &((*matRow)[i]);
            unsigned long size = variableList->size();

            long double cProb = 1;
            for(unsigned long j = 0; j < size; j++) {
                cProb *= 1 - mb.checkBuf[i][j];
            }

            // �ϐ��m�[�h�ւ̃��b�Z�[�W���Z�o����
            for(unsigned long j = 0; j < size; j++) {
                unsigned long variableId = (*variableList)[j] - 1;
                long double writeValue = 1 - cProb;
                if(1 - mb.getCheckMessage(i, variableId) != 0.0) {
                    writeValue = 1 - (cProb / (1 - mb.getCheckMessage(i, variableId)));
                }
                mb.checkToVariable(i, variableId, writeValue);
            }
        }

        //==== �ϐ��m�[�h���� & �������� ====//
        for(unsigned long i = 0; i < hc; i++) {
            //==== �ϐ��m�[�h���� ====//
            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            long double vProb = p;
            for(unsigned long j = 0; j < size; j++) {
                vProb *= mb.variableBuf[i][j];
            }

            // �`�F�b�N�m�[�h�ւ̃��b�Z�[�W���Z�o����
            for(unsigned long j = 0; j < size; j++) {
                unsigned long checkId = (*checkList)[j] - 1;
                long double writeValue = vProb;
                if(mb.getVariableMessage(i, checkId) != 0.0) {
                    writeValue = vProb / mb.getVariableMessage(i, checkId);
                }
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << "i = " << i << ", vProb = " << vProb << "\n";
            //==== ��������(���E�l�ȉ��Ȃ��1) ====//
            if(vProb < border) {
                borderChecker[i] = 1;
            } else {
                borderChecker[i] = 0;
            }
            probList[i] = vProb;
            //cout << (((int)borderChecker[i] == 1) ? "��" : "�~") << "vProb[" << i << "] = " << vProb << "\n";
        }

        //==== �������� ====//
        // borderChecker�̒l���S��1�Ȃ�Ύ���
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

    //==== �����񐔂̔��� ====//
    conv = loop + 1;

    //==== �����f�[�^�̐��` ====//
    for(unsigned long i = 0; i < hc; i++) {
        retObj.convProbList[i] = probList[i];
    }
    if(loop < maxLoop) {
        retObj.failFrag = false;
        //cout << "��������!!\n";
    } else {
        retObj.failFrag = true;
        //cout << "�������s�E�E�E\n";
    }
    retObj.convergenceNumber = conv;
    return retObj;
}

DensityEvolution::DEWSOutput DensityEvolution::execDEwithStepCheck(DEInput parameter, string fileName, unsigned long defConvNum) {
    //==== �ϐ��錾�E�O���� ====//
    unsigned long hr = H->row();    // �`�F�b�N�m�[�h��
    unsigned long hc = H->col();    // �ϐ��m�[�h��
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;
    double p = parameter.fripProbability;    // ���]�m��

    unsigned long loop = 0;    // ������
    const unsigned long maxLoop = parameter.convergenceFailNumber;    // �ő唽����(�������s����)[debug]
    unsigned long conv = 0; // ������

    long double border = parameter.borderProb;    // 0�Ƃ݂Ȃ��m��
    vector<int> borderChecker(hc);    // �ϐ��m�[�h���b�Z�[�W�����E�l�ȉ����ǂ������i�[
    vector<long double> probList(hc);    // �m�ϐ��m�[�h����̎����m����ێ�����x�N�g��(�ԋp�p)

    unsigned long maxWriteSize = maxLoop / 100 + 100;
    unsigned long writeCount = 0;    // �f�[�^���������
    DEWSOutput retObj = {vector<long double>(hc), false, 0, Matrix<long double>(maxWriteSize, hc), vector<long double>(maxWriteSize), vector<long double>(maxWriteSize)};    // �ԋp�\����

    vector<unsigned long> valueNumChecker(maxWriteSize);    // �ǂ̔������L�^��������ێ�
    Matrix<long double> valueChecker(maxWriteSize, hc);
    vector<long double> maxValueChecker(maxWriteSize);    // ��藦�̍ő�l��ێ�
    vector<long double> minValueChecker(maxWriteSize);    // ��藦�̍ŏ��l��ێ�
    vector<long double> aveValueChecker(maxWriteSize);    // ���ό�藦��ێ�

    for(unsigned long i = 0; i < maxWriteSize; i++) maxValueChecker[i] = -1.0;
    for(unsigned long i = 0; i < maxWriteSize; i++) minValueChecker[i] = -1.0;
    for(unsigned long i = 0; i < maxWriteSize; i++) aveValueChecker[i] = 0.0;

    //==== �S�Ă̕ϐ��m�[�h����אڂ���`�F�b�N�m�[�h��p�𑗂� ====//
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

        //==== �`�F�b�N�m�[�h���� ====//
        for(unsigned long i = 0; i < hr; i++) {
            const vector<int>* variableList = &((*matRow)[i]);
            unsigned long size = variableList->size();

            long double cProb = 1;
            for(unsigned long j = 0; j < size; j++) {
                cProb *= 1 - mb.checkBuf[i][j];
            }

            // �ϐ��m�[�h�ւ̃��b�Z�[�W���Z�o����
            for(unsigned long j = 0; j < size; j++) {
                unsigned long variableId = (*variableList)[j] - 1;
                long double writeValue = 1 - cProb;
                if(1 - mb.getCheckMessage(i, variableId) != 0.0) {
                    writeValue = 1 - (cProb / (1 - mb.getCheckMessage(i, variableId)));
                }
                mb.checkToVariable(i, variableId, writeValue);
            }
        }

        //==== �ϐ��m�[�h���� & �������� ====//
        for(unsigned long i = 0; i < hc; i++) {
            //==== �ϐ��m�[�h���� ====//
            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            long double vProb = p;
            for(unsigned long j = 0; j < size; j++) {
                vProb *= mb.variableBuf[i][j];
            }

            // �`�F�b�N�m�[�h�ւ̃��b�Z�[�W���Z�o����
            for(unsigned long j = 0; j < size; j++) {
                unsigned long checkId = (*checkList)[j] - 1;
                long double writeValue = vProb;
                if(mb.getVariableMessage(i, checkId) != 0.0) {
                    writeValue = vProb / mb.getVariableMessage(i, checkId);
                }
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << "i = " << i << ", vProb = " << vProb << "\n";
            //==== ��������(���E�l�ȉ��Ȃ��1) ====//
            if(vProb < border) {
                borderChecker[i] = 1;
            } else {
                borderChecker[i] = 0;
            }
            probList[i] = vProb;
            //cout << (((int)borderChecker[i] == 1) ? "��" : "�~") << "vProb[" << i << "] = " << vProb << "\n";
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

        //==== �������� ====//
        // borderChecker�̒l���S��1�Ȃ�Ύ���
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

    //==== �����񐔂̔��� ====//
    conv = loop + 1;

    //==== �����f�[�^�̐��` ====//
    for(unsigned long i = 0; i < hc; i++) {
        retObj.convProbList[i] = probList[i];
    }
    if(loop < maxLoop) {
        retObj.failFrag = false;
        cout << "��������!!\n";
    } else {
        retObj.failFrag = true;
        cout << "�������s�E�E�E\n";
    }
    retObj.convergenceNumber = conv;
    retObj.valueChecker = valueChecker;

    cout << "�����񐔁F" << conv << "��\n";
    cout << "�t�@�C�������o���J�n\n";
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
    cout << "�t�@�C�������o������\n";

    return retObj;
}

// ���[�v��Ԍ��������̌`��𗘗p���������̑т�ڍ�������Ԍ��������ɂ��āA��������܂ŉ�DE
// connectionArea : �ڍ����̈ʒu��\���x�N�g��
// blockEndPos : �e�т̏I�[����\���x�N�g��
DensityEvolution::DEOutput DensityEvolution::execSubDE(DEInput parameter, vector<int> connectionArea, vector<int> blockEndPos) {
    //==== �ϐ��錾�E�O���� ====//
    unsigned long hr = H->row();    // �`�F�b�N�m�[�h��
    unsigned long hc = H->col();    // �ϐ��m�[�h��
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;
    double p = parameter.fripProbability;    // ���]�m��

    unsigned long loop = 0;    // ������
    const unsigned long maxLoop = parameter.convergenceFailNumber;    // �ő唽����(�������s����)[debug]
    unsigned long conv = 0; // ������

    long double border = parameter.borderProb;    // 0�Ƃ݂Ȃ��m��
    vector<BinaryFiniteField> borderChecker(hc);    // �ϐ��m�[�h���b�Z�[�W�����E�l�ȉ����ǂ������i�[
    vector<long double> probList(hc);    // �m�ϐ��m�[�h����̎����m����ێ�����x�N�g��(�ԋp�p)

    unsigned long separateNum = 0;    // ��������������

    DEOutput retObj = {vector<long double>(hc), false, 0};    // �ԋp�\����

    //==== �S�Ă̕ϐ��m�[�h����אڂ���`�F�b�N�m�[�h��p�𑗂� ====//
    for(unsigned long i = 0; i < hc; i++) {
        const vector<int> checkList = (*matCol)[i];    // i�Ԗڂ̕ϐ��m�[�h�Ɛڑ����Ă���`�F�b�N�m�[�h�̔ԍ��̃��X�g
        int size = checkList.size();
        for(int j = 0; j < size; j++) {
            int checkId = checkList[j] - 1;    // alist�`���ł́A1���琔���n�߂�̂�-1���ēY�����x�N�g���ɍ��킹��
            mb.variableToCheck(i, checkId, p);    // i�Ԗڂ̕ϐ��m�[�h����AcheckId�Ԗڂ̃`�F�b�N�m�[�h�Ɋm��p�𑗐M
        }
    }

    for(loop = 0; loop < maxLoop; loop++) {
        //==== �`�F�b�N�m�[�h���� ====//
        for(unsigned long i = 0; i < hr; i++) {
            const vector<int>* variableList = &((*matRow)[i]);
            unsigned long size = variableList->size();

            long double cProb = 1;
            for(unsigned long j = 0; j < size; j++) {
                cProb *= 1 - mb.checkBuf[i][j];
            }

            // �ϐ��m�[�h�ւ̃��b�Z�[�W���Z�o����
            for(unsigned long j = 0; j < size; j++) {
                unsigned long variableId = (*variableList)[j] - 1;
                long double writeValue = 1 - cProb;
                if(1 - mb.getCheckMessage(i, variableId) != 0.0) {
                    writeValue = 1 - (cProb / (1 - mb.getCheckMessage(i, variableId)));
                }
                mb.checkToVariable(i, variableId, writeValue);
            }
        }

        //==== �ϐ��m�[�h���� & �������� ====//
        for(unsigned long i = 0; i < hc; i++) {
            //==== �ϐ��m�[�h���� ====//
            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            long double vProb = p;
            for(unsigned long j = 0; j < size; j++) {
                vProb *= mb.variableBuf[i][j];
            }

            // �`�F�b�N�m�[�h�ւ̃��b�Z�[�W���Z�o����
            for(unsigned long j = 0; j < size; j++) {
                unsigned long checkId = (*checkList)[j] - 1;
                long double writeValue = vProb;
                if(mb.getVariableMessage(i, checkId) != 0.0) {
                    writeValue = vProb / mb.getVariableMessage(i, checkId);
                }
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << "i = " << i << ", vProb = " << vProb << "\n";
            //==== ��������(���E�l�ȉ��Ȃ��1) ====//
            if(vProb < border) {
                borderChecker[i] = 1;
            } else {
                borderChecker[i] = 0;
            }
            probList[i] = vProb;
            //cout << (((int)borderChecker[i] == 1) ? "��" : "�~") << "vProb[" << i << "] = " << vProb << "\n";
        }

        //==== �������� ====//
        // connectionArea�̗v�f�ɑΉ�����ϐ��m�[�h���S�Ď������Ă��邩�𔻕�
        bool separateFrag = true;
        int caSize = connectionArea.size();
        for(int i = 0; i < caSize; i++) {
            if((int)borderChecker[connectionArea[i]] == 0) {
                separateFrag = false;
                break;
            }
        }
        if(separateFrag) {
            separateNum = loop;    // ��������������
            cout << loop << "��ڂŕ���!\n";
            break;
        }

        // borderChecker�̒l���S��1�Ȃ�Ύ���
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

    //==== �����񐔂̔��� ====//
    conv = loop + 1;

    //==== �����f�[�^�̐��` ====//
    for(unsigned long i = 0; i < hc; i++) {
        retObj.convProbList[i] = probList[i];
    }
    if(loop < maxLoop) {
        retObj.failFrag = false;
        //cout << "��������!!\n";
    } else {
        retObj.failFrag = true;
        //cout << "�������s�E�E�E\n";
    }
    retObj.convergenceNumber = conv;
    return retObj;
}

// �[��(depth) : ���]�m���̋��E�l�̐��x(1.0 * 10^(-1 * depth)�̐��x�ŋ��߂�)
long double DensityEvolution::borderProb(int valDim, int checkDim, int depth) const {
    long double border = 0.0;    // ���E�l�̌��
    for(int i = 1; i <= depth; i++) {
        int j;    // j * pow(0.1, i)�̒l��border�ɑ����B
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
    bool convergence = false;    // 0�Ɏ���������true
    long double p = fripProbability;
    unsigned long maxLoop = 100000;
    unsigned long loop = 1;
    long double border = pow((long double)0.1, 15);    // 0���Ɣ��肷�鋫�E�l
    for(loop = 1; loop <= maxLoop; loop++) {
        p = deeq(p, fripProbability, checkDim, valDim);
        if(p < border) {
            //cout << "0�Ɏ������܂���(fp = " << fripProbability <<", p = " << p << ", loop = " << loop << ")\n";
            convergence = true;
            break;
        }
    }
    if(loop >= maxLoop) {
        //cout << "0�Ɏ������܂���(fp = " << fripProbability <<", p = " << p << ", loop = " << loop << ")\n";
    }

    if(convergence) return true;
    else return false;
}

//==== �i���󋵂̕\�� ====//
void DensityEvolution::printProgress(double progressRate) const {
    int maxLength = 30;
    int length = (int)(maxLength * progressRate);
    system("cls");
    for(int i = 0; i < length; i++) {
        cout << "��";
    }
    for(int i = length; i < maxLength; i++) {
        cout << "��";
    }
    cout << " (" << progressRate * 100 << "%����)\n";
}

//==== �������ʂ�ۑ� ====//
void DensityEvolution::saveDEData(string fileName, double epsilon, unsigned long loop, unsigned long girth) const {
    ofstream ofs("deResult.txt", fstream::app);
    if(ofs.fail()) cout<< "�t�@�C�����J���܂���B\n";
    else {
        ofs.precision(10);
        ofs << fileName << " " << epsilon << " " << loop << " " << girth << "\n";
    }
    ofs.close();
}
