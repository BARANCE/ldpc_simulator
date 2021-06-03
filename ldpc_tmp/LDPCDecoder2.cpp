#include "LDPCDecoder2.h"
#include <math.h>

//==== �R���X�g���N�^ ====//
LDPCDecoder2::LDPCDecoder2(const SPMatrix* checkMatrix) : H(checkMatrix), mb(checkMatrix), channel(C_AWGN), SNR(0.0), fp(0.0) {
    //ESPInput espi = {vector<double>(4000, -1.0), 100};
    //sumProductDecode(espi);
}

//==== �ʐM�H�̑I�� ====//
void LDPCDecoder2::setChannel(EChannel channelName) {
    channel = channelName;
}

//==== �ʐM�H�Ɋւ���p�����[�^�̐ݒ� ====//
void LDPCDecoder2::setSNR(double rate) {
    SNR = rate;
}
void LDPCDecoder2::setFlipProbability(double flipProbability) {
    fp = flipProbability;
}

//==== ���� ====//
// �d����
LDPCDecoder2::ESPOutput LDPCDecoder2::hardDecision(ESPInput decodeParameter) {
    //==== �ϐ��錾�E�O���� ====//
    unsigned long hc = H->col();
    const vector<double> receptionWord = decodeParameter.receptionWord;    // ��M��
    vector<BinaryFiniteField> nHat(hc);    // �����

    for(int i = 0; i < hc; i++) {
        if(receptionWord[i] > 0.0) {
            nHat[i] = 0;
        } else if(receptionWord[i] < 0.0) {
            nHat[i] = 1;
        } else {    // �^�C�u���[�N
            double rnd = genrand_real2();
            if(rnd < 0.5) nHat[i] = 0;
            else nHat[i] = 1;
        }
    }

    //==== �p���e�B���� ====//
    bool failFrag = true;    // ���������s������true
    if(parityCheck(&nHat)) {
        failFrag = false;
    }

    //==== �����f�[�^�̐��` ====//
    ESPOutput retObj = {vector<BinaryFiniteField>(hc), false, 0, DECODE_SUCCESS};    // �ԋp�\����
    for(int i = 0; i < hc; i++) {
        retObj.decodeWord[i] = nHat[i];
    }
    retObj.convergenceNumber = 0;
    if(failFrag) {
        retObj.failFrag = true;
        retObj.errorType = OVER_MAX_CONVERGENCE;    // �������s
    } else {
        retObj.failFrag = false;
        retObj.errorType = DECODE_SUCCESS;    // ��������(���������ł͂Ȃ�)
    }

    return retObj;
}

// sum-product����
LDPCDecoder2::ESPOutput LDPCDecoder2::sumProductDecode(ESPInput decodeParameter) {
    //==== �ϐ��錾�E�O���� ====//
    unsigned long hr = H->row();    // �`�F�b�N�m�[�h��
    unsigned long hc = H->col();    // �ϐ��m�[�h��
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;

    unsigned long loop = 0;    // ������
    const unsigned long maxLoop = decodeParameter.convergenceFailNumber;    // �ő唽����(�������s����)
    unsigned long conv = 0;    // ������
    const vector<double> receptionWord = decodeParameter.receptionWord;    // ��M��

    const double codeRate = (double)(hc - hr) / hc;    // ��������
    const double variance = 0.5 * (1.0 / pow(10.0, SNR / 10.0)) * (double)hc / (double)(hc - hr);    // SN�䂩�狁�߂����U�l
    vector<double> logLikeRate(hc);    // �ϐ��m�[�h���Ƃ̑ΐ��ޓx��

    vector<double> renewChecker(hc, -100.0);    // �X�V���i��ł��邩�𒲍�����x�N�g��
    bool badConvergenceFrag = false;

    vector<BinaryFiniteField> nHat(hc);    // �ꎞ�����(�p���e�B������ʉ߂�����ԋp�\���̂Ɋi�[)

    ESPOutput retObj = {vector<BinaryFiniteField>(hc), false, 0, DECODE_SUCCESS};    // �ԋp�\����

    //==== �S�ϐ��m�[�h�̑ΐ��ޓx����v�Z ====//
    for(int i = 0; i < hc; i++) {
        logLikeRate[i] = 2.0 * receptionWord[i] / variance;
    }

    //==== �S�Ẵ`�F�b�N�m�[�h����אڂ���ϐ��m�[�h��0�𑗂� ====//
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

        //cout << "====�ϐ��m�[�h����====\n";

        //==== �ϐ��m�[�h���� & �ꎞ����r�b�g�̌��� ====//
        unsigned long renewCounter = 0;
        for(int i = 0; i < hc; i++) {
            double llr = logLikeRate[i];    // �ΐ��ޓx��
            //double vSum = sumMessageFromCheckNode(i);    // ���b�Z�[�W�̘a
            //double sendTmp = llr + vSum;

            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            double vSum = llr;
            for(int j = 0; j < size; j++) {
                vSum += mb.variableBuf[i][j];
            }

            for(int j = 0; j < size; j++) {    // �`�F�b�N�m�[�h�ւ̃��b�Z�[�W���Z�o����
                unsigned long checkId = (*checkList)[j] - 1;
                double writeValue = vSum - mb.getVariableMessage(i, checkId);
                //cout << "x" << i << "�� c" << checkId << " : " << writeValue << "\n";
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << sendTmp << " ";

            //==== ����r�b�g�̌��� ====//
            if(vSum < 0.0) {
                nHat[i] = 1;
            } else {
                nHat[i] = 0;
            }

            //==== �X�V���i��ł��邩 ====//
            if(renewChecker[i] == vSum) {
                renewCounter++;
            }
            renewChecker[i] = vSum;
        }

        //==== �p���e�B���� ====//
        /*for(int i = 0; i < hc; i++) {
            cout << nHat[i] << " ";
        }
        cout << "\n";*/
        if(parityCheck(&nHat)) {
        //if(isZero(&nHat)) {    // �S�[��������ɂ����ʗp���Ȃ�
            //cout << "��������!!\n";
            break;
        }
        /*if(renewCounter >= hr) {    // �����
            badConvergenceFrag = true;
            break;
        }*/
        // �������s�̓��C�����[�v�̃��[�v�񐔂Ŕ���

        //cout << "====�`�F�b�N�m�[�h����====\n";

        //==== �`�F�b�N�m�[�h���� ====//
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

            for(int j = 0; j < size; j++) {    // �ϐ��m�[�h�ւ̃��b�Z�[�W���Z�o����
                unsigned long variableId = (*variableList)[j] - 1;    // ���M��ϐ��m�[�h��ID

                //cout << "weiteValue = " << writeValue << ", writeValue2 = " << writeValue2 << "\n";
                double tmp2 = mb.getCheckMessage(i, variableId);    // ���M�悩��̃f�[�^
                double signValue = sign(tmp2);
                if(tmp2 != 0.0) {
                    signValue = cProd / signValue;
                } else {
                    signValue = cProd;
                }

                double fValue = gallagerf(absolute(tmp2));
                fValue = gallagerf(cSum - fValue);

                double writeValue = signValue * fValue;

                //cout << "c" << i << "�� x" << variableId << " : " << writeValue << "\n";
                mb.checkToVariable(i, variableId, writeValue);
            }
        }



        //break;

    }
    if(loop < maxLoop) conv = loop + 1;    // ������
    else conv = loop;

    //==== �����f�[�^�̐��` ====//
    for(int i = 0; i < hc; i++) {
        retObj.decodeWord[i] = nHat[i];
    }
    retObj.convergenceNumber = conv;
    if(conv >= maxLoop) {
        retObj.failFrag = true;
        retObj.errorType = OVER_MAX_CONVERGENCE;    // �������s
    } else {
        if(badConvergenceFrag) {
            retObj.failFrag = true;
            retObj.errorType = BAD_CONVERGENCE;    // �����
        } else {
            retObj.failFrag = false;
            retObj.errorType = DECODE_SUCCESS;    // ��������(���������ł͂Ȃ�)
        }
    }

    //printVector(&retObj.decodeWord);

    return retObj;
}

// min-sum����
LDPCDecoder2::ESPOutput LDPCDecoder2::minSumDecode(ESPInput decodeParameter, double scale) {
    //==== �ϐ��錾�E�O���� ====//
    unsigned long hr = H->row();    // �`�F�b�N�m�[�h��
    unsigned long hc = H->col();    // �ϐ��m�[�h��
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;

    unsigned long loop = 0;    // ������
    const unsigned long maxLoop = decodeParameter.convergenceFailNumber;    // �ő唽����(�������s����)
    unsigned long conv = 0;    // ������
    const vector<double> receptionWord = decodeParameter.receptionWord;    // ��M��

    const double codeRate = (double)(hc - hr) / hc;    // ��������
    const double variance = 0.5 * (1.0 / pow(10.0, SNR / 10.0)) * (double)hc / (double)(hc - hr);    // SN�䂩�狁�߂����U�l
    vector<double> logLikeRate(hc);    // �ϐ��m�[�h���Ƃ̑ΐ��ޓx��

    vector<double> renewChecker(hc, -100.0);    // �X�V���i��ł��邩�𒲍�����x�N�g��
    bool badConvergenceFrag = false;

    vector<BinaryFiniteField> nHat(hc);    // �ꎞ�����(�p���e�B������ʉ߂�����ԋp�\���̂Ɋi�[)

    ESPOutput retObj = {vector<BinaryFiniteField>(hc), false, 0, DECODE_SUCCESS};    // �ԋp�\����

    //----�e�X�g�p��������----
    //Matrix<double> hatCheck(maxLoop, hc);    // ����r�b�g�̕ϑJ���L�^����s��
    //----�e�X�g�p�����܂�----

    //==== �S�ϐ��m�[�h�̑ΐ��ޓx����v�Z ====//
    for(int i = 0; i < hc; i++) {
        logLikeRate[i] = 2.0 * receptionWord[i] / variance;
    }

    //==== �S�Ẵ`�F�b�N�m�[�h����אڂ���ϐ��m�[�h��0�𑗂� ====//
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

        //cout << "====�ϐ��m�[�h����====\n";

        //==== �ϐ��m�[�h���� & �ꎞ����r�b�g�̌��� ====//
        unsigned long renewCounter = 0;
        for(int i = 0; i < hc; i++) {
            double llr = logLikeRate[i];    // �ΐ��ޓx��
            //double vSum = sumMessageFromCheckNode(i);    // ���b�Z�[�W�̘a
            //double sendTmp = llr + vSum;

            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            double vSum = llr;
            for(int j = 0; j < size; j++) {
                vSum += mb.variableBuf[i][j];
            }

            for(int j = 0; j < size; j++) {    // �`�F�b�N�m�[�h�ւ̃��b�Z�[�W���Z�o����
                unsigned long checkId = (*checkList)[j] - 1;
                double writeValue = vSum - mb.getVariableMessage(i, checkId);
                //cout << "x" << i << "�� c" << checkId << " : " << writeValue << "\n";
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << sendTmp << " ";

            //---- �e�X�g�p�������� ----//
            //hatCheck[loop][i] = vSum;
            //---- �e�X�g�p�����܂� ----//
            if(vSum != vSum) {
                cout << "NaN���������Ă��܂��B\n";
            }

            //==== ����r�b�g�̌��� ====//
            if(vSum < 0.0) {
                nHat[i] = 1;
            } else if(vSum > 0.0) {
                nHat[i] = 0;
            } else {
                double rnd = genrand_real2();
                if(rnd < 0.5) nHat[i] = 0;
                else nHat[i] = 1;
            }

            //==== �X�V���i��ł��邩 ====//
            if(renewChecker[i] == vSum) {
                renewCounter++;
            }
            renewChecker[i] = vSum;
        }

        //==== �p���e�B���� ====//
        /*for(int i = 0; i < hc; i++) {
            cout << nHat[i] << " ";
        }
        cout << "\n";*/
        if(parityCheck(&nHat)) {
        //if(isZero(&nHat)) {    // �S�[��������ɂ����ʗp���Ȃ�
            //cout << "��������!!\n";
            break;
        }
        /*if(renewCounter >= hr) {    // �����
            badConvergenceFrag = true;
            break;
        }*/
        // �������s�̓��C�����[�v�̃��[�v�񐔂Ŕ���

        //cout << "====�`�F�b�N�m�[�h����====\n";

        //==== �`�F�b�N�m�[�h���� ====//
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

            for(int j = 0; j < size; j++) {    // �ϐ��m�[�h�ւ̃��b�Z�[�W���Z�o����
                unsigned long variableId = (*variableList)[j] - 1;    // ���M��ϐ��m�[�h��ID

                //cout << "weiteValue = " << writeValue << ", writeValue2 = " << writeValue2 << "\n";
                double tmp2 = mb.getCheckMessage(i, variableId);    // ���M�悩��̃f�[�^
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
                    if((*variableList)[k] - 1 != variableId) {    // ���鑊�������
                        if(absolute(mb.checkBuf[i][k]) < cMin) {
                            cMin = absolute(mb.checkBuf[i][k]);
                        }
                    }
                }

                double writeValue = scale * signValue * cMin;

                //cout << "c" << i << "�� x" << variableId << " : " << writeValue << "\n";
                mb.checkToVariable(i, variableId, writeValue);
            }
        }



        //break;

    }
    if(loop < maxLoop) conv = loop + 1;    // ������
    else conv = loop;

    //==== �����f�[�^�̐��` ====//
    for(int i = 0; i < hc; i++) {
        retObj.decodeWord[i] = nHat[i];
    }
    retObj.convergenceNumber = conv;
    if(conv >= maxLoop) {
        retObj.failFrag = true;
        retObj.errorType = OVER_MAX_CONVERGENCE;    // �������s
    } else {
        if(badConvergenceFrag) {
            retObj.failFrag = true;
            retObj.errorType = BAD_CONVERGENCE;    // �����
        } else {
            retObj.failFrag = false;
            retObj.errorType = DECODE_SUCCESS;    // ��������(���������ł͂Ȃ�)
        }
    }

    //---- �e�X�g�p�������� ----//
    //hatCheck.exportMatrix("hatCheck_RPEGSCCm400_n4000_wc4_b40.txt");
    //---- �e�X�g�p�����܂� ----//


    //printVector(&retObj.decodeWord);

    return retObj;
}

// ��������sum-product�����@
LDPCDecoder2::ERSPOutput LDPCDecoder2::sumProductDecodeOnBEC(ESPInput decodeParameter, ECDList researchParameter) {
    //==== �ϐ��錾�E�O���� ====//
    unsigned long hr = H->row();    // �`�F�b�N�m�[�h��
    unsigned long hc = H->col();    // �ϐ��m�[�h��
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;

    unsigned long loop = 0;    // ������
    const unsigned long maxLoop = decodeParameter.convergenceFailNumber;    // �ő唽����(�������s����)
    unsigned long conv = 0;    // ������
    const vector<double> receptionWord = decodeParameter.receptionWord;    // ��M��

    bool badConvergenceFrag = false;

    vector<double> nHat(hc);    // �ꎞ�����(�p���e�B������ʉ߂�����ԋp�\���̂Ɋi�[)

    int convCheckFrag = researchParameter.convCheckFrag;    // �����`�F�b�N�����s���邩�ǂ���
    int convCheckSize = researchParameter.convCheckSize;    // �����`�F�b�N���s�����߂̃o�b�t�@�̃T�C�Y
    vector<vector<double> > convCheck(convCheckSize, vector<double>(hc, 0.0));    // �����󋵂��`�F�b�N����2�����x�N�g��

    ERSPOutput retObj = {vector<double>(hc), false, 0, DECODE_SUCCESS, vector<vector<double> >(0)};    // �ԋp�\����

    //==== ��������l�͎�M�� ====//
    for(int i = 0; i < hc; i++) {
        nHat[i] = receptionWord[i];
        if(convCheckFrag) {
            convCheck[0][i] = receptionWord[i];
        }
    }

    //==== �S�Ẵ`�F�b�N�m�[�h����אڂ���ϐ��m�[�h��e�𑗂� ====//
    for(int i = 0; i < hr; i++) {
        const vector<int> variableList = (*matRow)[i];
        int size = variableList.size();
        for(int j = 0; j < size; j++) {
            int variableId = variableList[j] - 1;
            mb.checkToVariable(i, variableId, ERASURE);
        }
    }

    for(loop = 0; loop < maxLoop; loop++) {

        //==== �ϐ��m�[�h���� & �ꎞ����r�b�g�̌��� ====//
        unsigned long renewCounter = 0;
        for(int i = 0; i < hc; i++) {
            const vector<int>* checkList = &((*matCol)[i]);    // i�Ԗڂ̕ϐ��m�[�h�ɐڑ����Ă���`�F�b�N�m�[�h�ԍ��̃��X�g
            unsigned long size = checkList->size();    // i�Ԗڂ̗�̕ϐ��m�[�h�̎���(��d��)

            if(nHat[i] == ONE || nHat[i] == ZERO) {    // �ϐ��m�[�h�̐���l��0�܂���1�̏ꍇ
                // �אڂ���S�Ẵ`�F�b�N�m�[�h�ɐ���l�𑗐M����
                for(int j = 0; j < size; j++) {
                    unsigned long checkId = (*checkList)[j] - 1;
                    double writeValue = nHat[i];
                    mb.variableToCheck(i, checkId, writeValue);
                }

            } else if(nHat[i] == ERASURE) {    // �ϐ��m�[�h�̐���l��e�̏ꍇ
                // �����Ă����`�F�b�N�m�[�h���ϐ��m�[�h�̃��b�Z�[�W�ɒl�����邩���J�E���g����
                int valueCounter = 0;    // 0�܂���1�̃��b�Z�[�W�̐����J�E���g����
                long messagePointer = -1;    // 0�܂���1�̃��b�Z�[�W�������Ă������b�Z�[�W�̑��葊��(�������-1)
                double value = ERASURE;    // ����l

                for(int j = 0; j < size; j++) {
                    unsigned long checkId = (*checkList)[j] - 1;
                    double message = mb.getVariableMessage(i, checkId);
                    if(message == ZERO || message == ONE) {
                        valueCounter++;
                        messagePointer = checkId;
                        value = message;
                        if(valueCounter >= 2) break;    // ����ȏ���K�v�Ȃ�
                    } else if(message == ERASURE) {
                        // �������Ȃ�
                    } else {
                        cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]���b�Z�[�W�l�ɁA0�E1�Ee��3��ވȊO�̒l���g�p����Ă��܂��B\n";
                    }
                }

                // �J�E���^�[�̒l�ɂ���ď�������
                if(valueCounter == 0) {    // �S��ERASURE
                    // ���M���b�Z�[�W�͑S��ERASURE
                    for(int j = 0; j < size; j++) {
                        unsigned long checkId = (*checkList)[j] - 1;
                        double writeValue = ERASURE;
                        mb.variableToCheck(i, checkId, writeValue);
                    }
                } else if(valueCounter >= 2) {
                    // ���M���b�Z�[�W�͐���l
                    for(int j = 0; j < size; j++) {
                        unsigned long checkId = (*checkList)[j] - 1;
                        double writeValue = value;
                        mb.variableToCheck(i, checkId, writeValue);
                    }
                } else if(valueCounter == 1) {
                    // 0,1�̒l�𑗂��Ă����`�F�b�N�m�[�h�ɂ�ERASURE���A����ȊO��ERASURE�𑗂��Ă����`�F�b�N�m�[�h�ɂ͐���l�𑗂�
                    if(messagePointer == -1) cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]�����V���{���ȊO�̒l�𑗂��Ă����`�F�b�N�m�[�h�̒l���L�^����Ă��܂���B\n";
                    for(int j = 0; j < size; j++) {
                        unsigned long checkId = (*checkList)[j] - 1;
                        double writeValue = value;
                        if(checkId == messagePointer) {
                            writeValue = ERASURE;
                        }
                        mb.variableToCheck(i, checkId, writeValue);
                    }
                } else {
                    cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]�J�E���^�[�̒l���ُ�ł��B\n";
                    throw this;
                }

                //==== ����l�̍X�V ====//
                if(value != ERASURE) nHat[i] = value;
            } else {
                cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]����l�ɁA0�E1�Ee��3��ވȊO�̒l���g�p����Ă��܂��B\n";
                throw this;
            }
        }

        //==== �p���e�B���� ====//
        int parityCounter = 0;
        for(int i = 0; i <hc; i++) {
            if(convCheckFrag && loop < convCheckSize - 1) {
                convCheck[loop + 1][i] = nHat[i];
                if(convCheck[loop][i] != ERASURE && convCheck[loop + 1][i] == ERASURE) {
                    cout << "[loop : " << loop << ", i = " << i << "] ����������\n";
                }
            }
            if(nHat[i] == ERASURE) {
                parityCounter++;
                //break;
            }
        }
        if(parityCounter == 0) break;    // ��������

        // �f�o�b�O�p
        /*if(loop == 0) {
        for(int i = 0; i < hc; i++) {
            if(nHat[i] == ZERO) cout << "0";
            else if(nHat[i] == ONE) cout << "1";
            else if(nHat[i] == ERASURE) cout << "e";
            else cout << "?";
        }
        cout << "\n";
        }*/

        //==== �`�F�b�N�m�[�h���� ====//
        for(int i = 0; i < hr; i++) {
            const vector<int>* variableList = &((*matRow)[i]);
            unsigned long size = variableList->size();

            // �����Ă����ϐ��m�[�h���`�F�b�N�m�[�h�̃��b�Z�[�W�ɏ����V���{���������݂��邩���J�E���g����
            int erasureCounter = 0;
            long messagePointer = -1;    // erasure�𑗂��Ă����ϐ��m�[�h�̈ʒu(�������-1)
            double cSum = ZERO;    // �p���e�B�`�F�b�N�̍��v�l(���b�Z�[�W�Z�o�Ɏg�p)

            for(int j = 0; j < size; j++) {
                unsigned long variableId = (*variableList)[j] - 1;
                double message = mb.getCheckMessage(i, variableId);
                if(message == ZERO || message == ONE) {
                    if(message == ONE) {
                        if(cSum == ZERO) cSum = ONE;
                        else if(cSum == ONE) cSum = ZERO;
                        else cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]���b�Z�[�W���v�l�ɁA0�E1��2��ވȊO�̒l���g�p����Ă��܂��B\n";
                    }
                } else if(message == ERASURE) {
                    erasureCounter++;
                    messagePointer = variableId;
                } else {
                    cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]���b�Z�[�W�l�ɁA0�E1�Ee��3��ވȊO�̒l���g�p����Ă��܂��B\n";
                }
            }

            // �J�E���^�[�̒l�ɂ���ď�������
            if(erasureCounter == 0) {    // �S��0��1
                // ���M���b�Z�[�W�́A���M���肩��̃��b�Z�[�W����������M���b�Z�[�W�̘a
                for(int j = 0; j < size; j++) {
                    unsigned long variableId = (*variableList)[j] - 1;
                    double writeValue = ZERO;
                    if(mb.getCheckMessage(i, variableId) == ONE) {    // writeValue = cSum + mb.getCheckMessage(i, variableId)����낤�Ƃ��Ă��镔��
                        if(cSum == ONE) writeValue = ZERO;
                        else if(cSum == ZERO) writeValue = ONE;
                        else cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]���b�Z�[�W���v�l�ɁA0�E1��2��ވȊO�̒l���g�p����Ă��܂��B\n";
                    } else if(mb.getCheckMessage(i, variableId) == ZERO) {
                        if(cSum == ONE) writeValue = ONE;
                        else if(cSum == ZERO) writeValue = ZERO;
                        else cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]���b�Z�[�W���v�l�ɁA0�E1��2��ވȊO�̒l���g�p����Ă��܂��B\n";
                    } else {
                        cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]�����ł̓��b�Z�[�W��e�ł��邱�Ƃ�z�肵�Ă��܂���B\n";
                    }
                    mb.checkToVariable(i, variableId, writeValue);
                }
            } else if(erasureCounter >= 2) {
                // ���M���b�Z�[�W�͑S��ERASURE
                for(int j = 0; j < size; j++) {
                    unsigned long variableId = (*variableList)[j] - 1;
                    double writeValue = ERASURE;
                    mb.checkToVariable(i, variableId, writeValue);
                }
            } else if(erasureCounter == 1) {
                // ERASURE�𑗂��Ă����ϐ��m�[�h�ɂ͎�M���b�Z�[�W�̘a���A�����łȂ��ϐ��m�[�h�ɂ�ERASURE�𑗂�
                for(int j = 0; j < size; j++) {
                    unsigned long variableId = (*variableList)[j] - 1;
                    double writeValue = ERASURE;
                    if(variableId == messagePointer) {
                        writeValue = cSum;
                    }
                    mb.checkToVariable(i, variableId, writeValue);
                }
            } else {
                cout << "[LDPCDecoder2::sumProductDecodeOnBEC][error]�J�E���^�[�̒l���ُ�ł��B\n";
            }

        }


    }

    if(loop < maxLoop) conv = loop + 1;    // ������
    else conv = loop;

    //==== �����f�[�^�̐��` ====//
    for(int i = 0; i < hc; i++) {
        retObj.decodeWord[i] = nHat[i];
    }
    retObj.convergenceNumber = conv;
    if(conv >= maxLoop) {
        retObj.failFrag = true;
        retObj.errorType = OVER_MAX_CONVERGENCE;    // �������s
    } else {
        if(badConvergenceFrag) {
            retObj.failFrag = true;
            retObj.errorType = BAD_CONVERGENCE;    // �����
        } else {
            retObj.failFrag = false;
            retObj.errorType = DECODE_SUCCESS;    // ��������(���������ł͂Ȃ�)
        }
    }
    retObj.convCheckVec = convCheck;

    return retObj;
}

// �C��min-sum����
// assumptionVector : �V���[�g�������ӏ��̒l���i�[����x�N�g��
// assumptionVector[i] > 0 : i�Ԗڂ̕ϐ��m�[�h��0������
// assumptionVector[i] < 0 : i�Ԗڂ̕ϐ��m�[�h��1������
// assumptionVector[i] = 0 : �V���[�g�����Ȃ��ϐ��m�[�h
LDPCDecoder2::ESPOutput LDPCDecoder2::modMinSumDecode(ESPInput decodeParameter, double scale, vector<int>* assumptionVector) {
    //==== �ϐ��錾�E�O���� ====//
    unsigned long hr = H->row();    // �`�F�b�N�m�[�h��
    unsigned long hc = H->col();    // �ϐ��m�[�h��
    const vector<vector<int> >* matRow = H->matRow;
    const vector<vector<int> >* matCol = H->matCol;

    unsigned long loop = 0;    // ������
    const unsigned long maxLoop = decodeParameter.convergenceFailNumber;    // �ő唽����(�������s����)
    unsigned long conv = 0;    // ������
    const vector<double> receptionWord = decodeParameter.receptionWord;    // ��M��

    const double codeRate = (double)(hc - hr) / hc;    // ��������
    const double variance = 0.5 * (1.0 / pow(10.0, SNR / 10.0)) * (double)hc / (double)(hc - hr);    // SN�䂩�狁�߂����U�l
    vector<double> logLikeRate(hc);    // �ϐ��m�[�h���Ƃ̑ΐ��ޓx��

    vector<double> renewChecker(hc, -100.0);    // �X�V���i��ł��邩�𒲍�����x�N�g��
    bool badConvergenceFrag = false;

    vector<BinaryFiniteField> nHat(hc);    // �ꎞ�����(�p���e�B������ʉ߂�����ԋp�\���̂Ɋi�[)

    ESPOutput retObj = {vector<BinaryFiniteField>(hc), false, 0, DECODE_SUCCESS};    // �ԋp�\����

    //----�e�X�g�p��������----
    //Matrix<double> hatCheck(maxLoop, hc);    // ����r�b�g�̕ϑJ���L�^����s��
    //----�e�X�g�p�����܂�----

    //==== �S�ϐ��m�[�h�̑ΐ��ޓx����v�Z ====//
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
    //==== �S�Ẵ`�F�b�N�m�[�h����אڂ���ϐ��m�[�h��0�𑗂� ====//
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

        //cout << "====�ϐ��m�[�h����====\n";

        //==== �ϐ��m�[�h���� & �ꎞ����r�b�g�̌��� ====//
        unsigned long renewCounter = 0;
        for(int i = 0; i < hc; i++) {
            double llr = logLikeRate[i];    // �ΐ��ޓx��
            //double vSum = sumMessageFromCheckNode(i);    // ���b�Z�[�W�̘a
            //double sendTmp = llr + vSum;

            const vector<int>* checkList = &((*matCol)[i]);
            unsigned long size = checkList->size();

            double vSum = llr;
            for(int j = 0; j < size; j++) {
                vSum += mb.variableBuf[i][j];
            }

            for(int j = 0; j < size; j++) {    // �`�F�b�N�m�[�h�ւ̃��b�Z�[�W���Z�o����
                unsigned long checkId = (*checkList)[j] - 1;
                double writeValue = vSum - mb.getVariableMessage(i, checkId);
                //cout << "x" << i << "�� c" << checkId << " : " << writeValue << "\n";
                mb.variableToCheck(i, checkId, writeValue);
            }

            //cout << sendTmp << " ";

            //---- �e�X�g�p�������� ----//
            //hatCheck[loop][i] = vSum;
            //---- �e�X�g�p�����܂� ----//

            //==== ����r�b�g�̌��� ====//
            if(vSum < 0.0) {
                nHat[i] = 1;
            } else if(vSum > 0.0) {
                nHat[i] = 0;
            } else {
                double rnd = genrand_real2();
                if(rnd < 0.5) nHat[i] = 0;
                else nHat[i] = 1;
            }

            //==== �X�V���i��ł��邩 ====//
            if(renewChecker[i] == vSum) {
                renewCounter++;
            }
            renewChecker[i] = vSum;
        }

        //==== �p���e�B���� ====//
        /*for(int i = 0; i < hc; i++) {
            cout << nHat[i] << " ";
        }
        cout << "\n";*/
        if(parityCheck(&nHat)) {
        //if(isZero(&nHat)) {    // �S�[��������ɂ����ʗp���Ȃ�
            //cout << "��������!!\n";
            break;
        }
        /*if(renewCounter >= hr) {    // �����
            badConvergenceFrag = true;
            break;
        }*/
        // �������s�̓��C�����[�v�̃��[�v�񐔂Ŕ���

        //cout << "====�`�F�b�N�m�[�h����====\n";

        //==== �`�F�b�N�m�[�h���� ====//
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

            for(int j = 0; j < size; j++) {    // �ϐ��m�[�h�ւ̃��b�Z�[�W���Z�o����
                unsigned long variableId = (*variableList)[j] - 1;    // ���M��ϐ��m�[�h��ID

                //cout << "weiteValue = " << writeValue << ", writeValue2 = " << writeValue2 << "\n";
                double tmp2 = mb.getCheckMessage(i, variableId);    // ���M�悩��̃f�[�^
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
                    if((*variableList)[k] - 1 != variableId) {    // ���鑊�������
                        if(absolute(mb.checkBuf[i][k]) < cMin) {
                            cMin = absolute(mb.checkBuf[i][k]);
                        }
                    }
                }

                double writeValue = scale * signValue * cMin;

                //cout << "c" << i << "�� x" << variableId << " : " << writeValue << "\n";
                mb.checkToVariable(i, variableId, writeValue);
            }
        }



        //break;

    }
    if(loop < maxLoop) conv = loop + 1;    // ������
    else conv = loop;

    //==== �����f�[�^�̐��` ====//
    for(int i = 0; i < hc; i++) {
        retObj.decodeWord[i] = nHat[i];
    }
    retObj.convergenceNumber = conv;
    if(conv >= maxLoop) {
        retObj.failFrag = true;
        retObj.errorType = OVER_MAX_CONVERGENCE;    // �������s
    } else {
        if(badConvergenceFrag) {
            retObj.failFrag = true;
            retObj.errorType = BAD_CONVERGENCE;    // �����
        } else {
            retObj.failFrag = false;
            retObj.errorType = DECODE_SUCCESS;    // ��������(���������ł͂Ȃ�)
        }
    }

    //---- �e�X�g�p�������� ----//
    //hatCheck.exportMatrix("hoge.txt");
    //---- �e�X�g�p�����܂� ----//

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

//==== ���b�Z�[�W�̈ꊇ�v�Z ====//
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

//==== �A�[�N�n�C�p�{���b�N�^���W�F���g ====//
double LDPCDecoder2::atanh(double x) const {
    double retVal = 0.5;
    retVal *= log((1 + x) / (1 - x));
    return retVal;
}

//==== abs�֐� ====//
double LDPCDecoder2::absolute(double x) const {
    if(x >= 0) return x;
    else return -x;
}

//==== sign�֐� ====//
double LDPCDecoder2::sign(double x) const {
    if(x > 0) return 1.0;
    else if(x < 0) return -1.0;
    else return 0.0;
}

//==== gallager��f�֐� ====//
double LDPCDecoder2::gallagerf(double x) const {
    //if(x > 0.00001) {
        return log( (exp(x) + 1.0) / (exp(x) - 1.0) );
    //} else {
    //    return 12.21;
    //}
}

//==== �p���e�B���� ====//
bool LDPCDecoder2::parityCheck(const vector<BinaryFiniteField>* nHat) const {    // ��ʗp�BnHat��������ł����true��Ԃ�
    unsigned long hr = H->row();
    unsigned long hc = H->col();
    const vector<vector<int> >* matRow = H->matRow;
    for(int i = 0; i < hr; i++) {    // ���ꂼ��̍s�ɑ΂��Ď��s
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
bool LDPCDecoder2::isZero(const vector<BinaryFiniteField>* nHat) const {    // ���x�d���p�B����������o�ł��Ȃ���A�S�[��������ɂ����ʗp���Ȃ��̂Œ���
    unsigned long size = H->col();
    for(int i = 0; i < size; i++) {
        if((*nHat)[i] != 0) return false;
    }
    return true;
}
