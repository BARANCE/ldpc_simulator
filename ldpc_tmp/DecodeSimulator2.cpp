#include "DecodeSimulator2.h"
#include <math.h>

//==== �R���X�g���N�^ ====//
DecodeSimulator2::DecodeSimulator2(const SPMatrix* checkMatrix) : H(checkMatrix) {
    randFrag = 0;
}

//==== �������� ====//
DecodeSimulator2::DSOutput DecodeSimulator2::decodingSimulation(DSInput simulationParameter) {
    startClock = clock();    // ���Ԍv���J�n

    //==== �ϐ��錾�E�O���� ====//
    unsigned long hr = H->row();
    unsigned long hc = H->col();
    vector<BinaryFiniteField> codeWord = simulationParameter.cordWord;    // ������
    vector<double> transWord(hc);    // ���M��

    //==== �ʐM�H�̐ݒ� ====//
    EChannel channel = simulationParameter.channel;
    double SNR = simulationParameter.SNR;    // SN��[dB] or ���]�m���E�����m��
    const double codeRate = (double)(hc - hr) / (double)hc;    // ��������
    double variance = 0.5 * (1.0 / pow(10.0, SNR / 10.0)) * (double)hc / (double)(hc - hr);    // SN�䂩�狁�߂��G���̕��U
    double averageOfNoise = 0.0;    // �G���̕W�{����
    double varianceOfNoise = 0.0;    // �G���̕W�{���U

    unsigned long maxLoop = simulationParameter.maxLoop;    // �ő厎�s��
    unsigned long loop = 0;    // ���s��
    double scale = simulationParameter.scale;

    unsigned long convergenceFailNumber = simulationParameter.convergenceFailNumber;

    //==== ���ؗp�ϐ� ====//
    unsigned long maxError = simulationParameter.maxError;    // �ő�u���b�N����
    unsigned long numberOfSuccess = 0;    // ����������������
    unsigned long numberOfBlockError = 0;    // �u���b�N��肪����������
    unsigned long numberOfBitError = 0;    // �r�b�g��肪����������
    unsigned long error1 = 0;    // �������肪����������
    unsigned long error2 = 0;    // �������s�ɂ���肪����������
    unsigned long error3 = 0;    // ������ɂ���肪����������
    unsigned long maxConvergenceNumber = 0;    // �ő������
    unsigned long minConvergenceNumber = convergenceFailNumber;    // �ŏ�������
    unsigned long averageConvergenceNumber = 0;    // ���ώ�����
    //unsigned long averageOfBitError;    // �r�b�g���̕��ό�
    unsigned long minOfBitError = hc;    // �r�b�g���̍ŏ���
    unsigned long maxOfBitError = 0;    // �r�b�g���̍ő��
    vector<unsigned long> distributionOfNumberOfBitError(hc, 0);    // �r�b�g�����̕��z
    unsigned long averageOfBurstLength = 0;    // �o�[�X�g���̕���
    unsigned long errorFirst = hc;
    unsigned long errorEnd = 0;
    DSOutput retObj;    // �ԋp�\����
    int convCheckFrag = simulationParameter.convCheckFrag;    // �����`�F�b�N�����邩�ǂ���
    int convCheckSize = simulationParameter.convCheckSize;
    vector<vector<unsigned long> > convCheckCounter(convCheckSize, vector<unsigned long>(hc, 0));    // �����񐔂��Ƃ̌�����񐔂��L�^
    vector<vector<double> > convCheckProb(convCheckSize, vector<double>(hc, 0.0));    // �����񐔂��Ƃ̌�藦���L�^
    LDPCDecoder2::ECDList researchParameter = {convCheckFrag, convCheckSize};

    //==== ����ȕ��������p�̃p�����[�^ ====//
    double shortenRate = simulationParameter.shortenRate;

    //==== ������̐ݒ� ====//
    EDMethod decodingMethod = simulationParameter.decodingMethod;    // �����@
    LDPCDecoder2 dec(H);    // ������
    //dec.setChannel(LDPCDecoder2::EChannel::C_AWGN);    // AWGN�ʐM�H
    dec.setSNR(SNR);    // SN����Z�b�g
    LDPCDecoder2::ESPInput inputParameter = {vector<double>(hc, 0.0), convergenceFailNumber};    // ������ւ̓��̓f�[�^(���M��, �������s����)
    // �����M��͎G����t�����Đ������邽�߁A�����ł͏������̂ݍs���B
    LDPCDecoder2::ESPOutput outputParameter;    // �����킩��̏o�̓f�[�^
    LDPCDecoder2::ERSPOutput outputBECParameter;

    //==== �V���[�g��������p�ݒ� ====//
    vector<int> shortenIds(hc);    // �r�b�g���V���[�g������ӏ�
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

    //==== �o�C�i��-�o�C�|�[���ϊ�(BPSK���f��) ====//
    for(unsigned long i = 0; i < hc; i++) {
        if((int)codeWord[i] == 0) transWord[i] = 1.0;
        else transWord[i] = -1.0;
    }

    for(loop = 0; loop < maxLoop; loop++) {

        //==== �G���t��(�G���̓��v�������܂�) ====//
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
                double rand = genrand_real1(); //��l������[0,1] (32�r�b�g���x)
                if(SNR > rand) inputParameter.receptionWord[i] = ERASURE;
            }
        }

        //if(channel == EChannel::C_AWGN) cout << "AWGN\n";
        //else if(channel == EChannel::C_BEC) cout << "BEC\n";

        //==== ���� ====//
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

        //==== ���� ====//
        vector<BinaryFiniteField>* decodeWord = &outputParameter.decodeWord;
        unsigned long convn = outputParameter.convergenceNumber;
        if(channel == C_AWGN) {
        if(checkDecodeWord(&codeWord, decodeWord)) {    // ��������
            numberOfSuccess++;
            if(maxConvergenceNumber < convn) maxConvergenceNumber = convn;    // �ő�����񐔂��X�V
            if(minConvergenceNumber > convn) minConvergenceNumber = convn;    // �ŏ������񐔂��X�V
            averageConvergenceNumber += convn;    // ���ώ����񐔂��X�V

            if(loop % 10000 == 0) {
                //==== �o�ߎ��� ====//
                nowClock = clock();
                unsigned long span = nowClock - startClock;
                unsigned int tmpExecTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
                unsigned int tmpExecTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
                unsigned int tmpExecTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;

                //==== �\�� ====//
                double progressRate = (double)numberOfBlockError / maxError;
                printProgress(progressRate);
                cout << "�����s��F" << simulationParameter.fileName << "(" << hr << " �~ " << hc << ", R = " << codeRate << ")\n";
                cout << "�o�ߎ��ԁF" << tmpExecTimeHour << "����" << tmpExecTimeMin << "��" << tmpExecTimeSec << "�b\n";
                cout << "SN��F" << SNR << " (���U�F" << variance << ")\n";
                cout << "�G���W�{���ρF" << averageOfNoise / (loop + 1.0) << ", �W�{���U�F" << varianceOfNoise / (loop + 1.0) << "\n";
                cout << "���s�񐔁F" << loop + 1 << " (�ő�F" << maxLoop << ")\n";
                cout << "�u���b�N��藦�F" << numberOfBlockError / (loop + 1.0) << "\n";
                cout << "�r�b�g��藦�F" << numberOfBitError / ((loop + 1.0) * (double)hc) << "\n";
                cout << "�G���[�ʌ�藦�F" << (double)error1 / (loop + 1.0) << " / " << (double)error2 / (loop + 1.0) << " / " << (double)error3 / (loop + 1.0) << "\n";
                cout << "���ό��r�b�g���F" << (double)numberOfBitError / numberOfBlockError << "(�ŏ��F" << minOfBitError << " / �ő�F" << maxOfBitError <<")\n";
                cout << "���ό��o�[�X�g���F" << (double)averageOfBurstLength / numberOfBlockError << "(�擪�F" << errorFirst << " / �ŏI�F" << errorEnd << ")\n";
                cout << "���ώ����񐔁F" << (double)averageConvergenceNumber / numberOfSuccess << "(�ŏ��F" << minConvergenceNumber << " / �ő�F" << maxConvergenceNumber << ")\n";
            }
        } else {    // �������s
            //==== ���������� ====//
            double tmpSuccess = (double)numberOfSuccess / (loop + 1.0);    // �b�蕜��������

            //==== �u���b�N��藦�̑��� ====//
            numberOfBlockError++;
            double tmpBlockError = (double)numberOfBlockError / (loop + 1.0);    // �b��u���b�N��藦
            double progressRate = (double)numberOfBlockError / maxError;    // �i�s��

            //==== �r�b�g��萔�̑��� ====//
            unsigned long bitCount = 0;    // ���̌�
            unsigned long firstBitError = hc;    // �ŏ��̌��̈ʒu
            unsigned long finalBitError = 0;    // �Ō�̌��̈ʒu
            unsigned long burstLength = 0;    // ���̒���
            bool firstFrag = true;    // �ŏ��̌��ł��邩�ǂ���
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

            distributionOfNumberOfBitError[bitCount - 1]++;    // �r�b�g�����̕��z(0�v�f�ڂɌ��1�̌�������)

            double tmpBitError = (double)numberOfBitError / ((loop + 1.0) * (double)hc);    // �b��r�b�g��藦
            double tmpAveNumBitError = (double)numberOfBitError / numberOfBlockError;    // ���r�b�g�̕��ό�
            double tmpBurstLength = (double)averageOfBurstLength / numberOfBlockError;    // �o�[�X�g���̕���

            //==== ���̎�ނ̑��� ====//
            if(outputParameter.failFrag) {
                if(outputParameter.errorType == LDPCDecoder2::OVER_MAX_CONVERGENCE) {    // �������s
                    error2++;
                } else if(outputParameter.errorType == LDPCDecoder2::BAD_CONVERGENCE) {    // �����
                    error3++;
                } else {
                    cout << "Warning : �K�肳��Ă��Ȃ���肪�������Ă��܂��B\n";
                }
            } else {    // �����
                error1++;
            }

            double tmpError1 = (double)error1 / (loop + 1.0);
            double tmpError2 = (double)error2 / (loop + 1.0);
            double tmpError3 = (double)error3 / (loop + 1.0);

            //==== ������ ====//
            double tmpAveConvergence = (double)averageConvergenceNumber / numberOfSuccess;

            //==== ���s���� ====//
            nowClock = clock();
            unsigned long span = nowClock - startClock;
            unsigned int tmpExecTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
            unsigned int tmpExecTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
            unsigned int tmpExecTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;

            //==== �\�� ====//
            printProgress(progressRate);
            cout << "�����s��F" << simulationParameter.fileName << "(" << hr << " �~ " << hc << ", R = " << codeRate << ")\n";
            cout << "�o�ߎ��ԁF" << tmpExecTimeHour << "����" << tmpExecTimeMin << "��" << tmpExecTimeSec << "�b\n";
            cout << "SN��F" << SNR << " (���U�F" << variance << ")\n";
            cout << "�G���W�{���ρF" << averageOfNoise / (loop + 1.0) << ", �W�{���U�F" << varianceOfNoise / (loop + 1.0) << "\n";
            cout << "���s�񐔁F" << loop + 1 << " (�ő�F" << maxLoop << ")\n";
            cout << "�u���b�N��藦�F" << tmpBlockError << "\n";
            cout << "�r�b�g��藦�F" << tmpBitError << "\n";
            cout << "�G���[�ʌ�藦�F" << tmpError1 << " / " << tmpError2 << " / " << tmpError3 << "\n";
            cout << "���ό��r�b�g���F" << tmpAveNumBitError << "(�ŏ��F" << minOfBitError << " / �ő�F" << maxOfBitError <<")\n";
            cout << "���ό��o�[�X�g���F" << tmpBurstLength << "(�擪�F" << errorFirst << " / �ŏI�F" << errorEnd << ")\n";
            cout << "���ώ����񐔁F" << tmpAveConvergence << "(�ŏ��F" << minConvergenceNumber << " / �ő�F" << maxConvergenceNumber << ")\n";

            //==== �I������ ====//
            if(numberOfBlockError >= maxError) break;
        }
        } else if(channel == C_BEC) {
            vector<double> decodeWord = outputBECParameter.decodeWord;    // ������
            int errorCounter = 0;    // �����V���{���̐�
            for(unsigned long i = 0; i < hc; i++) {
                if(decodeWord[i] == ERASURE) errorCounter++;
            }
            vector<vector<double> > convCheckVec = outputBECParameter.convCheckVec;    // ������100��܂ŋL�^
            if(convCheckFrag) {
                for(int i = 0; i < convCheckSize; i++) {
                    for(unsigned long j = 0; j < hc; j++) {
                        if(convCheckVec[i][j] == ERASURE) {
                            convCheckCounter[i][j]++;
                        }
                    }
                }
            }


            if(errorCounter == 0) {    // ��������
                numberOfSuccess++;
                if(maxConvergenceNumber < convn) maxConvergenceNumber = convn;    // �ő�����񐔂��X�V
                if(minConvergenceNumber > convn) minConvergenceNumber = convn;    // �ŏ������񐔂��X�V
                averageConvergenceNumber += convn;    // ���ώ����񐔂��X�V

                if(loop % 10000 == 0) {
                    //==== �o�ߎ��� ====//
                    nowClock = clock();
                    unsigned long span = nowClock - startClock;
                    unsigned int tmpExecTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
                    unsigned int tmpExecTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
                    unsigned int tmpExecTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;

                    //==== �\�� ====//
                    double progressRate = (double)numberOfBlockError / maxError;
                    printProgress(progressRate);
                    cout << "�����s��F" << simulationParameter.fileName << "(" << hr << " �~ " << hc << ", R = " << codeRate << ")\n";
                    cout << "�o�ߎ��ԁF" << tmpExecTimeHour << "����" << tmpExecTimeMin << "��" << tmpExecTimeSec << "�b\n";
                    cout << "�����m���F" << SNR << " (���U�F" << variance << ")\n";
                    cout << "�G���W�{���ρF" << averageOfNoise / (loop + 1.0) << ", �W�{���U�F" << varianceOfNoise / (loop + 1.0) << "\n";
                    cout << "���s�񐔁F" << loop + 1 << " (�ő�F" << maxLoop << ")\n";
                    cout << "�u���b�N��藦�F" << numberOfBlockError / (loop + 1.0) << "\n";
                    cout << "�r�b�g��藦�F" << numberOfBitError / ((loop + 1.0) * (double)hc) << "\n";
                    cout << "�G���[�ʌ�藦�F" << (double)error1 / (loop + 1.0) << " / " << (double)error2 / (loop + 1.0) << " / " << (double)error3 / (loop + 1.0) << "\n";
                    cout << "���ό��r�b�g���F" << (double)numberOfBitError / numberOfBlockError << "(�ŏ��F" << minOfBitError << " / �ő�F" << maxOfBitError <<")\n";
                    cout << "���ό��o�[�X�g���F" << (double)averageOfBurstLength / numberOfBlockError << "(�擪�F" << errorFirst << " / �ŏI�F" << errorEnd << ")\n";
                    cout << "���ώ����񐔁F" << (double)averageConvergenceNumber / numberOfSuccess << "(�ŏ��F" << minConvergenceNumber << " / �ő�F" << maxConvergenceNumber << ")\n";
                }
            } else {    // �������s
                //==== ���������� ====//
                double tmpSuccess = (double)numberOfSuccess / (loop + 1.0);    // �b�蕜��������

                //==== �u���b�N��藦�̑��� ====//
                numberOfBlockError++;
                double tmpBlockError = (double)numberOfBlockError / (loop + 1.0);    // �b��u���b�N��藦
                double progressRate = (double)numberOfBlockError / maxError;    // �i�s��

                //==== �r�b�g��藦 ====//
                unsigned long bitCount = 0;    // ���̌�
                unsigned long firstBitError = hc;    // �ŏ��̌��̈ʒu
                unsigned long finalBitError = 0;    // �Ō�̌��̈ʒu
                unsigned long burstLength = 0;    // ���̒���
                bool firstFrag = true;    // �ŏ��̌��ł��邩�ǂ���
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

                distributionOfNumberOfBitError[bitCount - 1]++;    // �r�b�g�����̕��z(0�v�f�ڂɌ��1�̌�������)

                double tmpBitError = (double)numberOfBitError / ((loop + 1.0) * (double)hc);    // �b��r�b�g��藦
                double tmpAveNumBitError = (double)numberOfBitError / numberOfBlockError;    // ���r�b�g�̕��ό�
                double tmpBurstLength = (double)averageOfBurstLength / numberOfBlockError;    // �o�[�X�g���̕���

                //==== ���̎�ނ̑��� ====//
                if(outputParameter.failFrag) {
                    if(outputParameter.errorType == LDPCDecoder2::OVER_MAX_CONVERGENCE) {    // �������s
                        error2++;
                    } else if(outputParameter.errorType == LDPCDecoder2::BAD_CONVERGENCE) {    // �����
                        error3++;
                    } else {
                        cout << "Warning : �K�肳��Ă��Ȃ���肪�������Ă��܂��B\n";
                    }
                } else {    // �����
                    error1++;
                }

                double tmpError1 = (double)error1 / (loop + 1.0);
                double tmpError2 = (double)error2 / (loop + 1.0);
                double tmpError3 = (double)error3 / (loop + 1.0);

                //==== ������ ====//
                double tmpAveConvergence = (double)averageConvergenceNumber / numberOfSuccess;

                //==== ���s���� ====//
                nowClock = clock();
                unsigned long span = nowClock - startClock;
                unsigned int tmpExecTimeHour = span / (CLOCKS_PER_SEC * 60 * 60);
                unsigned int tmpExecTimeMin = (span % (CLOCKS_PER_SEC * 60 * 60)) / (CLOCKS_PER_SEC * 60);
                unsigned int tmpExecTimeSec = (span % (CLOCKS_PER_SEC * 60)) / CLOCKS_PER_SEC;

                //==== �\�� ====//
                printProgress(progressRate);
                cout << "�����s��F" << simulationParameter.fileName << "(" << hr << " �~ " << hc << ", R = " << codeRate << ")\n";
                cout << "�o�ߎ��ԁF" << tmpExecTimeHour << "����" << tmpExecTimeMin << "��" << tmpExecTimeSec << "�b\n";
                cout << "�����m���F" << SNR << " (���U�F" << variance << ")\n";
                cout << "�G���W�{���ρF" << averageOfNoise / (loop + 1.0) << ", �W�{���U�F" << varianceOfNoise / (loop + 1.0) << "\n";
                cout << "���s�񐔁F" << loop + 1 << " (�ő�F" << maxLoop << ")\n";
                cout << "�u���b�N��藦�F" << tmpBlockError << "\n";
                cout << "�r�b�g��藦�F" << tmpBitError << "\n";
                cout << "�G���[�ʌ�藦�F" << tmpError1 << " / " << tmpError2 << " / " << tmpError3 << "\n";
                cout << "���ό��r�b�g���F" << tmpAveNumBitError << "(�ŏ��F" << minOfBitError << " / �ő�F" << maxOfBitError <<")\n";
                cout << "���ό��o�[�X�g���F" << tmpBurstLength << "(�擪�F" << errorFirst << " / �ŏI�F" << errorEnd << ")\n";
                cout << "���ώ����񐔁F" << tmpAveConvergence << "(�ŏ��F" << minConvergenceNumber << " / �ő�F" << maxConvergenceNumber << ")\n";

                //==== �I������ ====//
                if(numberOfBlockError >= maxError) break;
            }
        }

        //break;
    }

    // �����󋵂̃`�F�b�N
    if(convCheckFrag) {
        for(int i = 0; i < convCheckSize; i++) {
            for(unsigned long j = 0; j < hc; j++) {
                convCheckProb[i][j] = convCheckCounter[i][j] / (double)(loop + 1.0);
            }
        }
    }

    //==== �������ʂ̐��` ====//
    retObj.loop = loop + 1;    // ���ۂ̎��s��
    retObj.blockErrorRate = (double)numberOfBlockError / (loop + 1.0);    // �u���b�N��藦
    retObj.bitErrorRate = (double)numberOfBitError / ((loop + 1.0) * (double)hc);    // �r�b�g��藦
    retObj.decodeSuccessRate = (double)numberOfSuccess / (loop + 1.0);    // ����������
    retObj.error1Rate = (double)error1 / (loop + 1.0);    // �������
    retObj.error2Rate = (double)error2 / (loop + 1.0);    // �������s��
    retObj.error3Rate = (double)error3 / (loop + 1.0);    // �������
    retObj.averageConvergenceNumber = (double)averageConvergenceNumber / numberOfSuccess;    // ���ώ�����
    retObj.minConvergenceNumber = minConvergenceNumber;    // �ŏ�������(��������x���������Ă��Ȃ��ꍇ�͕������ɓ�����)
    retObj.maxConvergenceNumber = maxConvergenceNumber;    // �ő������(��������x���������Ă��Ȃ��ꍇ��0)
    retObj.distributionOfNumberOfBitError = distributionOfNumberOfBitError;    // ���r�b�g���̕��z
    retObj.averageOfNumberOfBitError = (double)numberOfBitError / numberOfBlockError;    // ���ό��r�b�g��
    retObj.averageOfGaussianNoise = averageOfNoise / (loop + 1.0);    // �G���̕���
    retObj.varianceOfGaussianNoise = varianceOfNoise / (loop + 1.0);    // �G���̕��U
    retObj.averageOfBurstLength = (double)averageOfBurstLength / numberOfBlockError;    // ���ό��o�[�X�g��
    retObj.execTime = startClock - clock();    // ���s����

    //==== ���ʂ̏o�� ====//
    stringstream resultss;
    resultss << "result_" << simulationParameter.fileName;
    string resultFileName = resultss.str();    // ���ʂ̃t�@�C����
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

    stringstream resultss2;    // ��蕪�z�p
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
        stringstream resultss3;    // �����`�F�b�N�p
        resultss3 << "conv_SNR" << SNR << "_" << removeExtension(simulationParameter.fileName) << ".txt";
        string convFileName = resultss3.str();
        ofstream ofs3(convFileName.c_str(), fstream::out);    // �������ݗp
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

        stringstream resultss4;    // �����`�F�b�N�p2
        resultss4 << "conv_EP" << SNR << "_" << removeExtension(simulationParameter.fileName) << ".plt";
        string plotName = resultss4.str();
        ofstream ofs4(plotName.c_str(), fstream::out);    // gnuplot�p
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

//==== �G������ ====//
double DecodeSimulator2::getGaussianNoise(double average, double variance) const {
    double retVal;
    double r1 = 1.0 - genrand_real2();    // (0, 1]��l����
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

//==== ������Ƒ��M�����ꂪ���������ǂ��� ====//
bool DecodeSimulator2::checkDecodeWord(const vector<BinaryFiniteField>* codeWord, const vector<BinaryFiniteField>* recepWord) const {
    unsigned long size = codeWord->size();
    for(unsigned long i = 0; i < size; i++) {
        if((*codeWord)[i] != (*recepWord)[i]) return false;
    }
    return true;
}

//==== �i���󋵂̕\�� ====//
void DecodeSimulator2::printProgress(double progressRate) const {
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

//==== str����g���q����菜�� ====//
string DecodeSimulator2::removeExtension(string str) const {
    int pos = str.find_last_of(".");
    if(pos != -1) {
        return str.substr(0, pos);
    } else {
        return str;
    }
}
