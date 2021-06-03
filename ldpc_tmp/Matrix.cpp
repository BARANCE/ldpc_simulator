#include "SPMatrix.h"
#include "SMConstructor.h"
#include <iostream>
#include "Picture.h"
#include "GirthCount.h"
#include "DecodeSimulator2.h"
#include "LDPCEncoder2.h"
#include "HashOperator.h"
#include "DensityEvolution.h"
#include "LoopLDPCParser.h"

using namespace std;

void printProgress(double rate);

// ���͒l
// ldpc.exe [execMode] [[constructionData]or[SNR]] [fileName]
// execMode : c���w�肷��ƍ\�����[�h�Ae���w�肷��Ǝ������[�h
// ---- �\�����[�h ----
// constructionData : �\�����@���L�^�����t�@�C��
// fileName : �쐬���������s��̕ۑ���t�@�C��
// ---- �������[�h ----
// SNR : SN��
// fileName : �����s����L�^�����t�@�C����
// --------------------

int main(int argc, char *argv[]) {
    init_genrand((unsigned int)time(NULL));    // ����������

    /*int tmpHr = 10, tmpHc = 10;
    SPMatrix tmpmat(tmpHr, tmpHc);
    for(int i = 0; i < tmpHr; i++) {
        for(int j = 0; j < tmpHc; j++) {
            tmpmat.setValue(i,j,1);
        }
    }
    cout << GirthCount::loopNum(tmpmat) << "\n";*/


    //FuwaMinSum fms;
    //fms.exec("PEGSCCm400_n4000_wc4_b40.spmat");

    /*SPMatrix hoge1("LRLProtograph2_m6_n8_l3_r6_LHat4.txt");
    SPMatrix hoge2("LRLProtograph2_m9_n14_l3_r6_LHat7.txt");
    Matrix<BinaryFiniteField> hogeM1 = hoge1;
    Matrix<BinaryFiniteField> hogeM2 = hoge2;
    unsigned long splitId = 3;
    SMConstructor::connectLRL(&hogeM1, &hogeM2, splitId);*/

    cout << "----------------------------------------\n";
    cout << "�񌳑Ώ̒ʐM�H�ł�BP�����V�~�����[�^\n";
    cout << "Author : BARANCE 2011\n";
    cout << "----------------------------------------\n";

    //==== �ݒ�̓ǂݍ��� ====//
    string settingFileName = "setting.ini";
    bool settingExist = true;
    if(argc >= 2) {
        settingFileName = argv[1];
        settingExist = true;
    }
    HashOperator ho(settingFileName);

    //==== ���[�h�̑I�� ====//
    int constructionMode = 1;    // 0�Ȃ��藦�������[�h�A1�Ȃ猟���s��\�����[�h�A2�Ȃ疧�x���W�@���[�h�A3�Ȃ�s��ϊ����[�h
    /*if(argc >= 2) {
        string mode = argv[1];
        if(mode == "c") {
            constructionMode = 1;
        } else if(mode == "d") {
            constructionMode = 2;
        } else {
            constructionMode = 0;
        }
    }*/
    if(settingExist && ho.tryKey("mode")) {
        if(ho.getValueByKey<string>("mode") == "construct") {
            constructionMode = 1;
        } else if(ho.getValueByKey<string>("mode") == "experiment") {
            constructionMode = 0;
        } else if(ho.getValueByKey<string>("mode") == "densityEvolution") {
            constructionMode = 2;
        } else if(ho.getValueByKey<string>("mode") == "convert") {
            constructionMode = 3;
        } else {
            cout << "[main.cpp] WARNING: mode�ɂ�construct�܂���experiment���w�肵�ĉ������B\n";
        }
    }

    //vector<int> tmp1;
    //tmp1.size();

    if(constructionMode == 1) {
        cout << "�����s��\�����[�h\n";
        /*
            <�ݒ�t�@�C���̍���>
            ���l���w�肵�Ȃ��ꍇ�͊���l���g�p����܂��B
            [���ʍ���]
            saveFile : �����s���ۑ�����t�@�C�����B�g���q���ȗ������ꍇ�́A�t�@�C���t�H�[�}�b�g�ɏ]���Ċ���̊g���q���t�������B����l��"checkMatrix"�ł���B
            format : �t�@�C���`�����w�肷��B�ȉ���3��ނ��w��\�B���e�͑S�ăe�L�X�g�t�@�C���ł���B����l��"alist"�B
                - matrix : �����s����X�y�[�X��؂��0��1�̍s��`��(�����k)�ŕ\������B�擪�s�ɍs���Ɨ񐔂��L�q����B�����s��̑傫���ɂ���ẮA�t�@�C���T�C�Y�����ɑ傫���Ȃ邽�ߒ��ӂ��邱�ƁB�t���g���q��txt�B
                - alist : Mackay�̃f�[�^�x�[�X�Ɏg���Ă��錟���s��̕\���`���ł���B�������E�c�����Ƃ��ɑ������e�Ղł��邱�Ƃ���A�{�v���O�����̕�����ŗp���Ă���B�t���g���q��alist�B
                - spmat : �{�������̌����s��t�H�[�}�b�g�Ƃ��ė��p����Ă���`���ŁA�t�@�C���T�C�Y���������󂯓n���ɕ֗��B�������̎��Y�ƂȂ��Ă��鑽���̎����v���O�����ł�spmat�`���𗘗p����B�t���g���q��spmat�B
            constructMethod : �����s��̍\���@���w�肷��B��ɕ������̗��̂��g�p���Ă���B�ȉ��Ɏ�����ނ��w��\�ł���B����l��"GRC"�ł���B
                - GRC : gallager�̃����_���\���@�B���ȏ���LDPC�����̍\���@�Ƃ��ďЉ��Ă�����@�ł���B�����ɗ����A���S���Y�����܂ނ��߁A�o�͂����x�قȂ�B
                - PEG : Progressive Edge Growth�B���ݒm���Ă���LDPC�����̍\���@�ōł��������x���̌������\�͂����������\���ł���B�������Ԃ����ɒ����B
                - LDPCC : LDPC��ݍ��ݕ����B�v���g�O���t�̃R�s�[�ƌJ��Ԃ��z�u�ɂ��A�я�̌����s��𐶐�����B�v���g�O���t�̎w�肪�K�v�B
                - LRLProtograph : LRL��Ԍ��������ɗp����v���g�O���t���쐬����B�я��1��"��"�ɔz�u���ꂽ�s�񂪐��������B�����s��Ƃ��ė��p���邱�Ƃ͔񐄏��B
                - ProtographCode : �v���g�O���t�����B�v���g�O���t�̃R�s�[�ƁA���ʒu�̃m�[�h�Ԃł̕ӂ̌q���ς��ɂ���ċ���Ȍ����s��������ɐ����ł���B�v���g�O���t�̎w�肪�K�v�B�����ɗ����A���S���Y�����܂ނ��߁A�o�͂����x�قȂ�B
                - PEGSCC : PEG��Ԍ��������BPEG�̎�@�ɂ���āA�я��1�̑��݂��镔�������Ԍ����s��𐶐�����B�������x��PEG����⑁�����x�B
            [�\���@�ɂ���ė��p���鍀��]
            �������ł́A�����s��̍\���p�����[�^���L�q����B
            rowSize : �s���A�����L�����A�܂��̓`�F�b�N�m�[�h���B����l��600�B
            colSize : �񐔁A�ϐ��̐��A�܂��͕ϐ��m�[�h���B����l��4000�B
            wc : ��d�݁A�܂��͕ϐ��m�[�h�̎����B����l��4�B
            wr : �s�d�݁A�܂��̓`�F�b�N�m�[�h�̎����B�ꕔ�̍\���@�ł̂ݎg�p����B����l��60�B
            bandWidth : LDPC��ݍ��ݕ����A��Ԍ��������̑т̕��B�����Ŏw�肷��l�́A�c�����̕��ł���_�ɒ��ӁB����l��200�B
            copyNum : LDPC��ݍ��ݕ����A�v���g�O���t�����̃R�s�[�񐔂�\���l�B����l��4�B
        */

        //==== �����s��̍\�� ====//
        string constructMethod = "GRC";    // �\���@
        int shr = 2000;    // �s��
        int shc = 4000;    // ��
        vector<int> wc(shc);
        int weight = 3;    // ��d��
        for(int i = 0; i < shc; i++) {
            wc[i] = weight;
        }
        int rweight = 6;    // �s�d��(LRLProtograph�AMAC�Ȃǈꕔ�ł����g�p���Ȃ�)
        int bandWidth = 60;    // �т̕�(PEG��Ԍ��������ł̂ݗ��p)
        int copyNum = 50;    // �R�s�[�񐔁B�v���g�O���t�����ALDPC��ݍ��ݕ����Ŏg�p
        int blockSize = 40;    // �u���b�N�T�C�Y�BMAC�ŗ��p
        unsigned long connectNum = 40;    // �������B��Ԍ��������ŗ��p
        int connectSize = 1;    // �������̑傫��(���󌋍������ŗ��p)
        int CCMode = 12;    // �s�񒆉����ڑ����[�h�BCCLRLP�ŗ��p[1:��,2:�E,4:��,8:���̉��Z�l�ŕ\��]
        unsigned long CCCenterConnectWidth = 40;    // �������ɐڑ�����m�[�h��(�}�X�^�[)
        unsigned long CCBlockWidth = 10;    // �������ɐڑ�����m�[�h��(��)
        double rate = 0.0;    // �ڑ����̑傫����L�ɑ΂��銄��(���󌋍������ŗ��p)
        string loadMatrixName = "";    // �ǂݍ��ތ����s��̃t�@�C����

        int pictureSwit = 0;    // �摜�t�@�C���𐶐����邩�ǂ���(0�ŕW���ݒ�A1�ŃI���A2�ŃI�t)
        int pictureMaxCol = 8000;    // pictureSwit��0(�W���ݒ�)�̏ꍇ�ɁA�摜���o�͂���s��̗񐔂̍ő�l

        if(settingExist) {
			
            if(ho.tryKey("constructMethod")) constructMethod = ho.getValueByKey<string>("constructMethod");
            if(ho.tryKey("rowSize")) shr = ho.getValueByKey<unsigned long>("rowSize");
            if(ho.tryKey("colSize")) shc = ho.getValueByKey<unsigned long>("colSize");
            if(ho.tryKey("wc")) {
                int wcTmp = ho.getValueByKey<int>("wc");
				vector<int> wcTmp2(shc);	// bug
				wc = wcTmp2;
                weight = wcTmp;
                for(int i = 0; i < shc; i++) {
                    wc[i] = weight;
                }
            }
            if(ho.tryKey("wr")) rweight = ho.getValueByKey<int>("wr");
            if(ho.tryKey("bandWidth")) bandWidth = ho.getValueByKey<int>("bandWidth");
            if(ho.tryKey("copyNum")) copyNum = ho.getValueByKey<int>("copyNum");
            if(ho.tryKey("blockSize")) blockSize = ho.getValueByKey<int>("blockSize");
            if(ho.tryKey("connectNum")) connectNum = ho.getValueByKey<unsigned long>("connectNum");
            if(ho.tryKey("connectSize")) connectSize = ho.getValueByKey<int>("connectSize");
            if(ho.tryKey("rate")) rate = ho.getValueByKey<double>("rate");
            if(ho.tryKey("loadMatrixName")) loadMatrixName = ho.getValueByKey<string>("loadMatrixName");
            if(ho.tryKey("picture")) {    // �摜�̏o�͐ݒ�
                string picStr = ho.getValueByKey<string>("picture");
                if(picStr == "default") {pictureSwit = 0;}
                else if(picStr == "on") {pictureSwit = 1;}
                else if(picStr == "off") {pictureSwit = 2;}
                else {
                    cout << "[Matrix.cpp][warning]�摜�t�@�C���̏o�͎w�肪�s���ł��B\n";
                }
            }
            if(ho.tryKey("pictureMaxCol")) pictureMaxCol = ho.getValueByKey<int>("pictureMaxCol");    // �摜�̏o�͍ő啝
        }

        cout << "�\���@��" << constructMethod << "�ł��B\n";

        //==== �v���g�O���t ====//
        Matrix<BinaryFiniteField> P(20, 40);
        if(settingExist) {
            if(ho.tryKey("protograph")) {
                SPMatrix Ptmp(ho.getValueByKey<string>("protograph"));
                P = Ptmp;
            }
        }
        //P = SMConstructor::ReadAlist("CCLRLP_m42_n80_l3_r6_LHat40_mode1_CCW3_BW1.txt");
        //Matrix<BinaryFiniteField> P = SMConstructor::ReadAlist("PEG2_m200_n2000_wc6.txt");
        //Matrix<BinaryFiniteField> P = SMConstructor::PEGSCC2(shr, shc, wc, bandWidth);
        //Matrix<BinaryFiniteField> P = SMConstructor::LRL(6, 60, 120);

        //==== �\�� ====//
        //Matrix<BinaryFiniteField> H;
        SPMatrix H;
        if(constructMethod == "PEG") H = SMConstructor::PEG2(shr, shc, wc);
        else if(constructMethod == "PEGSCC") H = SMConstructor::PEGSCC2(shr, shc, wc, bandWidth);
        else if(constructMethod == "LDPCC") H = SMConstructor::SpatiallyCoupledCode(&P, copyNum);
        //else if(constructMethod == "LRLProtograph") H = SMConstructor::LRL(weight, rweight, shc);
        else if(constructMethod == "LRLProtograph") H = SMConstructor::LRLProtograph(weight, rweight, connectNum);
        else if(constructMethod == "CCLRLP") H = SMConstructor::CCLRLP(weight, rweight, connectNum, CCMode, CCCenterConnectWidth, CCBlockWidth);
        else if(constructMethod == "ProtographCode") H = SMConstructor::ProtographCode(&P, copyNum);
        else if(constructMethod == "GRC") H = SMConstructor::gallagerRandomConstruction(shc, rweight, weight);
        else if(constructMethod == "MAC") H = SMConstructor::MAC(weight, rweight, blockSize);
        else if(constructMethod == "LoadFile") H = SMConstructor::loadMatrixFile(loadMatrixName);
        else if(constructMethod == "ConnectLRLProtograph") H = SMConstructor::ConnectIndependentLRLProtograph(weight, rweight, connectNum, copyNum, connectSize);
        else if(constructMethod == "ConnectLRLProtographWithRate") H = SMConstructor::ConnectIndependentLRLProtographWithProportionalConnection (weight, rweight, connectNum, copyNum, rate);
        else if(constructMethod == "LoopLDPC") {
            LoopLDPCParser lp(ho);
            vector<LoopLDPC::LBandData> bd = lp.translate();
            //vector<LoopLDPC::LBandData> bd(3);
            H = SMConstructor::loopLDPC(bd);
        } else {
            cout << "[main.cpp] ERROR:�\���@���I������Ă��܂���B\n";
            return 0;
        }

        //==== �u���b�N��Ԍ��������p�ݒ� ====//
        /*int blockNum = 10;    // �u���b�N�̐�
        int conSize = 10;    // ������
        vector<Matrix<BinaryFiniteField>> blocks(blockNum);
        int hrEdge = 100, hrCenter = 36;    // �������ƒ[���̍s��
        vector<int> wcEdge(shc), wcCenter(shc);    // �������ƒ[���̗�d��
        for(int i = 0; i < shc; i++) {
            wcEdge[i] = 4;
            wcCenter[i] = 2;
        }
        Matrix<BinaryFiniteField> blockEdge = SMConstructor::ReadAlist("BlockSCC_m100_n400.txt");
        Matrix<BinaryFiniteField> blockCenter = SMConstructor::ReadAlist("BlockSCC_m36_n400.txt");
        for(int i = 0; i < blockNum; i++) {
            if(i == 0 || i == blockNum - 1) {
                blocks[i] = blockEdge;
            } else {
                blocks[i] = blockCenter;
            }
        }
        Matrix<BinaryFiniteField> H = SMConstructor::BlockSCC(blocks, conSize);*/
        SPMatrix HSP = H;

        //==== �����s��̕ۑ� ====//
        unsigned long hr = H.row();
        unsigned long hc = H.col();
        double codeRate = (double)(hc - hr) / (double)hc;

        string fileName = constructMethod;    // �ۑ���t�@�C��
        stringstream ss1;
        ss1 << fileName << "_m" << hr << "_n" << hc;    // �m����
        //ss1 << "_bn" << blockNum << "_cs" << conSize;    // �\���@�ɂ��⏕���
        if(constructMethod == "PEG") ss1 << "_wc" << weight;    // PEG�p
        else if(constructMethod == "PEGSCC") ss1 << "_wc" << weight << "_b" << bandWidth;    // PEGSCC�p
        else if(constructMethod == "LDPCC") ss1 << "_t" << copyNum;    // LDPC��ݍ��ݕ����p
        //else if(constructMethod == "LRLProtograph") ss1 << "_l" << weight << "_r" << rweight;    // LRL�p
        else if(constructMethod == "LRLProtograph") ss1 << "_l" << weight << "_r" << rweight << "_LHat" << connectNum;    // LRL�p
        else if(constructMethod == "CCLRLP") ss1 << "_l" << weight << "_r" << rweight << "_LHat" << connectNum << "_mode" << CCMode << "_CCW" << CCCenterConnectWidth << "_BW" << CCBlockWidth;    // LRL�p
        else if(constructMethod == "ProtographCode") ss1 << "_t" << copyNum;
        else if(constructMethod == "GRC") ss1 << "_wc" << weight << "_wr" << rweight;
        else if(constructMethod == "MAC") ss1 << "_j" << weight << "_k" << rweight << "_p" << blockSize;
        else if(constructMethod == "LoadFile") ss1 << "_name(" << loadMatrixName << ")";
        else if(constructMethod == "ConnectLRLProtograph") ss1 << "_l" << weight << "_r" << rweight << "_LHat" << connectNum << "_copyNum" << copyNum << "_conSize" << connectSize;    // �����Ԍ��������p
        else if(constructMethod == "ConnectLRLProtographWithRate") ss1 << "_l" << weight << "_r" << rweight << "_LHat" << connectNum << "_copyNum" << copyNum << "_conRate" << rate;    // �����Ԍ��������p
        else if(constructMethod == "LoopLDPC") ss1 << "_z" << copyNum;

        ss1 << ".txt";
        fileName = ss1.str();
        /*if (argc >= 4) {
            stringstream ss;
            ss << argv[3];
            fileName = ss.str();    // �R�}���h���C������t�@�C�������w�肷��ꍇ
        }*/
        if(settingExist && ho.tryKey("saveName")) {
            fileName = ho.getValueByKey<string>("saveFile");
        }
        HSP.writeAlist(fileName);
        //HSP.writeSpmat(fileName);
        if(pictureSwit == 1 || (pictureSwit == 0 && hc <= pictureMaxCol)) {    // �摜�̏o��
            Picture pic = Matrix<BinaryFiniteField>(H);
            stringstream ss2;
            ss2 << fileName << ".bmp";
            pic.writeBMP(ss2.str());
        }
        cout << "�����s��F" << fileName << "(" << hr << " �~ " << hc << ", R = " << codeRate << ")\n";
        if(hc < 8000) cout << "girthCount : " << GirthCount::loopNum(H) << "\n";
    } else if (constructionMode == 0) {
        //==== �����s��̓ǂݍ��� ====//
        string fileName = "PEGSCC_m800_n4000_wc6_b200.txt";    // �����s��t�@�C��
        /*if (argc >= 4) {
            stringstream ss;
            ss << argv[3];
            fileName = ss.str();    // �R�}���h���C������t�@�C�������w�肷��ꍇ
        }*/
        if(settingExist) if(ho.tryKey("fileName")) fileName = ho.getValueByKey<string>("fileName");
        SPMatrix H = SMConstructor::ReadAlist(fileName);
        //Matrix<BinaryFiniteField> H = HSP;

        //==== 4�T�C�N�����[�v�̐���\�� ====//
        if(H.col() <= 8000) {
            cout << "girthCount : " << GirthCount::loopNum(H) << "\n";
        }

        //==== ������ ====//
        bool encodeFrag = false;    // �������̗L��������
        Matrix<BinaryFiniteField> C(1, H.col());
        if(encodeFrag) {
            cout << "���������s���܂��B\n";
            vector<BinaryFiniteField> message = LDPCEncoder2::getRandomMessage(H.col() - H.row());
            LDPCEncoder2 enc(&H);
            LDPCEncoder2::LENInput lei = {message};
            LDPCEncoder2::LENOutput leo = enc.encode(lei);
            H = leo.checkMatrix;
            Matrix<BinaryFiniteField> C = leo.sendingWord;
        }

        //==== �����p�����[�^ ====//
        DecodeSimulator2::EDMethod decoding = DecodeSimulator2::C_MIN_SUM;    // �����@
        DecodeSimulator2::EChannel channel = DecodeSimulator2::C_AWGN;    // �ʐM�H
        double SNR = 3.0;    // SN��A�܂��͔��]�m���E�����m��
        unsigned long maxLoop = 10000000;    // �ő厎�s��
        unsigned long errorNum = 100;    // �u���b�N���̏��
        unsigned long maxConvergence = 1000;    // �������s����
        int convCheckFrag = 1;    // �����`�F�b�N�����邩�ǂ���
        int convCheckSize = 100;    // �����`�F�b�N�̃o�b�t�@�T�C�Y
        //if(argc >= 3) SNR = atof(argv[2]);
        double scale = 0.8;    // min-sum�̃X�P�[�����O�萔
        double shortenRate = 0.01;

        if(settingExist) {
            if(ho.tryKey("SNR")) SNR = ho.getValueByKey<double>("SNR");
            if(ho.tryKey("maxLoop")) maxLoop = ho.getValueByKey<unsigned long>("maxLoop");
            if(ho.tryKey("errorNum")) errorNum = ho.getValueByKey<unsigned long>("errorNum");
            if(ho.tryKey("maxConvergence")) maxConvergence = ho.getValueByKey<unsigned long>("maxConvergence");
            if(ho.tryKey("scale")) scale = ho.getValueByKey<double>("scale");

            if(ho.tryKey("decodingMethod")) {
                string dms = ho.getValueByKey<string>("decodingMethod");
                if(dms == "SumProduct") decoding = DecodeSimulator2::C_SUM_PRODUCT;
                else if(dms == "MinSum") decoding = DecodeSimulator2::C_MIN_SUM;
                else cout << "[Matrix.cpp] WARNING: decodingMethod���������w�肳��Ă��܂���BMin-Sum�������s���܂��B\n";
            }
            if(ho.tryKey("channel")) {
                string cas = ho.getValueByKey<string>("channel");
                if(cas == "AWGN") channel = DecodeSimulator2::C_AWGN;
                else if(cas == "BEC") channel = DecodeSimulator2::C_BEC;
                else if(cas == "BSC") channel = DecodeSimulator2::C_BSC;
                else cout << "[Matrix.cpp] WARNING: �ʐM�H���������w�肳��Ă��܂���BAWGN�ʐM�H���g�p���܂��B\n";
            }
            if(ho.tryKey("convergenceDataCheck")) {
                string cdc = ho.getValueByKey<string>("convergenceDataCheck");
                if(cdc == "true" || cdc == "TRUE") convCheckFrag = 1;
                else if(cdc == "false" || cdc == "FALSE") convCheckFrag = 0;
                else {
                    cout << "[Matrix.cpp][warning]convergenceDataCheck�͐^�U�l�Ŏw�肵�Ă��������B\n";
                }
            }
        }

        //==== �������� ====//
        DecodeSimulator2 ds2(&H);
        DecodeSimulator2::DSInput dsi = {decoding, channel, fileName, vector<BinaryFiniteField>(H.col(), 0), SNR, maxLoop, errorNum, maxConvergence, scale, shortenRate, convCheckFrag, convCheckSize};
        DecodeSimulator2::DSOutput result = ds2.decodingSimulation(dsi);
    } else if(constructionMode == 2) {
        //cout << "ok\n";

        //==== �����s��̓ǂݍ��� ====//
        string fileName = "LoopLDPC_m18_n28_z2.txt";    // �����s��t�@�C��
        /*if (argc >= 4) {
            stringstream ss;
            ss << argv[3];
            fileName = ss.str();    // �R�}���h���C������t�@�C�������w�肷��ꍇ
        }*/
        if(settingExist) if(ho.tryKey("fileName")) fileName = ho.getValueByKey<string>("fileName");
        SPMatrix HSP = SMConstructor::ReadAlist(fileName);
        Matrix<BinaryFiniteField> H = HSP;

        //==== 4�T�C�N�����[�v�̐���\�� ====//
        cout << "girthCount : " << GirthCount::loopNum(H) << "\n";

        //==== ���x���W�@ ====//
        DensityEvolution de(&HSP);
        long double zeroBorderProb = pow(0.1, 30);    // 0�Ƃ݂Ȃ��m��
        unsigned long maxConvergence = 1000000;    // �������s����
        unsigned long binarySearchLoop = 38;    // �񕪒T���̎��s��
        int convergenceDataCheck = 1;    // �����󋵃`�F�b�N�̐}���쐬���邩�ǂ���
        int separateBorderCheck = 0;
        int specialEPCheck = 0;
        string dataFileName = "deDataFile.txt";
        if(settingExist) {
            if(ho.tryKey("zeroBorderProb")) zeroBorderProb = ho.getValueByKey<double>("zeroBorderProb");
            if(ho.tryKey("maxConvergence")) maxConvergence = ho.getValueByKey<unsigned long>("maxConvergence");
            if(ho.tryKey("binarySearchLoop")) binarySearchLoop = ho.getValueByKey<unsigned long>("binarySearchLoop");
            if(ho.tryKey("convergenceDataCheck")) {
                string cdc = ho.getValueByKey<string>("convergenceDataCheck");
                if(cdc == "true" || cdc == "TRUE") convergenceDataCheck = 1;
                else if(cdc == "false" || cdc == "FALSE") convergenceDataCheck = 0;
            }
        }

        //==== tmp�̈� ====//
        //DensityEvolution::DEInput deip = {0.15740592, pow(0.1, 20), 1000000};
        //de.execDEwithStepCheck(deip, dataFileName, dataFileName);
        //cout << "�����܂�\n";
        //==== tmp�̈悱���܂� ====//

        //double epsilon = 0;    // �f�o�b�O�p�B�����Ă�������
        DensityEvolution::DESInput desi = {zeroBorderProb, maxConvergence, binarySearchLoop, fileName};
        long double epsilon = de.execDESimulation(desi);
        unsigned long girth = GirthCount::loopNum(H);
        cout << "girthCount : " << girth << "\n";

        if(convergenceDataCheck) {
            cout << "�����󋵃`�F�b�N�̂��߂̃t�@�C�����쐬���܂��B\n";
            DensityEvolution::DEInput dataMakeDEI = {epsilon , zeroBorderProb, maxConvergence};
            DensityEvolution::DEOutput dataMakeDEO = de.execDE(dataMakeDEI);
            unsigned long convergenceNumber = dataMakeDEO.convergenceNumber;
            cout << "������(!!)�F" << convergenceNumber << "\n";
            DensityEvolution::DEWSOutput dewso = de.execDEwithStepCheck(dataMakeDEI, dataFileName, convergenceNumber);

            de.saveDEData(fileName, epsilon, dewso.convergenceNumber, girth);
        }
        if(separateBorderCheck) {
            cout << "��������܂ł̔����񐔂��`�F�b�N���܂��B\n";
            //epsilon = 0.45;    // �f�o�b�O�p�B�����Ă�������
            DensityEvolution::DEInput dataMakeDEI = {epsilon, zeroBorderProb, maxConvergence};
            vector<int> connectionArea(0);
            int hc = HSP.col();
            int blocks = 3;
            for(int k = 0; k < blocks; k++) {
                connectionArea.push_back((hc / blocks) * (k + 1) - 4);
                connectionArea.push_back((hc / blocks) * (k + 1) - 3);
                connectionArea.push_back((hc / blocks) * (k + 1) - 2);
                connectionArea.push_back((hc / blocks) * (k + 1) - 1);
            }
            vector<int> blockEndPos;
            DensityEvolution::DEOutput dataMakeDEO = de.execSubDE(dataMakeDEI, connectionArea, blockEndPos);
            unsigned long separateNumber = dataMakeDEO.convergenceNumber;
            cout << "������(!!)�F" << separateNumber << "\n";
        }
        if(specialEPCheck) {
            cout << "����̏����m���ɂ���������񐔂��`�F�b�N���܂��B\n";
            double ep = 0.4;
            DensityEvolution::DEInput dataMakeDEI = {ep, zeroBorderProb, maxConvergence};
            DensityEvolution::DEOutput dataMakeDEO = de.execDE(dataMakeDEI);
            unsigned long convergenceNumber = dataMakeDEO.convergenceNumber;
            cout << "������(!!)�F" << convergenceNumber << "\n";
            DensityEvolution::DEWSOutput dewso = de.execDEwithStepCheck(dataMakeDEI, dataFileName, convergenceNumber);

            de.saveDEData(fileName, epsilon, dewso.convergenceNumber, girth);
        }

        //==== DE�̌��ʂ�ۑ� ====//

    } else if(constructionMode == 3) {    // spmat����alist�̕ϊ�
        string fileName = "PEG_m332_n1000_wc6.txt";    // �����s��t�@�C��
        string newFileName = "PEG_m332_n1000_wc6.spmat";    // �ύX��̃t�@�C��
        int oldFileType = 0;    // �t�@�C���^�C�v 0��alist�A1��spmat�A2�Ńv���[��
        int newFileType = 1;
        /*if (argc >= 4) {
            stringstream ss;
            ss << argv[3];
            fileName = ss.str();    // �R�}���h���C������t�@�C�������w�肷��ꍇ
        }*/
        if(settingExist) {
            if(ho.tryKey("oldFileName")) fileName = ho.getValueByKey<string>("oldFileName");
            if(ho.tryKey("newFileName")) newFileName = ho.getValueByKey<string>("newFileName");
            if(ho.tryKey("oldFileType")) {
                string fileTypeTmp = ho.getValueByKey<string>("oldFileType");
                if(fileTypeTmp == "alist") oldFileType = 0;
                else if(fileTypeTmp == "spmat") oldFileType = 1;
                else if(fileTypeTmp == "basic") oldFileType = 2;
            }
            if(ho.tryKey("newFileType")) {
                string fileTypeTmp = ho.getValueByKey<string>("newFileType");
                if(fileTypeTmp == "alist") newFileType = 0;
                else if(fileTypeTmp == "spmat") newFileType = 1;
                else if(fileTypeTmp == "basic") newFileType = 2;
            }
        }
        SPMatrix HSP;

        if(oldFileType == 0) {
            HSP = SMConstructor::ReadAlist(fileName);
        } else if(oldFileType == 1) {
            HSP = SMConstructor::ReadSpmat(fileName);
        } else {
            //Matrix<BinaryFiniteField> H(fileName);
            //HSP = H;
            exit(1);
        }

        stringstream ss;
        ss << fileName;

        if(newFileType == 0) {
            ss << ".txt";
            newFileName = ss.str();
            HSP.writeAlist(newFileName);
        } else if(newFileType == 1) {
            ss << ".spmat";
            newFileName = ss.str();
            HSP.writeSpmat(newFileName);
        } else {
            ss << ".txt";
            newFileName = ss.str();
            Matrix<BinaryFiniteField> H = HSP;
            H.exportMatrix(newFileName);
        }

    }

    cout << "�I��.\n";

    string last;
    cin>>last;

    return 0;
}

// �w�肵���lB���x�N�g��A�����菜���֐�
template <class Type> vector<Type> removeByValue(const vector<Type>& A, const Type& B) {
    int size = A.size();
    vector<Type> retVec(0);
    for(int i = 0; i < size; i++) {
        if(A[i] == B) {
        } else {
            retVec.push_back(A[i]);
        }
    }
    return retVec;
}
// �w�肵���W��A�ƏW��B�ɂ��āAA��B�|A��B���c���֐�(�r���I�_���a)

    /*int vecSizeTemp = 4000;
    vector<int> tempx(vecSizeTemp);
    for(int i = 0; i < vecSizeTemp; i++) tempx[i] = 3;
    Matrix<BinaryFiniteField> H = SMConstructor::PEGSCC(30, vecSizeTemp, tempx, 15);
    //Matrix<BinaryFiniteField> H = SMConstructor::RPEGSCC(400, vecSizeTemp, tempx, 40);
    //Matrix<BinaryFiniteField> H = SMConstructor::MPEGSCC(400, vecSizeTemp, tempx, 40);
    //Matrix<BinaryFiniteField> H = SMConstructor::RandomSCC(400, vecSizeTemp, tempx, 40);
    //cout << H;

    //Matrix<BinaryFiniteField> P = SMConstructor::ProgressiveEdgeGrowth(36, vecSizeTemp, tempx);
    //Matrix<BinaryFiniteField> H = SMConstructor::MAC(6, 20, 199);
    //Matrix<BinaryFiniteField> P = SMConstructor::MAC(3, 36, 37);
    //Matrix<BinaryFiniteField> P = SMConstructor::LRL2(40, 400, 4);
    //Matrix<BinaryFiniteField> H = SMConstructor::gallagerRandomConstruction(4000, 16, 8);
    //Matrix<BinaryFiniteField> P = SMConstructor::VMProtograph(2, 25, 41);
    //Matrix<BinaryFiniteField> H = SMConstructor::ProtographCode(&P, 10);
    //Matrix<BinaryFiniteField> H = SMConstructor::SpatiallyCoupledCode(&P, 10);
    //SPMatrix Htmp = SMConstructor::gallagerRandomConstruction(400, 50, 25);
    //Htmp.writeAlist("PEGm40_n400_wc4.txt");
    //Htmp.removeCycles();
    //Htmp.writeAlist("gal_0.5rm.txt");
    //Matrix<BinaryFiniteField> H = Htmp;
    //H.exportMatrix("gal_0.5rmd.txt");
    //SPMatrix foo = H;
    //foo.writeAlist("PEGSCC2m400_n4000_wc6.txt");
    //string fileName = "ProtographCode_0.9_2.txt";
    string fileName = "RPEGSCCm400_n4000_wc4_b40.txt";
    if (argc >= 3) {
        stringstream ss;
        ss << argv[2];
        fileName = ss.str();
    }
    //SPMatrix HSP = SMConstructor::ReadAlist(fileName);
    SPMatrix HSP = H;
    //Matrix<BinaryFiniteField> H = HSP;
    //Matrix<BinaryFiniteField> H = SMConstructor::ReadSpmat(fileName);
    ((SPMatrix)H).writeAlist("MPEGSCCm400_n4000_wc4_b40.txt");
    //((SPMatrix)H).writeSpmat("PEGSCCm30_n60_wc3_b6.spmat");
    Picture pic = H;
    pic.writeBMP("MPEGSCCm400_n4000_wc4_b40.bmp");*/

    // ������

    //H = H.getRENF_H();    // �W���`���擾
    //GHConverter ghc;
    //Matrix<BinaryFiniteField> G = ghc.ToGeneratorMatrix(H);    // �����s��
    //H.exportMatrix("renf_PEGSCCm400_n4000_wc4_b40.txt");
    //G.exportMatrix("gen_PEGSCCm400_n4000_wc4_b40.txt");


    //Matrix<BinaryFiniteField> C(1, H.col());

    /*stringstream fss;
    if(argc >= 4) {fss << argv[3];}
    else fss << "result_" << fileName;
    cout << fss.str() << "\n";*/

    //SPMatrix HSP = H;
    /*LDPCDecoder2 dc2(&HSP);
    vector<double> tempVec(6, 1.0);
    tempVec[0] = -0.01;
    LDPCDecoder2::ESPInput espi = {tempVec, 1};
    LDPCDecoder2::ESPOutput espo = dc2.sumProductDecode(espi);
    cout << espo.decodeWord;
    cout << "loop = " << espo.convergenceNumber << "\n";*/


    //DecodeSimulator ds(&H);
    //Matrix<double> result = ds.execAWGNSumProduct(C, SNR, maxLoop, errorNum, maxConvergence);
