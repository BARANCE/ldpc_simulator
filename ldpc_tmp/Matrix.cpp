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

// 入力値
// ldpc.exe [execMode] [[constructionData]or[SNR]] [fileName]
// execMode : cを指定すると構成モード、eを指定すると実験モード
// ---- 構成モード ----
// constructionData : 構成方法を記録したファイル
// fileName : 作成した検査行列の保存先ファイル
// ---- 実験モード ----
// SNR : SN比
// fileName : 検査行列を記録したファイル名
// --------------------

int main(int argc, char *argv[]) {
    init_genrand((unsigned int)time(NULL));    // 乱数初期化

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
    cout << "二元対称通信路でのBP復号シミュレータ\n";
    cout << "Author : BARANCE 2011\n";
    cout << "----------------------------------------\n";

    //==== 設定の読み込み ====//
    string settingFileName = "setting.ini";
    bool settingExist = true;
    if(argc >= 2) {
        settingFileName = argv[1];
        settingExist = true;
    }
    HashOperator ho(settingFileName);

    //==== モードの選択 ====//
    int constructionMode = 1;    // 0なら誤り率実験モード、1なら検査行列構成モード、2なら密度発展法モード、3なら行列変換モード
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
            cout << "[main.cpp] WARNING: modeにはconstructまたはexperimentを指定して下さい。\n";
        }
    }

    //vector<int> tmp1;
    //tmp1.size();

    if(constructionMode == 1) {
        cout << "検査行列構成モード\n";
        /*
            <設定ファイルの項目>
            ※値を指定しない場合は既定値が使用されます。
            [共通項目]
            saveFile : 検査行列を保存するファイル名。拡張子を省略した場合は、ファイルフォーマットに従って既定の拡張子が付加される。既定値は"checkMatrix"である。
            format : ファイル形式を指定する。以下の3種類が指定可能。内容は全てテキストファイルである。既定値は"alist"。
                - matrix : 検査行列をスペース区切りの0と1の行列形式(無圧縮)で表現する。先頭行に行数と列数を記述する。検査行列の大きさによっては、ファイルサイズが非常に大きくなるため注意すること。付加拡張子はtxt。
                - alist : Mackayのデータベースに使われている検査行列の表現形式である。横方向・縦方向ともに走査が容易であることから、本プログラムの復号器で用いている。付加拡張子はalist。
                - spmat : 本研究室の検査行列フォーマットとして利用されている形式で、ファイルサイズが小さく受け渡しに便利。研究室の資産となっている多くの実験プログラムではspmat形式を利用する。付加拡張子はspmat。
            constructMethod : 検査行列の構成法を指定する。主に符号名の略称を使用している。以下に示す種類が指定可能である。既定値は"GRC"である。
                - GRC : gallagerのランダム構成法。教科書でLDPC符号の構成法として紹介されている方法である。内部に乱択アルゴリズムを含むため、出力が毎度異なる。
                - PEG : Progressive Edge Growth。現在知られているLDPC符号の構成法で最も高いレベルの誤り訂正能力を持つ符号を構成できる。処理時間が非常に長い。
                - LDPCC : LDPC畳み込み符号。プロトグラフのコピーと繰り返し配置により、帯状の検査行列を生成する。プロトグラフの指定が必要。
                - LRLProtograph : LRL空間結合符号に用いるプロトグラフを作成する。帯状に1が"密"に配置された行列が生成される。検査行列として利用することは非推奨。
                - ProtographCode : プロトグラフ符号。プロトグラフのコピーと、同位置のノード間での辺の繋ぎ変えによって巨大な検査行列を高速に生成できる。プロトグラフの指定が必要。内部に乱択アルゴリズムを含むため、出力が毎度異なる。
                - PEGSCC : PEG空間結合符号。PEGの手法によって、帯状に1の存在する部分が並ぶ検査行列を生成する。生成速度はPEGよりやや早い程度。
            [構成法によって利用する項目]
            ※ここでは、検査行列の構成パラメータを記述する。
            rowSize : 行数、検査記号数、またはチェックノード数。既定値は600。
            colSize : 列数、変数の数、または変数ノード数。既定値は4000。
            wc : 列重み、または変数ノードの次数。既定値は4。
            wr : 行重み、またはチェックノードの次数。一部の構成法でのみ使用する。既定値は60。
            bandWidth : LDPC畳み込み符号、空間結合符号の帯の幅。ここで指定する値は、縦方向の幅である点に注意。既定値は200。
            copyNum : LDPC畳み込み符号、プロトグラフ符号のコピー回数を表す値。既定値は4。
        */

        //==== 検査行列の構成 ====//
        string constructMethod = "GRC";    // 構成法
        int shr = 2000;    // 行数
        int shc = 4000;    // 列数
        vector<int> wc(shc);
        int weight = 3;    // 列重み
        for(int i = 0; i < shc; i++) {
            wc[i] = weight;
        }
        int rweight = 6;    // 行重み(LRLProtograph、MACなど一部でしか使用しない)
        int bandWidth = 60;    // 帯の幅(PEG空間結合符号でのみ利用)
        int copyNum = 50;    // コピー回数。プロトグラフ符号、LDPC畳み込み符号で使用
        int blockSize = 40;    // ブロックサイズ。MACで利用
        unsigned long connectNum = 40;    // 結合数。空間結合符号で利用
        int connectSize = 1;    // 結合部の大きさ(叉状結合符号で利用)
        int CCMode = 12;    // 行列中央部接続モード。CCLRLPで利用[1:左,2:右,4:上,8:下の加算値で表現]
        unsigned long CCCenterConnectWidth = 40;    // 中央部に接続するノード数(マスター)
        unsigned long CCBlockWidth = 10;    // 中央部に接続するノード数(幅)
        double rate = 0.0;    // 接続部の大きさのLに対する割合(叉状結合符号で利用)
        string loadMatrixName = "";    // 読み込む検査行列のファイル名

        int pictureSwit = 0;    // 画像ファイルを生成するかどうか(0で標準設定、1でオン、2でオフ)
        int pictureMaxCol = 8000;    // pictureSwitが0(標準設定)の場合に、画像を出力する行列の列数の最大値

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
            if(ho.tryKey("picture")) {    // 画像の出力設定
                string picStr = ho.getValueByKey<string>("picture");
                if(picStr == "default") {pictureSwit = 0;}
                else if(picStr == "on") {pictureSwit = 1;}
                else if(picStr == "off") {pictureSwit = 2;}
                else {
                    cout << "[Matrix.cpp][warning]画像ファイルの出力指定が不正です。\n";
                }
            }
            if(ho.tryKey("pictureMaxCol")) pictureMaxCol = ho.getValueByKey<int>("pictureMaxCol");    // 画像の出力最大幅
        }

        cout << "構成法は" << constructMethod << "です。\n";

        //==== プロトグラフ ====//
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

        //==== 構成 ====//
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
            cout << "[main.cpp] ERROR:構成法が選択されていません。\n";
            return 0;
        }

        //==== ブロック空間結合符号用設定 ====//
        /*int blockNum = 10;    // ブロックの数
        int conSize = 10;    // 結合幅
        vector<Matrix<BinaryFiniteField>> blocks(blockNum);
        int hrEdge = 100, hrCenter = 36;    // 中央部と端部の行幅
        vector<int> wcEdge(shc), wcCenter(shc);    // 中央部と端部の列重み
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

        //==== 検査行列の保存 ====//
        unsigned long hr = H.row();
        unsigned long hc = H.col();
        double codeRate = (double)(hc - hr) / (double)hc;

        string fileName = constructMethod;    // 保存先ファイル
        stringstream ss1;
        ss1 << fileName << "_m" << hr << "_n" << hc;    // 確定情報
        //ss1 << "_bn" << blockNum << "_cs" << conSize;    // 構成法による補助情報
        if(constructMethod == "PEG") ss1 << "_wc" << weight;    // PEG用
        else if(constructMethod == "PEGSCC") ss1 << "_wc" << weight << "_b" << bandWidth;    // PEGSCC用
        else if(constructMethod == "LDPCC") ss1 << "_t" << copyNum;    // LDPC畳み込み符号用
        //else if(constructMethod == "LRLProtograph") ss1 << "_l" << weight << "_r" << rweight;    // LRL用
        else if(constructMethod == "LRLProtograph") ss1 << "_l" << weight << "_r" << rweight << "_LHat" << connectNum;    // LRL用
        else if(constructMethod == "CCLRLP") ss1 << "_l" << weight << "_r" << rweight << "_LHat" << connectNum << "_mode" << CCMode << "_CCW" << CCCenterConnectWidth << "_BW" << CCBlockWidth;    // LRL用
        else if(constructMethod == "ProtographCode") ss1 << "_t" << copyNum;
        else if(constructMethod == "GRC") ss1 << "_wc" << weight << "_wr" << rweight;
        else if(constructMethod == "MAC") ss1 << "_j" << weight << "_k" << rweight << "_p" << blockSize;
        else if(constructMethod == "LoadFile") ss1 << "_name(" << loadMatrixName << ")";
        else if(constructMethod == "ConnectLRLProtograph") ss1 << "_l" << weight << "_r" << rweight << "_LHat" << connectNum << "_copyNum" << copyNum << "_conSize" << connectSize;    // 叉状空間結合符号用
        else if(constructMethod == "ConnectLRLProtographWithRate") ss1 << "_l" << weight << "_r" << rweight << "_LHat" << connectNum << "_copyNum" << copyNum << "_conRate" << rate;    // 叉状空間結合符号用
        else if(constructMethod == "LoopLDPC") ss1 << "_z" << copyNum;

        ss1 << ".txt";
        fileName = ss1.str();
        /*if (argc >= 4) {
            stringstream ss;
            ss << argv[3];
            fileName = ss.str();    // コマンドラインからファイル名を指定する場合
        }*/
        if(settingExist && ho.tryKey("saveName")) {
            fileName = ho.getValueByKey<string>("saveFile");
        }
        HSP.writeAlist(fileName);
        //HSP.writeSpmat(fileName);
        if(pictureSwit == 1 || (pictureSwit == 0 && hc <= pictureMaxCol)) {    // 画像の出力
            Picture pic = Matrix<BinaryFiniteField>(H);
            stringstream ss2;
            ss2 << fileName << ".bmp";
            pic.writeBMP(ss2.str());
        }
        cout << "検査行列：" << fileName << "(" << hr << " × " << hc << ", R = " << codeRate << ")\n";
        if(hc < 8000) cout << "girthCount : " << GirthCount::loopNum(H) << "\n";
    } else if (constructionMode == 0) {
        //==== 検査行列の読み込み ====//
        string fileName = "PEGSCC_m800_n4000_wc6_b200.txt";    // 検査行列ファイル
        /*if (argc >= 4) {
            stringstream ss;
            ss << argv[3];
            fileName = ss.str();    // コマンドラインからファイル名を指定する場合
        }*/
        if(settingExist) if(ho.tryKey("fileName")) fileName = ho.getValueByKey<string>("fileName");
        SPMatrix H = SMConstructor::ReadAlist(fileName);
        //Matrix<BinaryFiniteField> H = HSP;

        //==== 4サイクルループの数を表示 ====//
        if(H.col() <= 8000) {
            cout << "girthCount : " << GirthCount::loopNum(H) << "\n";
        }

        //==== 符号化 ====//
        bool encodeFrag = false;    // 符号化の有無を決定
        Matrix<BinaryFiniteField> C(1, H.col());
        if(encodeFrag) {
            cout << "符号化を行います。\n";
            vector<BinaryFiniteField> message = LDPCEncoder2::getRandomMessage(H.col() - H.row());
            LDPCEncoder2 enc(&H);
            LDPCEncoder2::LENInput lei = {message};
            LDPCEncoder2::LENOutput leo = enc.encode(lei);
            H = leo.checkMatrix;
            Matrix<BinaryFiniteField> C = leo.sendingWord;
        }

        //==== 復号パラメータ ====//
        DecodeSimulator2::EDMethod decoding = DecodeSimulator2::C_MIN_SUM;    // 復号法
        DecodeSimulator2::EChannel channel = DecodeSimulator2::C_AWGN;    // 通信路
        double SNR = 3.0;    // SN比、または反転確率・消失確率
        unsigned long maxLoop = 10000000;    // 最大試行回数
        unsigned long errorNum = 100;    // ブロック誤りの上限
        unsigned long maxConvergence = 1000;    // 収束失敗条件
        int convCheckFrag = 1;    // 収束チェックをするかどうか
        int convCheckSize = 100;    // 収束チェックのバッファサイズ
        //if(argc >= 3) SNR = atof(argv[2]);
        double scale = 0.8;    // min-sumのスケーリング定数
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
                else cout << "[Matrix.cpp] WARNING: decodingMethodが正しく指定されていません。Min-Sum復号を行います。\n";
            }
            if(ho.tryKey("channel")) {
                string cas = ho.getValueByKey<string>("channel");
                if(cas == "AWGN") channel = DecodeSimulator2::C_AWGN;
                else if(cas == "BEC") channel = DecodeSimulator2::C_BEC;
                else if(cas == "BSC") channel = DecodeSimulator2::C_BSC;
                else cout << "[Matrix.cpp] WARNING: 通信路が正しく指定されていません。AWGN通信路を使用します。\n";
            }
            if(ho.tryKey("convergenceDataCheck")) {
                string cdc = ho.getValueByKey<string>("convergenceDataCheck");
                if(cdc == "true" || cdc == "TRUE") convCheckFrag = 1;
                else if(cdc == "false" || cdc == "FALSE") convCheckFrag = 0;
                else {
                    cout << "[Matrix.cpp][warning]convergenceDataCheckは真偽値で指定してください。\n";
                }
            }
        }

        //==== 復号実験 ====//
        DecodeSimulator2 ds2(&H);
        DecodeSimulator2::DSInput dsi = {decoding, channel, fileName, vector<BinaryFiniteField>(H.col(), 0), SNR, maxLoop, errorNum, maxConvergence, scale, shortenRate, convCheckFrag, convCheckSize};
        DecodeSimulator2::DSOutput result = ds2.decodingSimulation(dsi);
    } else if(constructionMode == 2) {
        //cout << "ok\n";

        //==== 検査行列の読み込み ====//
        string fileName = "LoopLDPC_m18_n28_z2.txt";    // 検査行列ファイル
        /*if (argc >= 4) {
            stringstream ss;
            ss << argv[3];
            fileName = ss.str();    // コマンドラインからファイル名を指定する場合
        }*/
        if(settingExist) if(ho.tryKey("fileName")) fileName = ho.getValueByKey<string>("fileName");
        SPMatrix HSP = SMConstructor::ReadAlist(fileName);
        Matrix<BinaryFiniteField> H = HSP;

        //==== 4サイクルループの数を表示 ====//
        cout << "girthCount : " << GirthCount::loopNum(H) << "\n";

        //==== 密度発展法 ====//
        DensityEvolution de(&HSP);
        long double zeroBorderProb = pow(0.1, 30);    // 0とみなす確率
        unsigned long maxConvergence = 1000000;    // 収束失敗条件
        unsigned long binarySearchLoop = 38;    // 二分探索の実行回数
        int convergenceDataCheck = 1;    // 収束状況チェックの図を作成するかどうか
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

        //==== tmp領域 ====//
        //DensityEvolution::DEInput deip = {0.15740592, pow(0.1, 20), 1000000};
        //de.execDEwithStepCheck(deip, dataFileName, dataFileName);
        //cout << "おしまい\n";
        //==== tmp領域ここまで ====//

        //double epsilon = 0;    // デバッグ用。消してください
        DensityEvolution::DESInput desi = {zeroBorderProb, maxConvergence, binarySearchLoop, fileName};
        long double epsilon = de.execDESimulation(desi);
        unsigned long girth = GirthCount::loopNum(H);
        cout << "girthCount : " << girth << "\n";

        if(convergenceDataCheck) {
            cout << "収束状況チェックのためのファイルを作成します。\n";
            DensityEvolution::DEInput dataMakeDEI = {epsilon , zeroBorderProb, maxConvergence};
            DensityEvolution::DEOutput dataMakeDEO = de.execDE(dataMakeDEI);
            unsigned long convergenceNumber = dataMakeDEO.convergenceNumber;
            cout << "収束回数(!!)：" << convergenceNumber << "\n";
            DensityEvolution::DEWSOutput dewso = de.execDEwithStepCheck(dataMakeDEI, dataFileName, convergenceNumber);

            de.saveDEData(fileName, epsilon, dewso.convergenceNumber, girth);
        }
        if(separateBorderCheck) {
            cout << "分離するまでの反復回数をチェックします。\n";
            //epsilon = 0.45;    // デバッグ用。消してください
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
            cout << "分離回数(!!)：" << separateNumber << "\n";
        }
        if(specialEPCheck) {
            cout << "特定の消失確率における収束回数をチェックします。\n";
            double ep = 0.4;
            DensityEvolution::DEInput dataMakeDEI = {ep, zeroBorderProb, maxConvergence};
            DensityEvolution::DEOutput dataMakeDEO = de.execDE(dataMakeDEI);
            unsigned long convergenceNumber = dataMakeDEO.convergenceNumber;
            cout << "収束回数(!!)：" << convergenceNumber << "\n";
            DensityEvolution::DEWSOutput dewso = de.execDEwithStepCheck(dataMakeDEI, dataFileName, convergenceNumber);

            de.saveDEData(fileName, epsilon, dewso.convergenceNumber, girth);
        }

        //==== DEの結果を保存 ====//

    } else if(constructionMode == 3) {    // spmat←→alistの変換
        string fileName = "PEG_m332_n1000_wc6.txt";    // 検査行列ファイル
        string newFileName = "PEG_m332_n1000_wc6.spmat";    // 変更後のファイル
        int oldFileType = 0;    // ファイルタイプ 0でalist、1でspmat、2でプレーン
        int newFileType = 1;
        /*if (argc >= 4) {
            stringstream ss;
            ss << argv[3];
            fileName = ss.str();    // コマンドラインからファイル名を指定する場合
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

    cout << "終了.\n";

    string last;
    cin>>last;

    return 0;
}

// 指定した値BをベクトルAから取り除く関数
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
// 指定した集合Aと集合Bについて、A∪B－A∩Bを残す関数(排他的論理和)

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

    // 符号化

    //H = H.getRENF_H();    // 標準形を取得
    //GHConverter ghc;
    //Matrix<BinaryFiniteField> G = ghc.ToGeneratorMatrix(H);    // 生成行列
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
