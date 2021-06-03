/*
    二元値を要素として持つ疎行列を構成するクラス
*/

#if !defined(___Class_SMConstructor)
#define ___Class_SMConstructor

#include <vector>
#include "Matrix.h"
#include "SPMatrix.h"
#include "BinaryFiniteField.h"
#include "LoopLDPC.h"
using namespace std;

template <class Type> ostream& operator<<(ostream& ss, const vector<Type>& M) {
    int vs = M.size();
    for(int i = 0; i < vs; i++) {
        ss << M[i];
        if(i < vs - 1) ss << " ";
    }
    ss << "\n";
    return ss;
}

class SMConstructor {
public:
    //==== 構成法を表す構造体 ====//
    /*typedef enum ___CONS_MODE {
        GallagerRandomConstruction, ModifyedArrayCode, VMProtograph, LRLCode, ProtographCode, SpatiallyCoupledCode, ProgressiveEdgeGrowth, PEGSCC, RingedPEGSCC, ModifyedPEGSCC, RandomSCC, blockBasedSCC
    } ConsMode;*/

    //==== 疎行列読み込みメソッド(クラスメソッド) ====//
    static SPMatrix ReadSpmat(string fileName) {
        SPMatrix spm(fileName, SPMatrix::WDYM);
        return spm;
    }
    static SPMatrix ReadAlist(string fileName) {
        SPMatrix spm(fileName);
        return spm;
    }

    //==== 疎行列構成メソッド(基本的にクラスメソッド) ====//
    //==== ギャラガーのランダム構成法 ====//
    // 列置換によって正則LDPC符号の検査行列を生成する。生成する行列はフルランクとは限らない。
    // colSize : 行列の列数
    // wr : 行重み
    // wc : 列重み
    static SPMatrix gallagerRandomConstruction(int colSize, int wr, int wc) {
        srand( (unsigned int)time( NULL ) );	// 種を指定
		if(colSize % wc != 0) {
			cout << "[SMConstructor::gallagerRandomConstruction][warning]colSizeがwcで割りきれません。\n";
		}
		SPMatrix H(colSize * wc / wr, colSize);	// 行列本体
		
		unsigned long hr = H.row();
		unsigned long hc = H.col();
		for(unsigned long i = 0; i < colSize / wr; i++) {
			unsigned long baseX = i;
			unsigned long baseY = wr * i;
			for(unsigned long j = 0 ; j < wr; j++) {
				H.setValue(baseX, baseY + j, 1);
			}
		}
		unsigned long iterationNumber = wc - 1;
		for(unsigned long i = 0; i < iterationNumber; i++) {
			vector<int> shuffler(colSize, 0);
			for(unsigned long j = 0; j < colSize; j++) {
				shuffler[j] = j;
			}
			std::random_shuffle(shuffler.begin(), shuffler.end());
			for(unsigned long j = 0; j < colSize; j++) {
				for(unsigned long k = 0; k < colSize / wr; k++) {
					BinaryFiniteField tmp = H.getValue(k, shuffler[j]);
					H.setValue(((i + 1) * colSize / wr) + k, j, tmp);
				}
			}
		}

		return H;
    }
    static Matrix<BinaryFiniteField> gallagerRandomConstruction(int rowSize, int colSize) {    // 自動生成版
        int tmp = gcd(rowSize, colSize);
        int wc = colSize / tmp;
        int wr = rowSize / tmp;
        return gallagerRandomConstruction(colSize, wr, wc);
    }

    // 最大公約数
    inline static int gcd( int m, int n ) {
        // 引数に０がある場合は０を返す
        if ( ( 0 == m ) || ( 0 == n ) )
            return 0;

        // ユークリッドの方法
        while( m != n ) {
            if ( m > n ) m = m - n;
            else         n = n - m;
        }
        return m;
    }//gcd

    //==== Modified Array Code(MAC) ====//
    // Array符号を発展させた符号の検査行列を生成する。高符号長時の性能が良い。
    // j : 縦方向のブロック数
    // k : 横方向のブロック数
    // p : 行列を構成するブロック(部分行列)サイズ(p×p)
    static Matrix<BinaryFiniteField> MAC(int paramJ, int paramK, int paramP) {
        Matrix<BinaryFiniteField> retMat(paramJ * paramP, paramK * paramP);    // jp×kp行列
        for(int j = 0; j < paramK; j++) {
            _matrixIdentBlockPut(0, j, paramP, &retMat);
        }
        for(int i = 1; i < paramJ; i++) {
            for(int j = 0; j < paramK; j++) {
                if(i < j) _matrixShiftBlockPut(i, j, paramP, j * (j - i), &retMat);
                else if(i == j) _matrixIdentBlockPut(i, j, paramP, &retMat);
            }
        }
        return retMat;
    }

    //==== ファンデルモンド行列に基づくプロトグラフ ====//
    // プロトグラフ符号で利用する小規模な検査行列を生成する。これ自体はArray符号による検査行列に一致する。
    // paramJ : 列重み
    // paramK : 行重み
    // paramP : ブロックサイズ
    static Matrix<BinaryFiniteField> VMProtograph(int paramJ, int paramK, int paramP) {
        Matrix<BinaryFiniteField> retMat(paramJ * paramP, paramK * paramP);
        for(int i = 0; i < paramJ; i++) {
            for(int j = 0; j < paramK; j++) {
                // 基準オフセットの決定
                int c = i * paramP;    // c行目
                int d = j * paramP;    // d列目

                // インタリープ数の決定
                int x = (i * j) % paramP;    // Iの各要素をxだけ右に巡回シフト

                for(int k = 0; k < x; k++) {
                    retMat[c + k + paramP - x][d + k] = 1;
                }
                for(int k = x; k < paramP; k++) {
                    retMat[c + k - x][d + k] = 1;
                }
            }
        }
        return retMat;
    }

    //==== (l, r, L')空間結合符号のプロトグラフ ====//
    // l : 列重み
    // r : 最大行重み
    // n : 符号長
    static Matrix<BinaryFiniteField> LRL(int l, int r, int n) {
        int hr = (int)((((double)n * l) / r) - 1 + l);
        int hc = n;
        double lsr = (double)l / r;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        //cout << "l = " << l << ", r = " << r << ", n = " << n << "\n";
        //cout << "[!]hr = " << hr << ", hc = " << hc << "\n";
        for(int i = 0; i < l; i++) {
            for(int j = 0; j < n; j++) {
                int baseRow = i + (int)(lsr * j);
                if(baseRow + l - i > hr) baseRow = hr - l + i;
                retMat[baseRow][j] = 1;
            }
        }
        return retMat;
    }

    // 列重みと行数の基準値、列数から(l, r, L')空間結合符号のプロトグラフを作成
    // mは必ず行数より小さくなる
    static Matrix<BinaryFiniteField> LRL2(int m, int n, int wc) {
        int r = (int)((double)(n * wc) / (m - wc + 1));
        return LRL(wc, r, n);
    }
    // l,r,LHatのパラメータから作成
    // l : 列重み
    // r : 最大行重み
    // LHat : 結合数
    static Matrix<BinaryFiniteField> LRLProtograph(unsigned long l, unsigned long r, unsigned long LHat) {
        unsigned long connectNum = LHat;
        unsigned long hr = l + (connectNum - 1);    // 行数
        unsigned long hc = r * connectNum / l;    // 列数
        Matrix<BinaryFiniteField> retMat(hr, hc);
        for(unsigned long i = 0; i < connectNum; i++) {    // 結合ブロックごと
            unsigned long xBase = i;
            unsigned long xBaseMax = xBase + l;
            unsigned long yBase = r * i / l;
            unsigned long yBaseMax = yBase + r / l;
            for(unsigned long j = xBase; j < xBaseMax; j++) {
                for(unsigned long k = yBase; k < yBaseMax; k++) {
                    retMat[j][k] = 1;
                }
            }
        }
        return retMat;
    }

    // LRLプロトグラフの端部を中央部に接続する
    // modeValue : 接続させる位置を2進数で決定[1:左、2:右、4:上、8:下](例)13なら左上下に接続
    // centerCoverWidth : 中央につなげる際の結合幅
    // blockWidth : 中央に繋げるブロックの大きさ(大きいほど多数のノードで中央に繋げる)
    static Matrix<BinaryFiniteField> CCLRLP(unsigned long l, unsigned long r, unsigned long LHat, int modeValue, unsigned long centerCoverWidth, unsigned long blockWidth) {
        Matrix<BinaryFiniteField> retMat = LRLProtograph(l, r, LHat);
        // モードの解析
        vector<int> mode(4, 0);
        if(modeValue % 2 > 0) mode[0] = 1;    // 左
        if((modeValue % 4) / 2 > 0) mode[1] = 1;    // 右
        if((modeValue % 8) / 4 > 0) mode[2] = 1;    // 上
        if(modeValue / 8 > 0) mode[3] = 1;    // 下

        unsigned long hr = retMat.row();
        unsigned long hc = retMat.col();
        // 定義域チェック
        if(mode[0] == 1 || mode[1] == 1) {    // 左右接続
            if(centerCoverWidth > hr) centerCoverWidth = hr;
            if(blockWidth > hc) blockWidth = hc;
        }
        if(mode[2] == 1 || mode[3] == 1) {    // 上下接続
            if(centerCoverWidth > hc) centerCoverWidth = hc;
            if(blockWidth > hr) blockWidth = hr;
        }

        // 接続ビットの配置
        if(mode[0] == 1) {    //左側
            unsigned long x1 = hr / 2 - centerCoverWidth / 2;
            unsigned long x2 = x1 + centerCoverWidth;
            unsigned long y1 = 0;
            unsigned long y2 = blockWidth;
            placeBitsOnMatrix(&retMat, x1, x2, y1, y2);
        }
        if(mode[1] == 1) {    // 右側
            unsigned long x1 = hr / 2 - centerCoverWidth / 2;
            unsigned long x2 = x1 + centerCoverWidth;
            unsigned long y1 = hc - blockWidth;
            unsigned long y2 = hc;
            placeBitsOnMatrix(&retMat, x1, x2, y1, y2);
        }
        if(mode[2] == 1) {    // 上側
            unsigned long x1 = 0;
            unsigned long x2 = blockWidth;
            unsigned long y1 = hc / 2 - centerCoverWidth / 2;
            unsigned long y2 = y1 + centerCoverWidth;
            placeBitsOnMatrix(&retMat, x1, x2, y1, y2);
        }
        if(mode[3] == 1) {    // 下側
            unsigned long x1 = hr - blockWidth;
            unsigned long x2 = hr;
            unsigned long y1 = hc / 2 - centerCoverWidth / 2;
            unsigned long y2 = y1 + centerCoverWidth;
            placeBitsOnMatrix(&retMat, x1, x2, y1, y2);
        }

        return retMat;
    }

    static void placeBitsOnMatrix(Matrix<BinaryFiniteField>* mat, unsigned long x1, unsigned long x2, unsigned long y1, unsigned long y2) {
        for(unsigned long i = x1; i < x2; i++) {
            for(unsigned long j = y1; j < y2; j++) {
                (*mat)[i][j] = 1;
            }
        }
    }


    //==== プロトグラフ符号(protograph code) ====//
    // プロトグラフを複数回コピーし、ノード間の置換を行うことによって、検査行列を生成する
    static SPMatrix ProtographCode(Matrix<BinaryFiniteField>* P, unsigned long T) {
        /*unsigned long pr = P->row();
        unsigned long pc = P->col();
        cout << "rowSize : " << pr * T << ", colSize = " << pc * T << "\n";

        Matrix<BinaryFiniteField> retMat(pr * T, pc * T);
        for(unsigned long i = 0; i < pr; i++) {
            for(unsigned long j = 0; j < pc; j++) {
                if((int)(*P)[i][j] == 1) {
                    // 基準オフセットの決定
                    unsigned long c = i * T;
                    unsigned long d = j * T;
                    vector<int> perm = getPerm(T);    // 置換ベクトルを取得
                    for(unsigned long k = 0; k < T; k++) {
                        retMat[c + perm[k]][d + k] = 1;
                    }
                }
            }
        }
        return retMat;*/
        unsigned long pr = P->row();
        unsigned long pc = P->col();
        cout << "rowSize : " << pr * T << ", colSize = " << pc * T << "\n";

        SPMatrix retMat(pr * T, pc * T);
        for(unsigned long i = 0; i < pr; i++) {
            for(unsigned long j = 0; j < pc; j++) {
                if((int)(*P)[i][j] == 1) {
                    // 基準オフセットの決定
                    unsigned long c = i * T;
                    unsigned long d = j * T;
                    vector<int> perm = getPerm(T);    // 置換ベクトルを取得
                    for(unsigned long k = 0; k < T; k++) {
                        retMat.setValue(c + perm[k], d + k, 1);
                    }
                }
            }
        }
        return retMat;
    }

    //==== 空間結合符号(spatially-coupled code) ====//
    // 小さな検査行列を基とし、対角線上に帯状にパリティシンボルが並ぶ検査行列を生成する。
    // P : 基となる小さなパリティ検査行列
    // T : コピー回数
    static Matrix<BinaryFiniteField> SpatiallyCoupledCode(Matrix<BinaryFiniteField>* P, int T) {
        int pr = P->row();
        int pc = P->col();
        Matrix<BinaryFiniteField> retMat(pr * (T + 1), pc * T);
        for(int k = 0; k < T; k++) {
            // 基準オフセットの決定
            int c = k * pr;    // c行目
            int d = k * pc;    // d列目

            for(int i = 0; i < pr; i++) {
                for(int j = 0; j < pc; j++) {
                    if((int)(*P)[i][j] == 1) {
                        if(j <= (pc * i) / pr) {
                            retMat[c + i][d + j] = 1;
                        } else {
                            retMat[c + i + pr][d + j] = 1;
                        }
                    }
                }
            }
        }
        return retMat;
    }

    //==== Progressive Edge Growth(PEG) ====//
    // タナーグラフの内径(girth)を可能な限り大きくするように、辺の接続を選択しながらグラフを構成する手法
    // rowSize : 検査行列の行数
    // colSize : 検査行列の列数
    // valDim : 変数ノードごとの次数を格納するベクトル。大きさはcolSizeと一致すること
    static Matrix<BinaryFiniteField> ProgressiveEdgeGrowth(int rowSize, int colSize, const vector<int>& valDim) {
        int hr = rowSize;
        int hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> checkDim(hr, 0);    // チェックノードの次数を格納するベクトル
        for(int i = 0; i < hc; i++) {
            // 1本目
            if(valDim[i] > 0) {
                int checkSelect = 0;    // 候補となるチェックノード
                for(int k = 1; k < hr; k++) {
                    if(checkDim[checkSelect] > checkDim[k]) {
                        checkSelect = k;    // 最も小さい次数のノードを選択
                    }
                }
                retMat[checkSelect][i] = 1;    // 辺を配置
                checkDim[checkSelect]++;    // チェックノードの次数上昇
            }

            // 2本目以降
            for(int j = 1; j < valDim[i]; j++) {
                vector<int> prevCheckVec(0);    // 一つ前のチェックノードの番号を格納するベクトル
                vector<int> nowCheckVec(0);    // 現在のチェックノードの番号を格納するベクトル
                vector<int> valList(0);    // 着目する変数ノードのリスト
                vector<int> nextValList(0);
                vector<int> valEnd(0);    // チェック済み変数ノードのリスト
                valList.push_back(i);

                //==== ここからループ処理 ====//
                while(1) {
                    int valNode = valList[valList.size() - 1];
                    contains(valEnd, valNode);
                    valList.pop_back();

                    // 接続するチェックノードのリスト
                    for(int k = 0; k < hr; k++) {
                        if((int)retMat[k][valNode] == 1) {
                            //cout << "nowCheckVec.size = " << nowCheckVec.size() << "→";
                            contains(nowCheckVec, k);
                            //cout << nowCheckVec.size() << "\n";
                        }
                    }

                    // チェックノードリストに接続する変数ノードのリスト
                    for(int k = 0; k < nowCheckVec.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[nowCheckVec[k]][l] == 1 && l != valNode) {    // ここが不明
                                if(!contains2(valEnd, l)) {
                                    contains(nextValList, l);
                                }
                            }
                        }
                    }

                    //cout << prevCheckVec.size() << " → " << nowCheckVec.size() << "\n";

                    //cout << "nowCheckVec = " << nowCheckVec;

                    if(nowCheckVec.size() == hr) {
                        //cout << "満タン";
                        nowCheckVec = prevCheckVec;
                        break;
                    }

                    if(valList.empty()) {

                        if(!nextValList.empty()) {
                            for(int k = 0; k < nextValList.size(); k++) {
                                valList.push_back(nextValList[k]);
                            }
                            nextValList.clear();

                            prevCheckVec.clear();
                            for(int k = 0; k < nowCheckVec.size(); k++) {
                                prevCheckVec.push_back(nowCheckVec[k]);
                            }
                        } else {
                            // チェックノードリストを更新
                            if(eqVec(prevCheckVec, nowCheckVec)) {
                                //cout << "収束";
                                break;
                            }
                            //cout << "空";
                            break;
                        }
                    }

                }
                //==== ループ処理ここまで ====//

                // 補集合を取得
                vector<int> comp = complement(nowCheckVec, hr);

                /*cout << "comp = [";
                for(int b = 0; b < comp.size(); b++) {
                    cout << comp[b] << " ";
                }
                cout << "]\n";*/

                int checkSelect = comp[0];
                for(int k = 1; k < comp.size(); k++) {
                    if(checkDim[checkSelect] > checkDim[comp[k]]) {
                        checkSelect = comp[k];    // 最も小さい次数のノードを選択
                    }
                }
                retMat[checkSelect][i] = 1;    // 辺を配置
                checkDim[checkSelect]++;    // チェックノードの次数上昇
            }
        }
        return retMat;
    }

    //==== PEG空間結合符号(PEG-SCC) ====//
    // PEGに辺を接続する際の制約を付加したもの
    // bandWidth : 帯の大きさ
    static Matrix<BinaryFiniteField> PEGSCC(int rowSize, int colSize, const vector<int>& valDim, int bandWidth) {
        int hr = rowSize;
        int hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> checkDim(hr, 0);    // チェックノードの次数を格納するベクトル
        for(int i = 0; i < hc; i++) {
            //cout << "i = " << i << "\n";
            // 基準オフセットの決定
            int c = ((double)(hr - bandWidth + 1) / hc) * i;
            //cout << "c = " << c << "\n";

            // 1本目
            if(valDim[i] > 0) {
                int checkSelect = c;    // 候補となるチェックノード
                for(int k = c; k < c + bandWidth; k++) {
                    if(checkDim[checkSelect] > checkDim[k]) {
                        checkSelect = k;    // 最も小さい次数のノードを選択
                    }
                }
                //cout << "checkSekect : " << checkSelect << "\n";
                retMat[checkSelect][i] = 1;    // 辺を配置
                checkDim[checkSelect]++;    // チェックノードの次数上昇
            }

            // 2本目以降
            for(int j = 1; j < valDim[i]; j++) {
                vector<int> prevCheckVec(0);    // 一つ前のチェックノードの番号を格納するベクトル
                vector<int> nowCheckVec(0);    // 現在のチェックノードの番号を格納するベクトル
                vector<int> valList(0);    // 着目する変数ノードのリスト
                vector<int> nextValList(0);
                vector<int> valEnd(0);    // チェック済み変数ノードのリスト
                valList.push_back(i);

                //==== ここからループ処理 ====//
                while(1) {
                    int valNode = valList[valList.size() - 1];
                    contains(valEnd, valNode);
                    valList.pop_back();

                    // 接続するチェックノードのリスト
                    for(int k = 0; k < hr; k++) {
                        if((int)retMat[k][valNode] == 1) {
                            //cout << "nowCheckVec.size = " << nowCheckVec.size() << "→";
                            contains(nowCheckVec, k);
                            //cout << nowCheckVec.size() << "\n";
                        }
                    }

                    // チェックノードリストに接続する変数ノードのリスト
                    for(int k = 0; k < nowCheckVec.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[nowCheckVec[k]][l] == 1 && l != valNode) {    // ここが不明
                                if(!contains2(valEnd, l)) {
                                    contains(nextValList, l);
                                }
                            }
                        }
                    }

                    //cout << prevCheckVec.size() << " → " << nowCheckVec.size() << "\n";

                    //cout << "nowCheckVec = " << nowCheckVec;

                    if(containMinToMax(nowCheckVec, c, c + bandWidth)) {
                    //if(nowCheckVec.size() == hr) {    // ここを変えるのかな？
                        //cout << "満タン";
                        nowCheckVec = prevCheckVec;
                        break;
                    }

                    if(valList.empty()) {

                        if(!nextValList.empty()) {
                            for(int k = 0; k < nextValList.size(); k++) {
                                valList.push_back(nextValList[k]);
                            }
                            nextValList.clear();

                            prevCheckVec.clear();
                            for(int k = 0; k < nowCheckVec.size(); k++) {
                                prevCheckVec.push_back(nowCheckVec[k]);
                            }
                        } else {
                            // チェックノードリストを更新
                            if(eqVec(prevCheckVec, nowCheckVec)) {
                                //cout << "収束";
                                break;
                            }
                            //cout << "空";
                            break;
                        }
                    }

                }
                //==== ループ処理ここまで ====//

                // 補集合を取得
                //vector<int> comp = complement(nowCheckVec, hr);
                vector<int> comp;
                for(int k = c; k < c + bandWidth; k++) {
                    if(!contains2(nowCheckVec, k)) comp.push_back(k);
                }
                //cout << "comp : " << comp;

                /*cout << "comp = [";
                for(int b = 0; b < comp.size(); b++) {
                    cout << comp[b] << " ";
                }
                cout << "]\n";*/

                int checkSelect = comp[0];
                for(int k = 1; k < comp.size(); k++) {
                    if(checkDim[checkSelect] > checkDim[comp[k]]) {
                        checkSelect = comp[k];    // 最も小さい次数のノードを選択
                    }
                }
                retMat[checkSelect][i] = 1;    // 辺を配置
                checkDim[checkSelect]++;    // チェックノードの次数上昇
            }
        }
        return retMat;
    }

    //==== Ringed-PEG空間結合符号(R-PEGSCC) ====//
    static Matrix<BinaryFiniteField> RPEGSCC(int rowSize, int colSize, const vector<int>& valDim, int bandWidth) {
        int hr = rowSize;
        int hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> checkDim(hr, 0);    // チェックノードの次数を格納するベクトル
        for(int i = 0; i < hc; i++) {
            //cout << "i = " << i << "\n";
            // 基準オフセットの決定
            int c = ((double)hr / hc) * i;
            //cout << "c = " << c << "\n";

            // 1本目
            if(valDim[i] > 0) {
                int checkSelect = c;    // 候補となるチェックノード
                for(int k = c; k < c + bandWidth; k++) {
                    int sel = k;
                    if(sel >= hr) sel -= hr;
                    if(checkDim[checkSelect] > checkDim[sel]) {
                        checkSelect = sel;    // 最も小さい次数のノードを選択
                    }
                }
                //cout << "checkSekect : " << checkSelect << "\n";
                retMat[checkSelect][i] = 1;    // 辺を配置
                checkDim[checkSelect]++;    // チェックノードの次数上昇
            }

            // 2本目以降
            for(int j = 1; j < valDim[i]; j++) {
                vector<int> prevCheckVec(0);    // 一つ前のチェックノードの番号を格納するベクトル
                vector<int> nowCheckVec(0);    // 現在のチェックノードの番号を格納するベクトル
                vector<int> valList(0);    // 着目する変数ノードのリスト
                vector<int> nextValList(0);
                vector<int> valEnd(0);    // チェック済み変数ノードのリスト
                valList.push_back(i);

                //==== ここからループ処理 ====//
                while(1) {
                    int valNode = valList[valList.size() - 1];
                    contains(valEnd, valNode);
                    valList.pop_back();

                    // 接続するチェックノードのリスト
                    for(int k = 0; k < hr; k++) {
                        if((int)retMat[k][valNode] == 1) {
                            //cout << "nowCheckVec.size = " << nowCheckVec.size() << "→";
                            contains(nowCheckVec, k);
                            //cout << nowCheckVec.size() << "\n";
                        }
                    }

                    // チェックノードリストに接続する変数ノードのリスト
                    for(int k = 0; k < nowCheckVec.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[nowCheckVec[k]][l] == 1 && l != valNode) {    // ここが不明
                                if(!contains2(valEnd, l)) {
                                    contains(nextValList, l);
                                }
                            }
                        }
                    }

                    //cout << prevCheckVec.size() << " → " << nowCheckVec.size() << "\n";

                    //cout << "nowCheckVec = " << nowCheckVec;
                    if(c + bandWidth < hr) {
                        if(containMinToMax(nowCheckVec, c, c + bandWidth)) {
                        //if(nowCheckVec.size() == hr) {    // ここを変えるのかな？
                            //cout << "満タン";
                            nowCheckVec = prevCheckVec;
                            break;
                        }
                    } else {
                        if(containMinToMax(nowCheckVec, 0, c + bandWidth - hr) && containMinToMax(nowCheckVec, c, hr)) {
                            nowCheckVec = prevCheckVec;
                            break;
                        }
                    }

                    if(valList.empty()) {

                        if(!nextValList.empty()) {
                            for(int k = 0; k < nextValList.size(); k++) {
                                valList.push_back(nextValList[k]);
                            }
                            nextValList.clear();

                            prevCheckVec.clear();
                            for(int k = 0; k < nowCheckVec.size(); k++) {
                                prevCheckVec.push_back(nowCheckVec[k]);
                            }
                        } else {
                            // チェックノードリストを更新
                            if(eqVec(prevCheckVec, nowCheckVec)) {
                                //cout << "収束";
                                break;
                            }
                            //cout << "空";
                            break;
                        }
                    }

                }
                //==== ループ処理ここまで ====//

                // 補集合を取得
                //vector<int> comp = complement(nowCheckVec, hr);
                vector<int> comp;
                for(int k = c; k < c + bandWidth; k++) {
                    int sel = k;
                    if(sel >= hr) sel -= hr;
                    if(!contains2(nowCheckVec, sel)) comp.push_back(sel);
                }
                //cout << "comp : " << comp;

                /*cout << "comp = [";
                for(int b = 0; b < comp.size(); b++) {
                    cout << comp[b] << " ";
                }
                cout << "]\n";*/

                int checkSelect = comp[0];
                //cout << "val:" << i << " dim:" << j << " c:" << c << " comp:" << comp;
                for(int k = 1; k < comp.size(); k++) {
                    if(checkDim[checkSelect] > checkDim[comp[k]]) {
                        checkSelect = comp[k];    // 最も小さい次数のノードを選択
                    }
                }
                retMat[checkSelect][i] = 1;    // 辺を配置
                checkDim[checkSelect]++;    // チェックノードの次数上昇
            }
        }
        return retMat;
    }

    //==== 辺配置順序を変更したPEG空間結合符号(PEG-SCC) ====//
    // PEGに辺を接続する際の制約を付加したもの
    // bandWidth : 帯の大きさ
    static Matrix<BinaryFiniteField> MPEGSCC(int rowSize, int colSize, const vector<int>& valDim, int bandWidth) {
        int hr = rowSize;
        int hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> checkDim(hr, 0);    // チェックノードの次数を格納するベクトル
        for(int a = 0; a < hc; a++) {
            int i;    // 配置位置を前→後ろ→前→…とする
            if(a % 2 == 0) {
                i = a / 2;
            } else {
                i = hc - (a / 2) - 1;
            }

            //cout << "i = " << i << "\n";
            // 基準オフセットの決定
            int c = ((double)(hr - bandWidth + 1) / hc) * i;
            //cout << "c = " << c << "\n";

            // 1本目
            if(valDim[i] > 0) {
                int checkSelect = c;    // 候補となるチェックノード
                for(int k = c; k < c + bandWidth; k++) {
                    if(checkDim[checkSelect] > checkDim[k]) {
                        checkSelect = k;    // 最も小さい次数のノードを選択
                    }
                }
                //cout << "checkSekect : " << checkSelect << "\n";
                retMat[checkSelect][i] = 1;    // 辺を配置
                checkDim[checkSelect]++;    // チェックノードの次数上昇
            }

            // 2本目以降
            for(int j = 1; j < valDim[i]; j++) {
                vector<int> prevCheckVec(0);    // 一つ前のチェックノードの番号を格納するベクトル
                vector<int> nowCheckVec(0);    // 現在のチェックノードの番号を格納するベクトル
                vector<int> valList(0);    // 着目する変数ノードのリスト
                vector<int> nextValList(0);
                vector<int> valEnd(0);    // チェック済み変数ノードのリスト
                valList.push_back(i);

                //==== ここからループ処理 ====//
                while(1) {
                    int valNode = valList[valList.size() - 1];
                    contains(valEnd, valNode);
                    valList.pop_back();

                    // 接続するチェックノードのリスト
                    for(int k = 0; k < hr; k++) {
                        if((int)retMat[k][valNode] == 1) {
                            //cout << "nowCheckVec.size = " << nowCheckVec.size() << "→";
                            contains(nowCheckVec, k);
                            //cout << nowCheckVec.size() << "\n";
                        }
                    }

                    // チェックノードリストに接続する変数ノードのリスト
                    for(int k = 0; k < nowCheckVec.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[nowCheckVec[k]][l] == 1 && l != valNode) {    // ここが不明
                                if(!contains2(valEnd, l)) {
                                    contains(nextValList, l);
                                }
                            }
                        }
                    }

                    //cout << prevCheckVec.size() << " → " << nowCheckVec.size() << "\n";

                    //cout << "nowCheckVec = " << nowCheckVec;

                    if(containMinToMax(nowCheckVec, c, c + bandWidth)) {
                    //if(nowCheckVec.size() == hr) {    // ここを変えるのかな？
                        //cout << "満タン";
                        nowCheckVec = prevCheckVec;
                        break;
                    }

                    if(valList.empty()) {

                        if(!nextValList.empty()) {
                            for(int k = 0; k < nextValList.size(); k++) {
                                valList.push_back(nextValList[k]);
                            }
                            nextValList.clear();

                            prevCheckVec.clear();
                            for(int k = 0; k < nowCheckVec.size(); k++) {
                                prevCheckVec.push_back(nowCheckVec[k]);
                            }
                        } else {
                            // チェックノードリストを更新
                            if(eqVec(prevCheckVec, nowCheckVec)) {
                                //cout << "収束";
                                break;
                            }
                            //cout << "空";
                            break;
                        }
                    }

                }
                //==== ループ処理ここまで ====//

                // 補集合を取得
                //vector<int> comp = complement(nowCheckVec, hr);
                vector<int> comp;
                for(int k = c; k < c + bandWidth; k++) {
                    if(!contains2(nowCheckVec, k)) comp.push_back(k);
                }
                //cout << "comp : " << comp;

                /*cout << "comp = [";
                for(int b = 0; b < comp.size(); b++) {
                    cout << comp[b] << " ";
                }
                cout << "]\n";*/

                int checkSelect = comp[0];
                for(int k = 1; k < comp.size(); k++) {
                    if(checkDim[checkSelect] > checkDim[comp[k]]) {
                        checkSelect = comp[k];    // 最も小さい次数のノードを選択
                    }
                }
                retMat[checkSelect][i] = 1;    // 辺を配置
                checkDim[checkSelect]++;    // チェックノードの次数上昇
            }
        }
        return retMat;
    }

    //==== ランダム空間結合符号 ====//
    // bandWidth : 帯の大きさ
    static Matrix<BinaryFiniteField> RandomSCC(int rowSize, int colSize, const vector<int>& valDim, int bandWidth) {
        int hr = rowSize;
        int hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);

        for(int i = 0; i < hc; i++) {
            int c = ((double)(hr - bandWidth + 1) / hc) * i;

            vector<int> randVec = getRandomSequence(c, c + bandWidth, valDim[i]);
            //cout << randVec;

            for(int j = 0; j < valDim[i]; j++) {
                retMat[randVec[j]][i] = 1;
            }
        }

        return retMat;
    }

    //==== ブロック配置空間結合符号 ====//
    // p×q行列を斜めに、少しだけチェックノードが重なるように配置して検査行列を構成する
    // 各行列ごとのサイズを可変とする
    static Matrix<BinaryFiniteField> BlockSCC(vector<Matrix<BinaryFiniteField> > blocks, unsigned long conSize) {
        int blockNum = blocks.size();    // 配置するブロックの数
        unsigned long hr = 0;
        unsigned long hc = 0;
        for(int i = 0; i < blockNum; i++) {
            hr += blocks[i].row();
            hc += blocks[i].col();
            if(blocks[i].row() < conSize) {
                cout << "ERROR : ブロック行列の行数が、結合幅より小さいため検査行列を構成できません。\n";
            }
        }
        hr -= conSize * (blockNum - 1); // 結合部分だけ小さくする
        cout << "hr = " << hr << ", hc = " << hc << "\n";
        Matrix<BinaryFiniteField> retMat(hr, hc);
        unsigned long hrSum = 0;
        unsigned long hcSum = 0;
        for(int i = 0; i < blockNum; i++) {
            unsigned long rbase = hrSum - (i * blockNum);    // 基準位置(行)
            unsigned long cbase = hcSum;    // 基準位置(列)

            unsigned long shr = blocks[i].row();
            unsigned long shc = blocks[i].col();
            for(int j = 0; j < shr; j++) {
                for(int k = 0; k < shc; k++) {
                    retMat[rbase + j][cbase + k] = blocks[i][j][k];
                }
            }

            //==== 基準位置を更新 ====//
            hrSum += shr;
            hcSum += shc;
        }
        return retMat;
    }

    //==== Progressive Edge Growth(PEG) 改訂版 ====//
    // タナーグラフの内径(girth)を可能な限り大きくするように、辺の接続を選択しながらグラフを構成する手法
    // rowSize : 検査行列の行数
    // colSize : 検査行列の列数
    // valDim : 変数ノードごとの次数を格納するベクトル。大きさはcolSizeと一致すること
    static Matrix<BinaryFiniteField> PEG2(int rowSize, int colSize, const vector<int>& valDim) {
        int hr = rowSize;
        int long hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> nowCheckDim(hr);    // 構成中の各チェックノードの現在の次数を表すベクトル
        vector<int> distance(hr);    // 任意の変数ノードの辺接続ステップにおいて、各チェックノードへの距離を一時格納するベクトル。

        //==== チェックノードの次数を0に初期化 ====//
        for(int i = 0; i < hr; i++) nowCheckDim[i] = 0;

        for(int i = 0; i < hc; i++) {
            int selectedVid = i;    // 選択された変数ノード。後に変更することを想定

            //==== i番目の変数ノードから各チェックノードへの距離を無限大(-1)に初期化 ====//
            for(int j = 0; j < hr; j++) distance[j] = -1;

            //==== チェックノードの選択可能範囲(PEG空間結合符号などで利用) ====//
            int startPos = 0;
            int endPos = hr;

            //==== 1本目は次数の一番低いチェックノードに接続 ====//
            int minDimId = startPos;
            int minDim = nowCheckDim[startPos];
            for(int j = startPos + 1; j < endPos; j++) {
                if(nowCheckDim[j] < minDim) {
                    minDimId = j;
                    minDim = nowCheckDim[j];
                }
            }
            retMat[minDimId][i] = 1;    // 1本目の辺を設置
            nowCheckDim[minDimId]++;

            //==== 距離リスト作成(探索)のための初期ノード ====//
            int depth = 0;    // 現在の深さ
            int maxDepth = 10;    // 最大深さ
            vector<int> ovList(0);
            ovList.push_back(i);    // 初期ノードは調査変数ノード

            //==== 変数ノードの次数を満たすまで繰り返し ====//
            for(int j = 1; j < valDim[i]; j++) {
                //==== 最も遠いチェックノードを発見するためのループ ====//
                for(depth = 0; depth < maxDepth; depth++) {
                    //==== チェックノードへの距離を更新 ====//
                    vector<int> children(0);
                    for(int l = 0; l < ovList.size(); l++) {
                        for(int k = 0; k < hr; k++) {
                            if((int)retMat[k][ovList[l]] == 1) {    // 接続しているチェックノードに対して
                                if(distance[k] == -1) {    // 距離が未設定(無限大)→更新
                                    children.push_back(k);
                                    distance[k] = depth;    // 距離を更新
                                }else if(distance[k] != -1 && depth < distance[k]) {    // 現在の距離(depth)より大きい距離が設定されている→更新
                                    children.push_back(k);
                                    distance[k] = depth;    // 距離を更新
                                }
                            }
                        }
                    }
                    //==== チェックノードリストchildrenの重複を取り除く ====//
                    vector<int> childrenRev(0);    // 重複を取り除いたchildren
                    for(int k = 0; k < children.size(); k++) {
                        bool frag = false;    // 重複があるとtrue
                        for(int l = 0; l < childrenRev.size(); l++) {
                            if(children[k] == childrenRev[l]) {
                                frag = true;
                                break;
                            }
                        }
                        if(frag == false) {
                            childrenRev.push_back(children[k]);
                        }
                    }

                    //==== チェックノードリストchildrenに接続する変数ノードを算出 ====//
                    vector<int> connectedValueList(0);
                    for(int k = 0; k < childrenRev.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[childrenRev[k]][l] == 1) {
                                connectedValueList.push_back(l);
                            }
                        }
                    }

                    //==== connectedValueListの重複を取り除く(ovListとの重複もだめ) ====//
                    vector<int> connectedRev(0);
                    for(int k = 0; k < connectedValueList.size(); k++) {
                        bool frag = false;    // 重複があるとtrue
                        for(int l = 0; l < connectedRev.size(); l++) {
                            if(connectedValueList[k] == connectedRev[l]) {
                                frag = true;
                                break;
                            }
                        }
                        for(int l = 0; l < ovList.size(); l++) {
                            if(connectedValueList[k] == ovList[l]) {
                                frag = true;
                                break;
                            }
                        }
                        if(connectedValueList[k] == i) {
                            frag = true;
                            break;
                        }
                        if(frag == false) {
                            connectedRev.push_back(connectedValueList[k]);
                        }
                    }

                    ovList = connectedRev;
                    if(ovList.empty()) break;    // 次に調べる変数ノードリストが空

                }
                //==== 距離が最も大きいノードリストを算出 ====//
                int maxValue = 1;    // 最大距離(無限大ならば-1)
                for(int j = startPos; j < endPos; j++) {
                    if(distance[j] == -1) {
                        maxValue = -1;
                        break;
                    } else if(distance[j] > maxValue) {
                        maxValue = distance[j];
                    }
                }
                vector<int> highGirthList(0);
                for(int j = startPos; j < endPos; j++) {
                    if(distance[j] == maxValue) {
                        highGirthList.push_back(j);
                    }
                }

                //==== その中で最も次数が小さいチェックノードを算出 ====//
                int cid = highGirthList[0];
                int cdim = nowCheckDim[highGirthList[0]];
                for(int j = 1; j < highGirthList.size(); j++) {
                    if(nowCheckDim[highGirthList[j]] < cdim) {
                        cid = highGirthList[j];    // 更新
                        cdim = nowCheckDim[highGirthList[j]];
                    }
                }

                //==== 辺を設置 ====//
                if((int)retMat[cid][i] == 1) {
                    cout << "distance[" << cid << "] = " << distance[cid] << "\n";
                    cout << "retMat[" << cid << "][" << i << "]に多重辺\n";

                }
                retMat[cid][i] = 1;
                nowCheckDim[cid]++;    // 次数上昇

                vector<int> ovList2(0);
                ovList2.push_back(i);
                ovList = ovList2;

            }
            int count = 0;
            for(int j = 0; j < hr; j++) {
                if((int)retMat[j][i] == 1) {
                    count++;
                }
            }
            if(count != valDim[i]) {
                cout << "valDim[" << i << "] = " << count << "\n";
            }

        }
        return retMat;
    }

    //==== PEG空間結合符号2 ====//
    static Matrix<BinaryFiniteField> PEGSCC2(int rowSize, int colSize, const vector<int>& valDim, int bandWidth) {
        int hr = rowSize;
        int long hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> nowCheckDim(hr);    // 構成中の各チェックノードの現在の次数を表すベクトル
        vector<int> distance(hr);    // 任意の変数ノードの辺接続ステップにおいて、各チェックノードへの距離を一時格納するベクトル。

        //==== チェックノードの次数を0に初期化 ====//
        for(int i = 0; i < hr; i++) nowCheckDim[i] = 0;

        for(int i = 0; i < hc; i++) {
            int selectedVid = i;    // 選択された変数ノード。後に変更することを想定

            //==== i番目の変数ノードから各チェックノードへの距離を無限大(-1)に初期化 ====//
            for(int j = 0; j < hr; j++) distance[j] = -1;

            //==== チェックノードの選択可能範囲(PEG空間結合符号などで利用) ====//
            int c = ((double)(hr - bandWidth + 1) / hc) * i;
            int startPos = c;
            int endPos = c + bandWidth;

            //==== 1本目は次数の一番低いチェックノードに接続 ====//
            int minDimId = startPos;
            int minDim = nowCheckDim[startPos];
            for(int j = startPos + 1; j < endPos; j++) {
                if(nowCheckDim[j] < minDim) {
                    minDimId = j;
                    minDim = nowCheckDim[j];
                }
            }
            retMat[minDimId][i] = 1;    // 1本目の辺を設置
            nowCheckDim[minDimId]++;

            //==== 距離リスト作成(探索)のための初期ノード ====//
            int depth = 0;    // 現在の深さ
            int maxDepth = 10;    // 最大深さ
            vector<int> ovList(0);
            ovList.push_back(i);    // 初期ノードは調査変数ノード

            //==== 変数ノードの次数を満たすまで繰り返し ====//
            for(int j = 1; j < valDim[i]; j++) {
                //==== 最も遠いチェックノードを発見するためのループ ====//
                for(depth = 0; depth < maxDepth; depth++) {
                    //==== チェックノードへの距離を更新 ====//
                    vector<int> children(0);
                    for(int l = 0; l < ovList.size(); l++) {
                        for(int k = 0; k < hr; k++) {
                            if((int)retMat[k][ovList[l]] == 1) {    // 接続しているチェックノードに対して
                                if(distance[k] == -1) {    // 距離が未設定(無限大)→更新
                                    children.push_back(k);
                                    distance[k] = depth;    // 距離を更新
                                }else if(distance[k] != -1 && depth < distance[k]) {    // 現在の距離(depth)より大きい距離が設定されている→更新
                                    children.push_back(k);
                                    distance[k] = depth;    // 距離を更新
                                }
                            }
                        }
                    }
                    //==== チェックノードリストchildrenの重複を取り除く ====//
                    vector<int> childrenRev(0);    // 重複を取り除いたchildren
                    for(int k = 0; k < children.size(); k++) {
                        bool frag = false;    // 重複があるとtrue
                        for(int l = 0; l < childrenRev.size(); l++) {
                            if(children[k] == childrenRev[l]) {
                                frag = true;
                                break;
                            }
                        }
                        if(frag == false) {
                            childrenRev.push_back(children[k]);
                        }
                    }

                    //==== チェックノードリストchildrenに接続する変数ノードを算出 ====//
                    vector<int> connectedValueList(0);
                    for(int k = 0; k < childrenRev.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[childrenRev[k]][l] == 1) {
                                connectedValueList.push_back(l);
                            }
                        }
                    }

                    //==== connectedValueListの重複を取り除く(ovListとの重複もだめ) ====//
                    vector<int> connectedRev(0);
                    for(int k = 0; k < connectedValueList.size(); k++) {
                        bool frag = false;    // 重複があるとtrue
                        for(int l = 0; l < connectedRev.size(); l++) {
                            if(connectedValueList[k] == connectedRev[l]) {
                                frag = true;
                                break;
                            }
                        }
                        for(int l = 0; l < ovList.size(); l++) {
                            if(connectedValueList[k] == ovList[l]) {
                                frag = true;
                                break;
                            }
                        }
                        if(connectedValueList[k] == i) {
                            frag = true;
                            break;
                        }
                        if(frag == false) {
                            connectedRev.push_back(connectedValueList[k]);
                        }
                    }

                    ovList = connectedRev;
                    if(ovList.empty()) break;    // 次に調べる変数ノードリストが空

                }
                //==== 距離が最も大きいノードリストを算出 ====//
                int maxValue = 1;    // 最大距離(無限大ならば-1)
                for(int j = startPos; j < endPos; j++) {
                    if(distance[j] == -1) {
                        maxValue = -1;
                        break;
                    } else if(distance[j] > maxValue) {
                        maxValue = distance[j];
                    }
                }
                vector<int> highGirthList(0);
                for(int j = startPos; j < endPos; j++) {
                    if(distance[j] == maxValue) {
                        highGirthList.push_back(j);
                    }
                }

                //==== その中で最も次数が小さいチェックノードを算出 ====//
                int cid = highGirthList[0];
                int cdim = nowCheckDim[highGirthList[0]];
                for(int j = 1; j < highGirthList.size(); j++) {
                    if(nowCheckDim[highGirthList[j]] < cdim) {
                        cid = highGirthList[j];    // 更新
                        cdim = nowCheckDim[highGirthList[j]];
                    }
                }

                //==== 辺を設置 ====//
                if((int)retMat[cid][i] == 1) {
                    cout << "distance[" << cid << "] = " << distance[cid] << "\n";
                    cout << "retMat[" << cid << "][" << i << "]に多重辺\n";

                }
                retMat[cid][i] = 1;
                nowCheckDim[cid]++;    // 次数上昇

                vector<int> ovList2(0);
                ovList2.push_back(i);
                ovList = ovList2;

            }
            int count = 0;
            for(int j = 0; j < hr; j++) {
                if((int)retMat[j][i] == 1) {
                    count++;
                }
            }
            if(count != valDim[i]) {
                cout << "valDim[" << i << "] = " << count << "\n";
            }

        }
        return retMat;
    }

    //==== LRLを結合して新しいプロトグラフを構成 ====//
    // lrl1 : 外側のLRL空間結合符号プロトグラフを指定
    // lrl2 : 内側のLRL空間結合符号プロトグラフを指定
    // splitId : 結合箇所の変数ノードid
    static Matrix<BinaryFiniteField> connectLRL(Matrix<BinaryFiniteField>* lrl1, Matrix<BinaryFiniteField>* lrl2, unsigned long splitId) {
        unsigned long hr1 = lrl1->row();
        unsigned long hr2 = lrl2->row();
        unsigned long hc1 = lrl1->col();
        unsigned long hc2 = lrl2->col();
        if(splitId < 0 || splitId > hc1) throw ConstructExeption();

        //==== LRL1の切り口の左右のデータを取得 ====//
        // 左側
        BandData left1 = {0, 0, 0};
        if(splitId != 0) {
            left1 = getBandData(lrl1, splitId - 1);
        }
        // 右側
        BandData right1 = {0, 0, 0};
        if(splitId != hc1) {
            right1 = getBandData(lrl1, splitId);
        }

        //==== 接続するLRL2の切り口のデータを取得 ====//
        // 左側
        BandData left2 = getBandData(lrl2, 0);
        // 右側
        BandData right2 = getBandData(lrl2, hc2 - 1);

    }

    //==== Loop空間結合符号 ====//
    static Matrix<BinaryFiniteField> loopLDPC(const vector<LoopLDPC::LBandData>& bandData) {
        LoopLDPC lldpc(bandData);
        return lldpc.generate();
    }

    //==== 検査行列ファイルを読み込む ====//
    static Matrix<BinaryFiniteField> loadMatrixFile(const string fileName) {
        SPMatrix H(fileName);
        Matrix<BinaryFiniteField> retMat = H;

        //==== 叉状空間結合符号のテスト(デバッグ用) ====//
        // ここから
        /*unsigned long hr = 31;
        unsigned long hc = 58;
        Matrix<BinaryFiniteField> PosMat(hr, hc);*/

        // 近いところは近くと接続
        /*unsigned long d = 2;
        for(unsigned long i = 0; i < d; i++) {
            PosMat[30 - i][56 - i * 2] = 1;
            PosMat[30 - i][57 - i * 2] = 1;
        }*/

        // 遠いところほど近くと接続
        /*unsigned long e = 2;
        for(unsigned long i = 0; i < e; i++) {
            PosMat[30 - i][58 - e * 2 + i * 2] = 1;
            PosMat[30 - i][59 - e * 2 + i * 2] = 1;
        }*/

        // 三角形塗りつぶし
        /*unsigned long f = 10;
        for(unsigned long i = 0; i < f; i++) {
            for(unsigned long j = 0; j < f * 2 - i * 2; j++) {
                PosMat[30 - i][58 - f * 2 + i * 2 + j] = 1;
            }
        }*/

        // EXチェックノードを伸ばす接続
        /*unsigned long g = 29;
        for(unsigned long i = 0; i < g; i++) {
            PosMat[30][58 - g * 2 + i * 2] = 1;
            PosMat[30][59 - g * 2 + i * 2] = 1;
        }*/

        // EX変数ノードを伸ばす接続
        /*unsigned long h = 29;
        for(unsigned long i = 0; i < h; i++) {
            PosMat[30 - i][56] = 1;
            PosMat[30 - i][57] = 1;
        }*/

        // 中央部結合
        /*PosMat[14][28] = 1;
        PosMat[14][29] = 1;
        PosMat[15][26] = 1;
        PosMat[15][27] = 1;*/

        /*for(unsigned long i = 0; i < hr; i++) {
            for(unsigned long j = 0; j < hc; j++) {
                retMat[i][j + hc] = PosMat[i][j];
                retMat[i + hr][j + hc * 2] = PosMat[i][j];
                retMat[i + hr * 2][j] = PosMat[i][j];
            }
        }*/
        // ここまで

        return retMat;
    }

    //==== プロトグラフ叉状空間結合符号の作成 ====//
    // 独立ケースプロトグラフの行列ファイルを読み込み、接続部を追加
    // fileName : [string]読み込むファイルのパス＋名前
    // l : [int]列重み
    // r : [int]最大行重み
    // L : [int]結合数
    // copyNum : [int]鎖の本数
    // connectType : [int]接続方法
    // connectSize : [int]接続部の大きさ
    static Matrix<BinaryFiniteField> ConnectIndependentLRLProtograph (int l, int r, int L, int copyNum, int connectSize) {
        //Matrix<BinaryFiniteField> retMat = SMConstructor::loadMatrixFile(fileName);
        unsigned long protHr = l + L - 1;    // lrLプロトグラフの行数
        unsigned long protHc = r * L / l;    // lrLプロトグラフの列数
        unsigned long hr = protHr * copyNum;    // retMatの行数
        unsigned long hc = protHc * copyNum;    // retMatの列数
        Matrix<BinaryFiniteField> retMat(hr, hc);

        // 一貫性エラーチェック
        /*if(hr != retMat.row() || hc != retMat.col()) {    // 指定した内容とプロトグラフが異なる
            cout << "[SMConstructor::loadAndConnectIndependentLRLProtograph][error]読み込んだプロトグラフの行数・列数と、指定されたパラメータ(l,r,L,copyNum)から算出される行数・列数が一致しません。プロトグラフの指定が誤っているか、パラメータを正しく指定してください。\n";
            Matrix<BinaryFiniteField> nu(0,0);
            return nu;
        }
        if(connectSize > L) {
            cout << "[SMConstructor::loadAndConnectIndependentLRLProtograph][error]接続部の大きさconnectSizeは、結合数L以下でなければなりません。\n";
            Matrix<BinaryFiniteField> nu(0,0);
            return nu;
        }*/

        // 各ブロック基準位置計算
        // ※注意：接続部Cの基準位置は以下の通り
        // a番目のプロトグラフを(a + 1 mod copyNum)番目のプロトグラフと接続する場合、接続部Cは
        // 基準行：rowStartPow[a] ... rowStartPow[a + 1] - 1
        // 基準列：colStartPow[(a + 1) % copyNum] ... colStartPow[(a + 1) % copyNum + 1] - 1
        // となる点に注意。
        vector<unsigned long> rowStartPos(copyNum + 1, 0);
        vector<unsigned long> colStartPos(copyNum + 1, 0);
        for(int i = 0 ; i < copyNum + 1; i++) {
            rowStartPos[i] = protHr * i;
            colStartPos[i] = protHc * i;
        }

        // lrL空間結合符号のプロトグラフを配置
        for(int i = 0; i < copyNum; i++) {
            int baseX = rowStartPos[i];
            int baseY = colStartPos[i];
            int bHeight = l;
            int bWidth = r / l;    // 各ブロックの横幅
            for(int j = 0; j < L; j++) {
                int bbaseX = baseX + j;
                int bbaseY = baseY + (r * j / l);
                for(int x = 0; x < bHeight; x++) {
                    for(int y = 0; y < bWidth; y++) {
                        retMat[bbaseX + x][bbaseY + y] = 1;
                    }
                }
            }
        }

        // 接続部の作成
        // i番目のブロックと、(i + 1) % copyNum 番目のブロックを接続
        for(int i = 0; i < copyNum; i++) {
            // i番目の接続部の配置位置決定
            unsigned long baseRow = rowStartPos[i];
            unsigned long baseCol = colStartPos[(i + 1) % copyNum];
            unsigned long baseRowEnd = baseRow + protHr;
            unsigned long baseColEnd = baseCol + protHc;

            int e = connectSize;
            for(int j = 0; j < e; j++) {
                for(int k = 0; k < r / l; k++) {
                    retMat[baseRowEnd - 1 - j][baseColEnd - e * (r/l) + j * (r/l) + k] = 1;
                }
            }
        }


        return retMat;
    }

    //==== プロトグラフ叉状空間結合符号の作成(結合数に応じて接続部を自動的に大きくするVersion) ====//
    // ※注意：LoopLDPCメソッドとは異なり、ファイルから読み込む点に注意すること
    // rate : [double]接続部の大きさのLに対する割合(L=8でrate=0.5なら接続部の大きさは4)
    static Matrix<BinaryFiniteField> ConnectIndependentLRLProtographWithProportionalConnection (int l, int r, int L, int copyNum, double rate) {
        return ConnectIndependentLRLProtograph(l, r, L, copyNum, (int)(L * rate));
    }


private :
    typedef struct ___BANDDATA {
        unsigned long top;
        unsigned long bottom;
        unsigned long bandWidth;
    } BandData;

    // 指定した位置の帯の最上部・最下部・縦幅を取得する関数
    static BandData getBandData(const Matrix<BinaryFiniteField>* mat, unsigned long variableId) {
        BandData retObj = {0, 0, 0};
        unsigned long hr = mat->row();
        unsigned long hc = mat->col();
        if(variableId < 0 || variableId >= hc) throw ConstructExeption();
        unsigned long j = variableId;
        for(unsigned long i = 0; i < hr; i++) {
            if((int)(*mat)[i][j] == 1) {
                retObj.top = j;
                break;
            }
        }
        for(unsigned long i = hr - 1; i >= 0; i--) {
            if((int)(*mat)[i][j] == 1) {
                retObj.bottom = j;
                break;
            }
        }
        if(retObj.top > retObj.bottom) retObj.top = retObj.bottom;
        retObj.bandWidth = retObj.bottom - retObj.top;

        return retObj;
    }

    // b-a個の中からn個ランダムに取ったベクトルを生成
    static vector<int> getRandomSequence(int a, int b, int n) {
        if(a > b) {
            int temp = a;
            a = b;
            b = temp;
        }
        int c = b - a;
        if(n > c) n = c;    // nは最大b-cまで
        vector<int> retVec(n);
        RandomValue rv;
        //cout << "bandWidth = " << c << "\n";
        for(int i = 0; i < n; i++) {
            retVec[i] = rv.getRand(c - i) + a;
            //if(retVec[i] < a) cout << "おかしい\n";
        }
        for(int i = n - 1; i >= 0; i--) {
            for(int j = i - 1; j >= 0; j--) {
                if(retVec[i] >= retVec[j]) {
                    retVec[i]++;
                    //if(retVec[i] >= b) cout << "おかしい\n";
                }
            }
        }
        std::sort(retVec.begin(), retVec.end());
        return retVec;
    }

    static void _matrixIdentBlockPut(int row, int col, int paramP, Matrix<BinaryFiniteField>* M) {    // 単位行列を設置
        for(int i = 0; i < paramP; i++) {
            (*M)[i + (row * paramP)][i + (col * paramP)] = 1;
        }
    }
    static void _matrixShiftBlockPut(int row, int col, int paramP, int shiftNum, Matrix<BinaryFiniteField>* M) {
        if(shiftNum > 0) shiftNum %= paramP;
        else {
            shiftNum = (-1) * (((-1) * shiftNum) % paramP);
            shiftNum += paramP;
        }
        for(int i = 0; i < paramP; i++) {
            if(i < shiftNum) {
                (*M)[(paramP - shiftNum) + i + (row * paramP)][i + (col * paramP)] = 1;
            } else {
                (*M)[i - shiftNum + (row * paramP)][i + (col * paramP)] = 1;
            }
        }
    }

    // 置換ベクトルを取得
    static inline vector<int> getPerm(unsigned long size) {
        vector<int> retVec;
        retVec.reserve(size);
        vector<int> tmpVec(size);
        for(unsigned long i = 0; i < size; i++) {
            tmpVec[i] = i;
        }
        for(unsigned long i = 0; i < size; i++) {
            RandomValue rv;
            unsigned long r = rv.getRand(size - i);
            retVec.push_back(tmpVec[r]);

            vector<int> tmp2;
            tmp2.reserve(size - i - 1);
            for(unsigned long j = 0; j < r; j++) {
                tmp2.push_back(tmpVec[j]);
            }
            for(unsigned long j = r + 1; j < size - i; j++) {
                tmp2.push_back(tmpVec[j]);
            }
            tmpVec = tmp2;
        }
        return retVec;
    }

    // vectorの値にvalが含まれているかどうかを調べ、なければ追加しソートする
    static void contains(vector<int>& vec, int val) {
        int vs = vec.size();
        int check = 0;
        for(int i = 0; i < vs; i++) {
            if(vec[i] == val) {
                check = 1;
                break;
            }
        }
        if(check == 0) vec.push_back(val);

        std::sort(vec.begin(), vec.end());
    }

    static int eqVec(const vector<int>& A, const vector<int>& B) {
        int as = A.size();
        int bs = B.size();
        if(as != bs) return 0;
        for(int i = 0; i < as; i++) {
            if(A[i] != B[i]) return 0;
        }
        return 1;
    }
    static int contains2(const vector<int>& vec, int val) {
        for(int i = 0; i < vec.size(); i++) {
            if(vec[i] == val) return 1;
        }
        return 0;
    }
    // vectorの中にmin以上max未満の値が全て含まれている場合は1を返す
    inline static int containMinToMax(const vector<int>& vec, int min, int max) {
        int counter = 0;
        for(int i = min; i < max; i++) {
            if(contains2(vec, i)) counter++;
        }
        if(counter == max - min) {
            return 1;
        } else return 0;
    }

    // 補集合を取得
    static vector<int> complement(const vector<int>& vec, int size) {
        vector<int> retVec(0);
        int vsize = vec.size();
        for(int i = 0; i < size; i++) {
            int check = 0;
            for(int j = 0; j < vsize; j++) {
                if(vec[j] == i) check = 1;
            }
            if(check == 0) retVec.push_back(i);
        }
        return retVec;
    }

public:
    class ConstructExeption{};    // 例外クラス

};

#endif
