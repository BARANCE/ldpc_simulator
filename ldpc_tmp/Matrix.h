/*
    行列を扱うクラステンプレートMatrix
    【機能】
    ・テンプレートなので、任意の可算クラスを行列の成分として扱える
    ・ファイルからの行列オブジェクト作成
    ・単位行列の生成
    ・行列同士の加減乗算
    ・定数倍
    ・行列の切り出しや、行列同士の接合
    ・転置行列
    ・列置換、行置換
*/


#if !defined(___Class_Matrix)
#define ___Class_Matrix

//#include "MatrixFileOperator.h"
//#include "IOOError.h"    // エラー処理
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include "RandomValue.h"
#include "BinaryFiniteField.h"

using namespace std;

template < class Type > class MOperationErr;


//==== 行列クラステンプレート =====//
template <class Type> class Matrix {
    unsigned long rowSize;    // 行列の行数
    unsigned long colSize;    // 行列の列数
    /*int maxRowSize;    // 確保済みの行数
    int maxColSize;    // 確保済みの列数
    Type* initializer;    // 行列の初期化子
    bool initExist;    // initializerが存在するかどうか*/

    vector<vector<Type> >* mat;    // 行列本体
public :
    //===== コンストラクタ =====//
    Matrix() : rowSize(0), colSize(0){    // デフォルトコンストラクタ
        regionSet(rowSize, colSize);
    }
    Matrix(unsigned long r, unsigned long c) : rowSize(r), colSize(c) {    // r×c行列を生成
        regionSet(rowSize, colSize);
    }
    Matrix(const Matrix& x) : rowSize(x.row()), colSize(x.col()) {    // コピーコンストラクタ
        regionSet(rowSize, colSize);
        //cout << "rowSize = " << rowSize << ", colSize = " << colSize << "\n";
        matCopy(x);
    }
    explicit Matrix(const string fileName) {    // ファイルから
        /*MatrixFileOperator fo(fileName);
        Matrix<Type> nMat = fo.readMatrix<Type>();
        fo.close();
        rowSize = nMat.row();
        colSize = nMat.col();
        regionSet(rowSize, colSize);
        matCopy(nMat);*/

        ifstream ifs(fileName.c_str(), fstream::in);
        if (ifs.fail()) throw IOExceptions<Type>(this);
        ifs >> rowSize >> colSize;
        if (ifs.eof()) throw FileFormatErr<Type>(this, "FileFormatErr:ファイルの行数または列数が読み込めません。\n");
        regionSet(rowSize, colSize);
        for(int i = 0; i < rowSize; i++) {
            for(int j = 0; j < colSize; j++) {
                if (ifs.eof()) throw FileFormatErr<Type>(this, "FileFormatErr:指定された行数・列数に達する前にファイルの終端に達しました。\n");
                ifs >> (*mat)[i][j];
            }
        }
        ifs.close();
    }

    template <class Type2> class IOExceptions {
    public:
        IOExceptions(const Matrix<Type2>* ident) {}
    };
    template <class Type2> class FileFormatErr {
    public:
        FileFormatErr(const Matrix<Type2>* ident) {}
        FileFormatErr(const Matrix<Type2>* ident, const string str) {
            cout << str;
        }
    };

    //===== 変換コンストラクタ =====//
    Matrix(const vector<Type> x) : rowSize(1), colSize(x.size()) {
        regionSet(rowSize, colSize);
        vecCopy(x);
    }

    //===== デストラクタ =====//
    ~Matrix() {
        delete mat;
        mat = NULL;
    }

    //===== 値取得用メソッド =====//
    unsigned long row() const {return rowSize;};    // 行数を返す
    unsigned long col() const {return colSize;};    // 列数を返す
    stringstream getStringstream() const {    // 他関数で使用する文字列ストリーム取得用関数
        stringstream ss;
        unsigned long i, j;
        for(i = 0; i < rowSize; i++) {
            for(j = 0; j < colSize; j++) {
                ss << (*mat)[i][j];
                if(j < colSize - 1) ss << " ";
            }
            ss << "\n";
        }
        return ss;
    }
    string toString() const {return getStringstream().str();}
    void print() const {cout << toString();}

    //===== 行列を切り出す =====//
    // index行を抜き出す
    Matrix<Type> getRow(unsigned long index) const {
        Matrix<Type> valMat(1, colSize);
        for(unsigned long i = 0; i < colSize; i++) {
            valMat[0][i] = (*mat)[index][i];
        }
        //cout << "rowSize = " << valMat.row() << ", colSize = " << valMat.col() << "\n";
        return valMat;
    }
    // index列を抜き出す
    Matrix<Type> getCol(unsigned long index) const {
        Matrix<Type> valMat(rowSize, 1);
        for(unsigned long i = 0; i < rowSize; i++) {
            valMat[i][0] = (*mat)[i][index];
        }
        return valMat;
    }
    // 部分行列を得る
    Matrix<Type> getSubMatrix(unsigned long row1, unsigned long row2, unsigned long col1, unsigned long col2) const {
        if(row1 > row2) {
            unsigned long tmp = row1;
            row1 = row2;
            row2 = tmp;
        }
        if(col1 > col2) {
            unsigned long tmp = col1;
            col1 = col2;
            col2 = tmp;
        }
        Matrix<Type> valMat(row2 - row1, col2 - col1);
        for(unsigned long i = row1; i < row2; i++) {
            for(unsigned long j = col1; j < col2; j++) {
                valMat[i - row1][j - col1] = (*mat)[i][j];
            }
        }
        return valMat;
    }

    //==== 行または列の削除 ====//
    Matrix<Type> deleteRow(unsigned long index) const {    // index行を削除する
        Matrix<Type> valMat(rowSize - 1, colSize);
        for(unsigned long i = 0; i < index; i++) {    // index行の前まで
            for(unsigned long j = 0; j < colSize; j++) {
                valMat[i][j] = (*mat)[i][j];
            }
        }
        for(unsigned long i = index + 1; i < rowSize; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                valMat[i - 1][j] = (*mat)[i][j];
            }
        }
        return valMat;
    }

    Matrix<Type> deleteCol(unsigned long index) const {    // index列を削除する
        Matrix<Type> valMat(rowSize, colSize - 1);
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < index; j++) {
                valMat[i][j] = (*mat)[i][j];
            }
            for(unsigned long j = index + 1; j < colSize; j++) {
                valMat[i][j - 1] = (*mat)[i][j];
            }
        }
        return valMat;
    }
    Matrix<Type> deleteRowCol(unsigned long row, unsigned long col) const {    // row行とcol列を削除
        Matrix<Type> valMat(rowSize - 1, colSize - 1);
        if(row >= rowSize) cout << "添字エラー\n";
        if(col >= colSize) cout << "添字エラー\n";
        for(unsigned long i = 0; i < row; i++) {
            for(unsigned long j = 0; j < col; j++) {
                valMat[i][j] = (*mat)[i][j];
            }
            for(unsigned long j = col + 1; j < colSize; j++) {
                valMat[i][j - 1] = (*mat)[i][j];
            }
        }
        for(unsigned long i = row + 1; i < rowSize; i++) {
            for(unsigned long j = 0; j < col; j++) {
                valMat[i - 1][j] = (*mat)[i][j];
            }
            for(unsigned long j = col + 1; j < colSize; j++) {
                valMat[i - 1][j - 1] = (*mat)[i][j];
            }
        }
        return valMat;
    }

    //===== 行列を変換する =====//
    // 転置行列を得る
    Matrix<Type> transposed() {
        Matrix<Type> valMat(colSize, rowSize);
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                valMat[j][i] = (*mat)[i][j];
            }
        }
        return valMat;
    }

    //==== 行列の接合 ====//
    // 上方向に行列を接合する
    Matrix<Type> unionTop(const Matrix<Type>& M) {
        unsigned long mr = M.row();
        Matrix<Type> valMat(rowSize + mr, colSize);
        for(unsigned long i = 0; i < mr; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                valMat[i][j] = M[i][j];
            }
        }
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                valMat[i + mr][j] = (*mat)[i][j];
            }
        }
        return valMat;
    }

    // 下方向に行列を接合する
    Matrix<Type> unionBottom(const Matrix<Type>& M) {
        unsigned long mr = M.row();
        Matrix<Type> valMat(rowSize + mr, colSize);
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                valMat[i][j] = (*mat)[i][j];
            }
        }
        for(unsigned long i = 0; i < mr; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                valMat[i + rowSize][j] = M[i][j];
            }
        }
        return valMat;
    }
    void unionBottomOwn(const Matrix<Type>& M) {
        unsigned long mr = M.row();
        vector<vector<Type> >* matTmp = new vector<vector<Type> >(rowSize + mr, vector<Type>(colSize));
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                (*matTmp)[i][j] = (*mat)[i][j];
            }
        }
        for(unsigned long i = 0; i < mr; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                (*matTmp)[i + rowSize][j] = M[i][j];
            }
        }
        rowSize = rowSize + mr;
        delete mat;
        mat = matTmp;
    }
    // 左方向に行列を接合する
    Matrix<Type> unionLeft(const Matrix<Type>& M) {
        unsigned long mc = M.col();
        Matrix<Type> valMat(rowSize, colSize + mc);
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < mc; j++) {
                valMat[i][j] = M[i][j];
            }
        }
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                valMat[i][j + mc] = (*mat)[i][j];
            }
        }
        return valMat;
    }
    // 右方向に行列を接合する
    Matrix<Type> unionRight(const Matrix<Type>& M) {
        unsigned long mc = M.col();
        Matrix<Type> valMat(rowSize, colSize + mc);
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                valMat[i][j] = (*mat)[i][j];
            }
        }
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < mc; j++) {
                valMat[i][j + colSize] = M[i][j];
            }
        }
        return valMat;
    }

    //==== 行または列の置換(交換) ====//
    // col1列とcol2列を列置換
    void colSubstitution(const unsigned long col1, const unsigned long col2) {
        vector<Type> tmpVec(rowSize);
        for(unsigned long i = 0; i < rowSize; i++) {
            tmpVec[i] = (*mat)[i][col1];
            (*mat)[i][col1] = (*mat)[i][col2];
            (*mat)[i][col2] = tmpVec[i];
        }
    }
    // row1行とrow2行を行置換
    void rowSubstitution(const unsigned long row1, const unsigned long row2) {
        vector<Type> tmpVec(colSize);
        for(unsigned long j = 0; j < colSize; j++) {
            tmpVec[j] = (*mat)[row1][j];
            (*mat)[row1][j] = (*mat)[row2][j];
            (*mat)[row2][j] = tmpVec[j];
        }
    }
    // ランダム列置換を1回だけ行う
    void randomColSubstitution() {
        RandomValue rv;
        unsigned long r1 = rv.getRand(colSize);
        unsigned long r2 = rv.getRand(colSize);
        while(r1 == r2) {r2 = rv.getRand(colSize);}
        //cout << "r1 = " << r1 << ", r2 = " << r2 << "\n";
        colSubstitution(r1, r2);
    }

    //==== 行基本変形 ====//
    // A行目をB行目に足す
    void AddRowFromAToB(const unsigned long& A, const unsigned long& B) {
        for(unsigned long j = 0; j < colSize; j++) {
            (*mat)[B][j] += (*mat)[A][j];
        }
    }
    // A行目をB行目から引く
    void SubRowFromAToB(const unsigned long& A, const unsigned long& B) {
        for(unsigned long j = 0; j < colSize; j++) {
            (*mat)[B][j] -= (*mat)[A][j];
        }
    }

    //==== 行列の標準形を取得(検査行列版) ====//
    // ※行基本変形のみで標準形に至らない場合は、列置換も行う
    Matrix<Type> getRENF_H() {
        Matrix<Type> valMat = (*this);    // 生成行列
        unsigned long hr = rowSize;
        unsigned long hc = colSize;
        bool deleteCheck = false;
        for(unsigned long s = 0; s < hr; s++) {
            unsigned long n = hc - hr;
            //cout << s << "行" << s+n << "列目の処理：";
            if((int)valMat[s][s + n] == 0) {
                int checker = 0;
                for(unsigned long i = s + 1; i < hr; i++) {
                    //cout << valMat[i][s + n] << "\n";
                    if((int)valMat[i][s + n] != 0) {
                        valMat.rowSubstitution(s, i);    // [s][s+n]要素が1となるように行置換
                        checker = 1;
                        break;
                    }
                }
                if(checker == 0) {    // 列置換が必要
                    //cout << "[*info*] " << s+n << "列目には1のある行がありません。";
                    //cout << "\n";
                    for(unsigned long j = 0; j < n; j++) {
                        if((int)valMat[s][j] != 0) {
                            valMat.colSubstitution(j, s + n);    // [s][s+n]要素が1となるように列置換
                            (*this).colSubstitution(j, s + n);
                            //cout << "[*info*] " << j << "列目と" << s + n << "列目を置換しました。検査行列の列置換を行いました。\n";
                            checker = 2;
                            break;
                        }
                    }
                    if(checker == 0) {    // 列置換もできない
                        //cout << "[*info*] " << s << "行目が全零です。標準形に変形できません。\n";
                        //valMat[s][s+n] = 1;
                        //valMat.exportMatrix("tmp1.txt");
                        valMat = valMat.deleteRowCol(s, s+n);
                        hr--;
                        hc--;
                        n = hc - hr;
                        //valMat.exportMatrix("tmp2.txt");
                        deleteCheck = true;
                        //cout << "[*info*] " << s << "行、および" << s+n << "列目を削除しました。\n";
                    }
                } else {}//cout << "\n";
            }
            if(deleteCheck == false && (int)valMat[s][s + n] != 0) {
                // s+1行～rowSize行でs+n列目に1が立っている行にs行目を足す
                for(unsigned long j = 0; j < s; j++) {
                    if((int)valMat[j][s + n] != 0) valMat.AddRowFromAToB(s, j);
                }
                for(unsigned long j = s+1; j < rowSize; j++) {
                    if((int)valMat[j][s + n] != 0) valMat.AddRowFromAToB(s, j);
                }

            }
            //cout << "s = " << s << "\n";
        }
        //cout << valMat << "\n";

        return valMat;
    }

    //===== オブジェクト操作用メソッド =====//
    void regionSet(const unsigned long r, const unsigned long c) {    // r×c行列を設置(既存のものは消去)
        try {
            mat = new vector<vector<Type> >(r, vector<Type>(c));
        } catch(std::bad_alloc ba) {    // 領域確保失敗
            //throw MatrixBadAlloc(&ba, r, c);
            throw this;
        }
    }
    void matCopy(const Matrix<Type>& M) {    // Mをコピー(既存のものは消去)
        unsigned long i, j;
        //cout << "rowSize = " << rowSize << ", colSize = " << colSize << "\n";
        for(i = 0; i < rowSize; i++) {
            for(j = 0; j < colSize; j++) {
                (*mat)[i][j] = M[i][j];
            }
        }
    }
    void vecCopy(const vector<Type>& V) {
        unsigned long i;
        unsigned long size = V.size();
        for(i = 0; i < size; i++) {
            (*mat)[0][i] = V[i];
        }
    }
    void reserve(unsigned long r, unsigned long c) {    // r×c行列の領域を確保(既存のものを維持)
        vector<vector<Type> >* tmpMat;    // コピー処理用の配列
        Type* a;
        tmpMat = getMatrixRegion(r, c, a);
        for(unsigned long i = 0; i < rowSize || i < r; i++) {
            for(unsigned long j = 0; j < colSize || j < c; j++) {
                (*tmpMat)[i][j] = (*mat)[i][j];    // 行列コピー
            }
        }
        rowSize = r;
        colSize = c;
        delete mat;
        mat = NULL;
        mat = tmpMat;
    }
    // 単位行列化(既存のものは消去)
    void setIdentity() {
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                if(i == j) (*mat)[i][j] = 1;
                else (*mat)[i][j] = 0;
            }
        }
    }

    // ファイルに保存
    void exportMatrix(string fileName) {
        ofstream ofs(fileName.c_str(), fstream::out);
        if(ofs.fail()) throw IOExceptions<Type>(this);
        ofs << rowSize << " " << colSize << "\n";
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                ofs << (*mat)[i][j];
                if(j < colSize - 1) ofs << " ";
            }
            ofs << "\n";
        }
        ofs.close();
    }

    //===== 演算子の二重定義[1] =====//
    Matrix<Type>& operator=(const Matrix<Type>& M) {    // 代入演算子
        if(&M != this) {    // 自分自身では無い
            rowSize = M.row();
            colSize = M.col();
            delete mat;
            mat = NULL;
            regionSet(rowSize, colSize);
            matCopy(M);
        }
        return *this;
    }
    vector<Type>& operator[](const int i) const {return (*mat)[i];}    // 添字演算子
    vector<Type>& operator()(const int i) const {return (*mat)[i];}    // 上と同じ
    Type& operator()(const int i, const int j) const {return (*mat)[i][j];}    // アクセス
};

//===== クラスで使う関数 =====//
// rowSize×colSizeの行列領域を確保する
// 型Typeにデフォルトコンストラクタが無い場合は未定義
// rowSize : 行数
// colSize : 列数
// a : 型推測用変数(内部値は使用しない)
template <class Type>
vector<vector<Type> >* getMatrixRegion(const unsigned long rowSize, const unsigned long colSize, const Type* a) {
    vector<vector<Type> >* mat;
    try{
        mat = new vector<vector<Type> >(rowSize, vector<Type>(colSize));
    } catch(std::bad_alloc ba) {    // 領域確保失敗
        cout << "メモリエラーが発生しました。";
        /* throw MatrixBadAlloc(&ba, rowSize, colSize); */
        throw ba;
    }
    return mat;
}

//===== 演算子の二重定義[2] =====//
// <<演算子
template <class Type> ostream& operator<<(ostream& ss, const Matrix<Type>& M) {return ss << M.toString();}
// 加算
template <class Type> Matrix<Type> operator+(const Matrix<Type>& M, const Matrix<Type>& N) {
    unsigned long mr = M.row();
    unsigned long mc = M.col();
    unsigned long nr = N.row();
    unsigned long nc = N.col();
    if(mr != nr || mc != nc) throw MOperationErr<Type>(&M, &N);
    Matrix<Type> valMat(mr, mc);    // Mと同じ大きさの行列
    for(unsigned long i = 0; i < mr; i++) {
        for(unsigned long j = 0; j < mc; j++) {
            valMat[i][j] = M[i][j] + N[i][j];
        }
    }
    return valMat;
}
// 加算
template <class Type> Matrix<Type> operator+=(Matrix<Type>& M, const Matrix<Type>& N) {
    unsigned long mr = M.row();
    unsigned long mc = M.col();
    unsigned long nr = N.row();
    unsigned long nc = N.col();
    if(mr != nr || mc != nc) throw MOperationErr<Type>(&M, &N);
    for(unsigned long i = 0; i < mr; i++) {
        for(unsigned long j = 0; j < mc; j++) {
            M[i][j] = M[i][j] + N[i][j];
        }
    }
    return M;
}
// 減算
template <class Type> Matrix<Type> operator-(const Matrix<Type>& M, const Matrix<Type>& N) {
    return M + (-1) * N;
}
// 減算
template <class Type> Matrix<Type> operator-(Matrix<Type>& M, const Matrix<Type>& N) {
    unsigned long mr = M.row();
    unsigned long mc = M.col();
    unsigned long nr = N.row();
    unsigned long nc = N.col();
    if(mr != nr || mc != nc) throw MOperationErr<Type>(&M, &N);
    for(unsigned long i = 0; i < mr; i++) {
        for(unsigned long j = 0; j < mc; j++) {
            M[i][j] = M[i][j] - N[i][j];
        }
    }
    return M;
}
// 乗算
template <class Type> Matrix<Type> operator*(const Matrix<Type>& M, const Matrix<Type>& N) {
    unsigned long mr = M.row();
    unsigned long mc = M.col();
    unsigned long nr = N.row();
    unsigned long nc = N.col();
    if(mc != nr) throw MOperationErr<Type>(&M, &N);
    Matrix<Type> valMat(mr, nc);
    for(unsigned long i = 0; i < mr; i++) {
        for(unsigned long j = 0; j < nc; j++) {
            valMat[i][j] = 0;
            for(unsigned long k = 0; k < mc; k++) {
                valMat[i][j] += M[i][k] * N[k][j];
            }
        }
    }
    return valMat;
}
template <class Type> Matrix<Type> operator*(const Matrix<Type>& M, const double n) {
    unsigned long mr = M.row();
    unsigned long mc = M.col();
    Matrix<Type> valMat(mr, mc);
    for(unsigned long i = 0; i < mr; i++) {
        for(unsigned long j = 0; j < mc; j++) {
            valMat[i][j] = M[i][j] * n;
        }
    }
    return valMat;
}
template <class Type> Matrix<Type> operator*(const double m, const Matrix<Type>& N) {return N * m;}
// 乗算
template <class Type> Matrix<Type> operator*=(Matrix<Type>& M, const Matrix<Type>& N) {
    Matrix<Type> valMat(M.row(), N.col());
    valMat = M * N;
    M = valMat;
    return M;
}
template <class Type> Matrix<Type> operator*=(Matrix<Type>& M, const double n) {
    Matrix<Type> valMat(M.row(), M.col());
    valMat = M * n;
    M = valMat;
    return M;
}
// 等価
template <class Type> int operator==(const Matrix<Type>& M, const Matrix<Type>& N) {
    unsigned long mr = M.row();
    unsigned long mc = M.col();
    unsigned long nr = N.row();
    unsigned long nc = N.col();
    if(mr != nr) return 0;
    if(mc != nc) return 0;
    for(unsigned long i = 0; i < mr; i++) {
        for(unsigned long j = 0; j < mc; j++) {
            if(M[i][j] != N[i][j]) return 0;
        }
    }
    return 1;
}

//===== 行列計算エラー =====//
// 【発生要因】
// ・サイズが異なる行列同士の加算・減算を行った
// ・列数と行数が一致しない2つの行列をこの順番で乗算しようとした
template < class Type >
class MOperationErr{
  const Matrix<Type>* identA;
  const Matrix<Type>* identB;
public:
  MOperationErr(const Matrix<Type>* a, const Matrix<Type>* b) : identA(a), identB(b) {}
};


#endif
