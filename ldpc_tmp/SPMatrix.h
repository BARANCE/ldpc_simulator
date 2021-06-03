/*
    alistとspmat形式の疎行列を扱うクラス
*/

#if !defined(___Class_SPMatrix)
#define ___Class_SPMatrix

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include "Matrix.h"
#include "BinaryFiniteField.h"
#include "SetCalc.h"
using namespace std;

const int ps = 1;

class SPMatrix{
public :
    int rowSize;    // 行数
    int colSize;    // 列数
    int maxWr;    // 行重み
    int maxWc;    // 列重み
    vector<int> wr;    // 各行の行重みを記録したベクトル
    vector<int> wc;    // 各行の列重みを記録したベクトル
    vector<vector<int> >* matCol;
    vector<vector<int> >* matRow;
    string fn;    // ファイル名


    static const int WDYM = 0;    // spmat利用時のフラグ
    //==== コンストラクタ ====//
    SPMatrix();
    explicit SPMatrix(string fileName);    // ファイルから読み込み
    SPMatrix(string fileName, const int wdym);    // SPMat形式を読み込む
    SPMatrix(int r, int c);
    SPMatrix(const Matrix<BinaryFiniteField>& M);    // 変換
    SPMatrix(const SPMatrix& M);    // コピー
    ~SPMatrix(){
        delete matCol;
        delete matRow;
    }    // デストラクタ

    //==== 取得メソッド ====//
    Matrix<BinaryFiniteField> getMatrix() const;
    string getFileName() const {return fn;}
    int row() const {return rowSize;}
    int col() const {return colSize;}
    int getMaxWr() const {return maxWr;}
    int getMaxWc() const {return maxWc;}
    vector<int> getWrVector() const {return wr;}
    vector<int> getWcVector() const {return wc;}
    vector<vector<int> >* getSPMatrixCol() const {return matCol;}
    vector<vector<int> >* getSPMatrixColCopy() const;    // コピーを受け取るバージョン
    vector<vector<int> >* getSPMatrixRow() const {return matRow;}
    vector<vector<int> >* getSPMatrixRowCopy() const;
    int isZero() const;    // 行列が零行列かどうか
    void printProgress(int n) const;

    //==== 行列を切り出す ====//
    SPMatrix getRow(int index) const;
    SPMatrix getCol(int index) const;

    //==== 行列を変換する ====//
    SPMatrix transposed() const;    // 転置行列を得る

    //==== 行データ・列データの変換 ====//
    void replaceMatRowFromMatCol();    // matColからmatRowを生成
    void replaceMatColFromMatRow();    // matRowからmatColを生成

    //==== 設定メソッド ====//
    void readAlist();    // ファイルから読み込んでプロパティにセット
    void readAlist(string fileName);
    void readSpmat();    // ファイルから読み込んでプロパティにセット
    void readSpmat(string fileName);
    void writeAlist();
    void writeAlist(string fileName);
    void writeSpmat();    // プロパティを元にファイルに書き込み
    void writeSpmat(string fileName);
    void setFileName(string fileName) {fn = fileName;}
    void setMatrix(const Matrix<BinaryFiniteField>& M);
    void setMatrix(const SPMatrix& M);

    void ___setMaxWr(int n) {maxWr = n;}
    void ___setMaxWc(int n) {maxWc = n;}

    BinaryFiniteField getValue(int x, int y) const;    // x行y列の値を取得する。
    void setValue(int x, int y, const BinaryFiniteField value);    // x行y列にvalueを配置する。

    void removeCycles(); // 4サイクルループの除去

    //==== 内部計算用メソッド ====//
    int vecMatchCounter(const vector<int>& v1, const vector<int>& v2) const;    // v1とv2で一致している要素の数を数える
    int checkValueInVec(const vector<int>& v, const int& value) const;    // vにvalueが含まれているかを調べる
    vector<int> vecExclusiveOr(const vector<int>& v1, const vector<int>& v2) const;    // v1とv2の片方だけに含まれる要素を取り出す
    vector<int> removeValueByVec(const vector<int>& v, const int& value) const;    // vからvalueを取り除く

    //==== 演算子の二重定義 ====//
    operator Matrix<BinaryFiniteField>() const;
    SPMatrix& operator=(const Matrix<BinaryFiniteField>& M);
    SPMatrix& operator=(const SPMatrix& M);
    friend SPMatrix operator+(const SPMatrix& M, const SPMatrix& N);    // 加算
    friend SPMatrix operator-(const SPMatrix& M, const SPMatrix& N);
    friend SPMatrix operator*(const SPMatrix& M, const SPMatrix& N);    // 乗算
    friend SPMatrix operator*(const SPMatrix& M, const BinaryFiniteField& N);
    friend SPMatrix operator*(const BinaryFiniteField& M, const SPMatrix& N);
    friend int operator==(const SPMatrix& M, const SPMatrix& N);    // 比較

    static const int ZERO = 0;    // 零行列
    static const int E = 1;    // 単位行列
    static const int ONE = 2;    // 全一行列
    friend int operator==(const SPMatrix& M, const int& identifier);    // 零行列・単位行列・全一行列かどうかをチェック
    friend int operator!=(const SPMatrix& M, const int& identifier);

    //vector<int>& operator[](const int& i) {return (*matRow)[i];};


    //==== エラークラス ====//
    class UndefErr{    // 未定義値エラー
        SPMatrix* p;
    public:
        UndefErr(SPMatrix* ident) : p(ident) {}
    };

    class IOExceptions{    // ファイルオープンエラー
        SPMatrix* p;
    public:
        IOExceptions(SPMatrix* ident) : p(ident) {}
    };

    class EOFExceptions{    // 予定していないところでEOFが出現した
    };

    class FileFormatErr {    // ファイルのフォーマット(spmat形式)が正しくない
    };

    //==== 文字列変換 ====//
    stringstream& getStringstream() const {    // 他関数で使用する文字列ストリーム取得用関数
        stringstream ss;
        ss << colSize << " " << rowSize << "\n";
        ss << maxWc << " " << maxWr << "\n";
        for(int i = 0; i < colSize; i++) {
            //cout << wc[i] << " ";
            ss << wc[i];
            if(i < colSize - 1) ss << " ";
        }
        //cout << "\n";
        ss << "\n";
        for(int i = 0; i < rowSize; i++) {
            //cout << wr[i] << " ";
            ss << wr[i];
            if(i < rowSize - 1) ss << " ";
        }
        ss << "\n";
        //ss << "\n";
        for(int i = 0; i < colSize; i++) {
            for(int j = 0; j < wc[i]; j++) {
                ss << (*matCol)[i][j];
                if(j < wc[i] - 1) ss << " ";
            }
            ss << "\n";
        }
        for(int i = 0; i < rowSize; i++) {
            for(int j = 0; j < wr[i]; j++) {
                ss << (*matRow)[i][j];
                if(j < wr[i] - 1) ss << " ";
            }
            ss << "\n";
        }
        return ss;
    }
    string toString() const {return getStringstream().str();}
    void print() const {cout << toString();}
    friend ostream& operator<<(ostream& ss, const SPMatrix& v) {return ss << v.toString();}
};

#endif
