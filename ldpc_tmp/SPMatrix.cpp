/*
    spmat形式の疎行列を扱うクラス(実装)
*/

#include "SPMatrix.h"
using namespace std;

SPMatrix::SPMatrix() : rowSize(0), colSize(0), maxWr(0), maxWc(0), wr(0), wc(0), fn(""), matCol(new vector<vector<int> >), matRow(new vector<vector<int> >) {
}
SPMatrix::SPMatrix(int r, int c) : rowSize(r), colSize(c), maxWr(0), maxWc(0), wr(r,0), wc(c,0), fn(""), matCol(new vector<vector<int> >), matRow(new vector<vector<int> >) {
    vector<vector<int> >* matColTmp = new vector<vector<int> >(colSize);
    vector<vector<int> >* matRowTmp = new vector<vector<int> >(rowSize);
    delete matCol;
    delete matRow;
    matCol = matColTmp;
    matRow = matRowTmp;
}
SPMatrix::SPMatrix(string fileName) : rowSize(0), colSize(0), maxWr(0), maxWc(0), wr(0), wc(0), fn(fileName), matCol(new vector<vector<int> >), matRow(new vector<vector<int> >) {
    readAlist();
}
SPMatrix::SPMatrix(string fileName, const int wdym) : rowSize(0), colSize(0), maxWr(0), maxWc(0), wr(0), wc(0), fn(fileName), matCol(new vector<vector<int> >), matRow(new vector<vector<int> >) {
    readSpmat();
}
SPMatrix::SPMatrix(const Matrix<BinaryFiniteField>& M) : rowSize(0), colSize(0), maxWr(0), maxWc(0), wr(0), wc(0), fn(""), matCol(new vector<vector<int> >), matRow(new vector<vector<int> >) {
    setMatrix(M);
}
SPMatrix::SPMatrix(const SPMatrix& M) : rowSize(0), colSize(0), maxWr(0), maxWc(0), wr(0), wc(0), fn(""), matCol(new vector<vector<int> >), matRow(new vector<vector<int> >) {
    setMatrix(M);
}

//==== 取得メソッド ====//
Matrix<BinaryFiniteField> SPMatrix::getMatrix() const {
    Matrix<BinaryFiniteField> retMat(rowSize, colSize);
    for(int i = 0; i < colSize; i++) {
        for(int j = 0; j < wc[i]; j++) {
            //cout << (*mat)[i][j] << " ";
            retMat[(*matCol)[i][j] - 1][i] = 1;
        }
        //cout << "\n";
    }
    return retMat;
}
void SPMatrix::printProgress(int n) const {
    for(int i = 0; i < n; i++) {
        cout << "■";
    }
}
// コピーを受け取るバージョン
vector<vector<int> >* SPMatrix::getSPMatrixRowCopy() const {
    vector<vector<int> >* retMatRow = new vector<vector<int> >(rowSize);
    for(int i = 0; i < rowSize; i++) {
        for(int j = 0; j < wr[i]; j++) {
            int tmp = (*matRow)[i][j];
            (*retMatRow)[i].push_back(tmp);
        }
    }
    return retMatRow;
}
vector<vector<int> >* SPMatrix::getSPMatrixColCopy() const {
    vector<vector<int> >* retMatCol = new vector<vector<int> >(colSize);
    for(int i = 0; i < colSize; i++) {
        for(int j = 0; j < wc[i]; j++) {
            int tmp = (*matCol)[i][j];
            (*retMatCol)[i].push_back(tmp);
        }
    }
    return retMatCol;
}
// 行列が零行列ならば1を返す
int SPMatrix::isZero() const {
    if(maxWr == 0 && maxWc == 0) return 1;
    else return 0;
}

//==== 行列を切り出す ====//
SPMatrix SPMatrix::getRow(int index) const {
    Matrix<BinaryFiniteField> matTmp(1, colSize);
    for(int i = 0; i < wr[index]; i++) {
        matTmp[0][(*matRow)[index][i] - 1] = 1;
    }
    SPMatrix retMat = matTmp;
    return retMat;
}
SPMatrix SPMatrix::getCol(int index) const {
    Matrix<BinaryFiniteField> matTmp(rowSize, 1);
    for(int i = 0; i < rowSize; i++) {
        for(int j = 0; j < wr[i]; j++) {
            if(((*matRow)[i][j] - 1) == index) matTmp[i][0] = 1;
        }
    }
    SPMatrix retMat = matTmp;
    return retMat;
}

//==== 行列を変換する ====//
SPMatrix SPMatrix::transposed() const {
    SPMatrix retMat;
    retMat.fn = fn;
    retMat.rowSize = colSize;
    retMat.colSize = rowSize;
    retMat.maxWr = maxWc;
    retMat.maxWc = maxWr;
    retMat.wr = wc;
    retMat.wc = wr;
    delete retMat.matRow;
    retMat.matRow = new vector<vector<int> >(retMat.rowSize);
    for(int i = 0; i < retMat.rowSize; i++) {
        for(int j = 0; j < retMat.wr[i]; j++) {
            (*retMat.matRow)[i].push_back((*matCol)[i][j]);
        }
    }
    delete retMat.matCol;
    retMat.matCol = new vector<vector<int> >(retMat.colSize);
    for(int i = 0; i < retMat.colSize; i++) {
        for(int j = 0; j < retMat.wc[i]; j++) {
            (*retMat.matCol)[i].push_back((*matRow)[i][j]);
        }
    }
    return retMat;
}
//==== 行データ・列データの変換 ====//
void SPMatrix::replaceMatRowFromMatCol(){    // matColからmatRowを生成
    delete matRow;
    matRow = new vector<vector<int> >(rowSize, vector<int>(0));
    for(int i = 0; i < rowSize; i++) {
        (*matRow)[i].reserve(wr[i]);    // 領域を確保
    }
    for(int i = 0; i < colSize; i++) {
        for(int j = 0; j < wc[i]; j++) {
            int pos = (*matCol)[i][j] - 1;
            (*matRow)[pos].push_back(i + 1);
        }
    }
    for(int i = 0; i < rowSize; i++) {
        std::sort((*matRow)[i].begin(), (*matRow)[i].end());
    }
}
void SPMatrix::replaceMatColFromMatRow(){    // matRowからmatColを生成
    delete matCol;
    matCol = new vector<vector<int> >(colSize, vector<int>(0));
    for(int i = 0; i < colSize; i++) {
        (*matCol)[i].reserve(wc[i]);    // 領域を確保
    }
    for(int i = 0; i < rowSize; i++) {
        for(int j = 0; j < wr[i]; j++) {
            int pos = (*matRow)[i][j] - 1;
            (*matCol)[pos].push_back(i + 1);
        }
    }
    for(int i = 0; i < colSize; i++) {
        std::sort((*matCol)[i].begin(), (*matCol)[i].end());
    }
}

//==== 設定メソッド ====//
// ファイルから読み込んでプロパティにセット
void SPMatrix::readAlist() {
    if (fn.empty()) throw UndefErr(this);    // ファイル名が空白
    ifstream fs(fn.c_str(), fstream::in);
    if (fs.fail()) throw IOExceptions(this);    // ファイルが存在しない

    if(fs.peek() == '\n' || fs.eof()) {throw FileFormatErr();}

    fs >> colSize >> rowSize;
    fs >> maxWc >> maxWr;
    vector<int> wcVec(colSize);
    for(int i = 0; i < colSize; i++) {
        fs >> wcVec[i];
    }

    wc = wcVec;
    vector<int> wrVec(rowSize);
    for(int i = 0; i < rowSize; i++) {
        fs >> wrVec[i];
    }
    wr = wrVec;


    vector<vector<int> >* matTmp = new vector<vector<int> >(colSize);
    for(int i = 0; i < colSize; i++) {
        for(int j = 0; j < wc[i]; j++) {
            int tmp;
            fs >> tmp;
            //cout << tmp << " ";
            (*matTmp)[i].push_back(tmp);
        }
        //cout << "\n";
    }

    vector<vector<int> >* matRowTmp = new vector<vector<int> >(rowSize);
    for(int i = 0; i < rowSize; i++) {
        for(int j = 0; j < wr[i]; j++) {
            int tmp;
            fs >> tmp;
            (*matRowTmp)[i].push_back(tmp);
        }
    }

    fs.close();

    delete matCol;
    delete matRow;
    matCol = matTmp;
    matRow = matRowTmp;
}
void SPMatrix::readAlist(string fileName) {
    fn = fileName;
    readAlist();
}
void SPMatrix::readSpmat() {
    if (fn.empty()) throw UndefErr(this);    // ファイル名が空白
    ifstream fs(fn.c_str(), fstream::in);
    if (fs.fail()) throw IOExceptions(this);    // ファイルが存在しない

    if(fs.peek() == '\n' || fs.eof()) {throw FileFormatErr();}
    fs >> colSize >> rowSize;
    fs >> maxWr >> maxWc;
    vector<int> wrVec(rowSize);
    for(int i = 0; i < rowSize; i++) {
        fs >> wrVec[i];
    }
    wr = wrVec;
    vector<int> wcVec(colSize);
    for(int i = 0; i < colSize; i++) {
        fs >> wcVec[i];
    }
    wc = wcVec;

    vector<vector<int> >* matRowTmp = new vector<vector<int> >(rowSize);
    vector<vector<int> >* matColTmp = new vector<vector<int> >(colSize);
    for(int i = 0; i < rowSize; i++) {
        (*matRowTmp)[i].reserve(wr[i]);
        for(int j = 0; j < wr[i]; j++) {
            int tmp;
            fs >> tmp;
            (*matRowTmp)[i].push_back(tmp);
            (*matColTmp)[tmp - 1].push_back(i + 1);
        }
    }
    fs.close();

    delete matRow;
    delete matCol;
    matRow = matRowTmp;
    matCol = matColTmp;
}
void SPMatrix::readSpmat(string fileName) {
    fn = fileName;
    readSpmat();
}

void SPMatrix::setMatrix(const Matrix<BinaryFiniteField>& M) {
    rowSize = M.row();
    colSize = M.col();
    int wrCounter = 0;
    int wcCounter = 0;
    int mwr = 0; // 最大行重み・列重みは測定する
    int mwc = 0;
    vector<int> wrTmp(0);
    wrTmp.reserve(rowSize);
    vector<int> wcTmp(0);
    wcTmp.reserve(colSize);

    vector<vector<int> >* matRowTmp = new vector<vector<int> >(rowSize);
    vector<vector<int> >* matColTmp = new vector<vector<int> >(colSize);
    for(int i = 0; i < rowSize; i++) {
        wrCounter = 0;
        for(int j = 0; j < colSize; j++) {
            if((int)M[i][j] > 0) {
                (*matRowTmp)[i].push_back(j + 1);
                wrCounter++;
            }
        }
        wrTmp.push_back(wrCounter);
        if(wrCounter > mwr) mwr = wrCounter;
    }

    for(int j = 0; j < colSize; j++) {
        wcCounter = 0;
        for(int i = 0; i < rowSize; i++) {
            if((int)M[i][j] > 0) {
                (*matColTmp)[j].push_back(i + 1);
                wcCounter++;
            }
        }
        wcTmp.push_back(wcCounter);
        if(wcCounter > mwc) mwc = wcCounter;
    }

    maxWr = mwr;
    maxWc = mwc;
    wr = wrTmp;
    wc = wcTmp;

    delete matCol;
    delete matRow;
    matCol = matColTmp;
    matRow = matRowTmp;

}
void SPMatrix::setMatrix(const SPMatrix& M) {
    rowSize = M.row();
    colSize = M.col();
    maxWr = M.getMaxWr();
    maxWc = M.getMaxWc();

    vector<int> wrTmp(rowSize);
    vector<int> wrCopy = M.getWrVector();
    vector<int> wcTmp(colSize);
    vector<int> wcCopy = M.getWcVector();
    for(int i = 0; i < rowSize; i++) {
        wrTmp[i] = wrCopy[i];
    }
    for(int j = 0; j < colSize; j++) {
        wcTmp[j] = wcCopy[j];
    }
    wr = wrTmp;
    wc = wcTmp;

    vector<vector<int> >* matRowTmp = new vector<vector<int> >(rowSize);
    vector<vector<int> >* matColTmp = new vector<vector<int> >(colSize);
    vector<vector<int> >* matRowCopy = M.getSPMatrixRow();
    vector<vector<int> >* matColCopy = M.getSPMatrixCol();
    for(int i = 0; i < rowSize; i++) {
        (*matRowTmp)[i].reserve(wrTmp[i]);
        for(int j = 0; j < wrTmp[i]; j++) {
            (*matRowTmp)[i].push_back((*matRowCopy)[i][j]);
        }
    }
    for(int j = 0; j < colSize; j++) {
        (*matColTmp)[j].reserve(wcTmp[j]);
        for(int i = 0; i < wcTmp[j]; i++) {
            (*matColTmp)[j].push_back((*matColCopy)[j][i]);
        }
    }

    delete matRow;
    delete matCol;
    matRow = matRowTmp;
    matCol = matColTmp;

    fn = M.getFileName();
}

// プロパティを元にファイルに書き込み
void SPMatrix::writeAlist() {
    if (fn.empty()) throw UndefErr(this);    // ファイル名が空白
    ofstream fs(fn.c_str(), ios::out);
    if (fs.fail()) throw IOExceptions(this);    // ファイルに書き込めない(使用中orアクセス権なし)

    fs << colSize << " " << rowSize << "\n";
    fs << maxWc << " " << maxWr << "\n";
    for(int i = 0; i < colSize; i++) {
        fs << wc[i];
        if(i < colSize - 1) fs << " ";
    }
    fs << "\n";
    for(int i = 0; i < rowSize; i++) {
        fs << wr[i];
        if(i < rowSize - 1) fs << " ";
    }
    fs << "\n";
    for(int i = 0; i < colSize; i++) {
        for(int j = 0; j < wc[i]; j++) {
            fs << (*matCol)[i][j];
            if(j < wc[i] - 1) fs << " ";
        }
        fs << "\n";
    }
    for(int i = 0; i < rowSize; i++) {
        for(int j = 0; j < wr[i]; j++) {
            fs << (*matRow)[i][j];
            if(j < wr[i] - 1) fs << " ";
        }
        fs << "\n";
    }
    fs.close();
}
void SPMatrix::writeAlist(string fileName) {
    fn = fileName;
    writeAlist();
}
void SPMatrix::writeSpmat() {
    if (fn.empty()) throw UndefErr(this);    // ファイル名が空白
    ofstream fs(fn.c_str(), ios::out);
    if (fs.fail()) throw IOExceptions(this);    // ファイルに書き込めない(使用中orアクセス権なし)

    fs << colSize << " " << rowSize << "\n";
    fs << maxWr << " " << maxWc << "\n";
    for(int i = 0; i < rowSize; i++) {
        fs << wr[i];
        if(i < rowSize - 1) fs << " ";
    }
    fs << "\n";
    for(int i = 0; i < colSize; i++) {
        fs << wc[i];
        if(i < colSize - 1) fs << " ";
    }
    fs << "\n";
    for(int i = 0; i < rowSize; i++) {
        for(int j = 0; j < wr[i]; j++) {
            fs << (*matRow)[i][j];
            if(j < wr[i] - 1) fs << " ";
        }
        fs << "\n";
    }
    fs.close();
}
void SPMatrix::writeSpmat(string fileName) {
    fn = fileName;
    writeSpmat();
}

// x行y列の値を取得する。
BinaryFiniteField SPMatrix::getValue(int x, int y) const {
    if(x < 0 || x >= rowSize) throw UndefErr(NULL);
    if(y < 0 || y >= colSize) throw UndefErr(NULL);
    for(int i = 0; i < wr[x]; i++) {
        if((*matRow)[x][i] == y + 1) {
            return 1;
        }
    }
    return 0;
}
// x行y列にvalueを配置する。
void SPMatrix::setValue(int x, int y, const BinaryFiniteField value){
    if(x < 0 || x >= rowSize) throw UndefErr(this);
    if(y < 0 || y >= colSize) throw UndefErr(this);
    if(value == 1) {
        int check = 0;
        for(int i = 0; i < wr[x]; i++) {
            if((*matRow)[x][i] == y + 1) {
                check = 1;
                break;
            }
        }
        if(check == 0) {
            wr[x]++;
            wc[y]++;
            if(wr[x] > maxWr) maxWr = wr[x];
            if(wc[y] > maxWc) maxWc = wc[y];
            (*matRow)[x].push_back(y + 1);
            (*matCol)[y].push_back(x + 1);
            std::sort((*matRow)[x].begin(), (*matRow)[x].end());
            std::sort((*matCol)[y].begin(), (*matCol)[y].end());
        }
    } else if(value == 0) {    // 激遅なので注意！！
        int pos = -1;
        for(int i = 0; i < wr[x]; i++) {
            if((*matRow)[x][i] == y + 1) {
                pos = i;
                break;
            }
        }
        if(pos >= 0) {
            vector<int> tmpMatRowX(wr[x] - 1);
            for(int i = 0; i < pos; i++) {
                tmpMatRowX[i] = (*matRow)[x][i];
            }
            for(int i = pos + 1; i < wr[x]; i++) {
                tmpMatRowX[i - 1] = (*matRow)[x][i];
            }
            (*matRow)[x] = tmpMatRowX;
            wr[x]--;
            wc[y]--;
            int mr = 0, mc = 0;
            for(int i = 0; i < rowSize; i++) {
                if(wr[i] > mr) mr = wr[i];
            }
            maxWr = mr;
            for(int i = 0; i < colSize; i++) {
                if(wc[i] > mc) mc = wc[i];
            }
            maxWc = mc;
            replaceMatColFromMatRow();
        }
    } else {
        throw UndefErr(this);
    }
}

// 4サイクルループの除去
void SPMatrix::removeCycles() {
    SetCalc<int> sc;
    for(int i = 0; i < colSize; i++) {
        for(int j = i + 1; j < colSize; j++) {
            /*int mis = (*matCol)[i].size();    // debug
            cout << i << " [";    //debug
            for(int a = 0; a < mis; a++) {
                cout << (*matCol)[i][a] << " ";
            }
            cout << "]\n";
            int mjs = (*matCol)[j].size();    // debug
            cout << j << " [";    //debug
            for(int a = 0; a < mjs; a++) {
                cout << (*matCol)[j][a] << " ";
            }
            cout << "] → ";*/

            (*matCol)[j] = sc.spor(&(*matCol)[i], &(*matCol)[j]);
            wc[j] = (*matCol)[j].size();

            /*mjs = (*matCol)[j].size();
            cout << j << " [";    //debug
            for(int a = 0; a < mjs; a++) {
                cout << (*matCol)[j][a] << " ";
            }
            cout << "]\n";
            //break;*/
        }
        //break;    //debug
    }

    vector<vector<int> >* matRowTmp = new vector<vector<int> >(rowSize);
    vector<int> wrTmp(rowSize, 0);
    for(int i = 0; i < colSize; i++) {
        for(int j = 0; j < wc[i]; j++) {
            (*matRowTmp)[(*matCol)[i][j] - 1].push_back(i + 1);
            wrTmp[(*matCol)[i][j] - 1]++;
        }
    }
    wr = wrTmp;
    delete matRow;
    matRow = matRowTmp;

    int max = 0;
    for(int i = 0; i < rowSize; i++) {
        if(wr[i] > max) max = wr[i];
    }
    maxWr = max;
    max = 0;
    for(int i = 0; i < colSize; i++) {
        if(wc[i] > max) max = wc[i];
    }
    maxWc = max;
}

//==== 内部計算用メソッド ====//
int SPMatrix::vecMatchCounter(const vector<int>& v1, const vector<int>& v2) const {
    int counter = 0;
    int size = v1.size();
    for(int i = 0; i < size; i++) {
        if(checkValueInVec(v2, v1[i])) counter++;
    }
    return counter;
}
int SPMatrix::checkValueInVec(const vector<int>& v, const int& value) const {
    int size = v.size();
    for(int i = 0; i < size; i++) {
        if(v[i] == value) return 1;
    }
    return 0;
}
vector<int> SPMatrix::vecExclusiveOr(const vector<int>& v1, const vector<int>& v2) const {
    vector<int> retVec = v1;
    int size = v2.size();
    for(int i = 0; i < size; i++) {
        retVec = removeValueByVec(retVec, v2[i]);
    }
    return retVec;
}
vector<int> SPMatrix::removeValueByVec(const vector<int>& v, const int& value) const {
    vector<int> retVec(0);
    int size = v.size();
    for(int i = 0; i < size; i++) {
        if(v[i] != value) retVec.push_back(v[i]);
    }
    return retVec;
}

//==== 演算子の多重定義 ====//
SPMatrix::operator Matrix<BinaryFiniteField>() const {
    return getMatrix();
}
SPMatrix& SPMatrix::operator=(const Matrix<BinaryFiniteField>& M) {
    setMatrix(M);
    return (*this);
}
SPMatrix& SPMatrix::operator=(const SPMatrix& M) {
    setMatrix(M);
    return (*this);
}
SPMatrix operator+(const SPMatrix& M, const SPMatrix& N) {
    int mr = M.row();
    int mc = M.col();
    SPMatrix retMat;
    int maxWr = 0;
    int maxWc = 0;
    vector<int> wrTmp(mr);
    vector<int> wcTmp(mc);
    vector<vector<int> >* matRowTmp = new vector<vector<int> >(mr);
    vector<vector<int> >* matColTmp = new vector<vector<int> >(mc);
    for(int i = 0; i < mr; i++) {
        (*matRowTmp)[i] = retMat.vecExclusiveOr((*M.matRow)[i], (*N.matRow)[i]);
        wrTmp[i] = (*matRowTmp)[i].size();
        if(wrTmp[i] > maxWr) maxWr = wrTmp[i];
    }
    for(int i = 0; i < mc; i++) {
        (*matColTmp)[i] = retMat.vecExclusiveOr((*M.matCol)[i], (*N.matCol)[i]);
        wcTmp[i] = (*matColTmp)[i].size();
        if(wcTmp[i] > maxWc) maxWc = wcTmp[i];
    }

    retMat.rowSize = mr;
    retMat.colSize = mc;
    retMat.maxWr = maxWr;
    retMat.maxWc = maxWc;
    retMat.wr = wrTmp;
    retMat.wc = wcTmp;
    delete retMat.matRow;
    retMat.matRow = matRowTmp;
    delete retMat.matCol;
    retMat.matCol = matColTmp;
    retMat.setFileName(M.getFileName());
    return retMat;
}
SPMatrix operator-(const SPMatrix& M, const SPMatrix& N) {
    return M + N;
}
SPMatrix operator*(const SPMatrix& M, const SPMatrix& N) {

    // 糞
    int mr = M.row();
    int mc = M.col();
    int nr = N.row();
    int nc = N.col();
    vector<int> mwr = M.getWrVector();
    vector<int> nwc = N.getWcVector();
    vector<vector<int> >* matRowM = M.getSPMatrixRow();
    vector<vector<int> >* matColN = N.getSPMatrixCol();
    vector<vector<int> >* matRowTmp = new vector<vector<int> >(mr);
    vector<vector<int> >* matColTmp = new vector<vector<int> >(nc);

    // 書き込み用情報
    SPMatrix retMat;
    retMat.rowSize = mr;
    retMat.colSize = nc;
    vector<int> wr(mr, 0);
    vector<int> wc(nc, 0);
    int maxWr = 0;
    int maxWc = 0;

    //cout << "mr = " << mr << ", nc = " << nc << "\n";

    for(int i = 0; i < mr; i++) {
        for(int j = 0; j < nc; j++) {
            int match = retMat.vecMatchCounter((*matRowM)[i], (*matColN)[j]);
            if(match % 2 != 0) {
                (*matRowTmp)[i].push_back(j + 1);
                wr[i]++;
                (*matColTmp)[j].push_back(i + 1);
                wc[j]++;
            }
            //cout << "[" << i << "][" << j << "]match = " << match << "\n";
            //cout << "wr[" << i << "] = " << wr[i] << ", wc[" << j << "] = " << wc[j] << "\n";
            if(wc[j] > maxWc) maxWc = wc[j];
        }
        if(wr[i] > maxWr) maxWr = wr[i];
    }
    //cout << "maxWr = " << maxWr << ", maxWc = " << maxWc << "\n";
    retMat.maxWr = maxWr;
    retMat.maxWc = maxWc;
    retMat.wr = wr;
    retMat.wc = wc;
    delete retMat.matCol;
    retMat.matCol = matColTmp;
    delete retMat.matRow;
    retMat.matRow = matRowTmp;

    return retMat;
}
SPMatrix operator*(const SPMatrix& M, const BinaryFiniteField& N) {
    BinaryFiniteField tmp = N;
    SPMatrix tmp2 = M;
    if((int)tmp == 1) {
    } else {
        tmp2.maxWr = 0;
        tmp2.maxWc = 0;
        delete tmp2.matRow;
        tmp2.matRow = new vector<vector<int> >(tmp2.rowSize);
        delete tmp2.matCol;
        tmp2.matCol = new vector<vector<int> >(tmp2.colSize);
        for(int i = 0; i < tmp2.rowSize; i++) {
            tmp2.wr[i] = 0;
        }
        for(int i = 0; i < tmp2.colSize; i++) {
            tmp2.wc[i] = 0;
        }
    }
    return tmp2;
}
SPMatrix operator*(const BinaryFiniteField& M, const SPMatrix& N) {
    return N * M;
}
int operator==(const SPMatrix& M, const int& identifier) {
    if(identifier == SPMatrix::ZERO) return M.isZero();
    else return 0;
}
int operator!=(const SPMatrix& M, const int& identifier) {
    if(M == identifier) return 0;
    else return 1;
}
