/*
    �s��������N���X�e���v���[�gMatrix
    �y�@�\�z
    �E�e���v���[�g�Ȃ̂ŁA�C�ӂ̉Z�N���X���s��̐����Ƃ��Ĉ�����
    �E�t�@�C������̍s��I�u�W�F�N�g�쐬
    �E�P�ʍs��̐���
    �E�s�񓯎m�̉�����Z
    �E�萔�{
    �E�s��̐؂�o����A�s�񓯎m�̐ڍ�
    �E�]�u�s��
    �E��u���A�s�u��
*/


#if !defined(___Class_Matrix)
#define ___Class_Matrix

//#include "MatrixFileOperator.h"
//#include "IOOError.h"    // �G���[����
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include "RandomValue.h"
#include "BinaryFiniteField.h"

using namespace std;

template < class Type > class MOperationErr;


//==== �s��N���X�e���v���[�g =====//
template <class Type> class Matrix {
    unsigned long rowSize;    // �s��̍s��
    unsigned long colSize;    // �s��̗�
    /*int maxRowSize;    // �m�ۍς݂̍s��
    int maxColSize;    // �m�ۍς݂̗�
    Type* initializer;    // �s��̏������q
    bool initExist;    // initializer�����݂��邩�ǂ���*/

    vector<vector<Type> >* mat;    // �s��{��
public :
    //===== �R���X�g���N�^ =====//
    Matrix() : rowSize(0), colSize(0){    // �f�t�H���g�R���X�g���N�^
        regionSet(rowSize, colSize);
    }
    Matrix(unsigned long r, unsigned long c) : rowSize(r), colSize(c) {    // r�~c�s��𐶐�
        regionSet(rowSize, colSize);
    }
    Matrix(const Matrix& x) : rowSize(x.row()), colSize(x.col()) {    // �R�s�[�R���X�g���N�^
        regionSet(rowSize, colSize);
        //cout << "rowSize = " << rowSize << ", colSize = " << colSize << "\n";
        matCopy(x);
    }
    explicit Matrix(const string fileName) {    // �t�@�C������
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
        if (ifs.eof()) throw FileFormatErr<Type>(this, "FileFormatErr:�t�@�C���̍s���܂��͗񐔂��ǂݍ��߂܂���B\n");
        regionSet(rowSize, colSize);
        for(int i = 0; i < rowSize; i++) {
            for(int j = 0; j < colSize; j++) {
                if (ifs.eof()) throw FileFormatErr<Type>(this, "FileFormatErr:�w�肳�ꂽ�s���E�񐔂ɒB����O�Ƀt�@�C���̏I�[�ɒB���܂����B\n");
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

    //===== �ϊ��R���X�g���N�^ =====//
    Matrix(const vector<Type> x) : rowSize(1), colSize(x.size()) {
        regionSet(rowSize, colSize);
        vecCopy(x);
    }

    //===== �f�X�g���N�^ =====//
    ~Matrix() {
        delete mat;
        mat = NULL;
    }

    //===== �l�擾�p���\�b�h =====//
    unsigned long row() const {return rowSize;};    // �s����Ԃ�
    unsigned long col() const {return colSize;};    // �񐔂�Ԃ�
    stringstream getStringstream() const {    // ���֐��Ŏg�p���镶����X�g���[���擾�p�֐�
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

    //===== �s���؂�o�� =====//
    // index�s�𔲂��o��
    Matrix<Type> getRow(unsigned long index) const {
        Matrix<Type> valMat(1, colSize);
        for(unsigned long i = 0; i < colSize; i++) {
            valMat[0][i] = (*mat)[index][i];
        }
        //cout << "rowSize = " << valMat.row() << ", colSize = " << valMat.col() << "\n";
        return valMat;
    }
    // index��𔲂��o��
    Matrix<Type> getCol(unsigned long index) const {
        Matrix<Type> valMat(rowSize, 1);
        for(unsigned long i = 0; i < rowSize; i++) {
            valMat[i][0] = (*mat)[i][index];
        }
        return valMat;
    }
    // �����s��𓾂�
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

    //==== �s�܂��͗�̍폜 ====//
    Matrix<Type> deleteRow(unsigned long index) const {    // index�s���폜����
        Matrix<Type> valMat(rowSize - 1, colSize);
        for(unsigned long i = 0; i < index; i++) {    // index�s�̑O�܂�
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

    Matrix<Type> deleteCol(unsigned long index) const {    // index����폜����
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
    Matrix<Type> deleteRowCol(unsigned long row, unsigned long col) const {    // row�s��col����폜
        Matrix<Type> valMat(rowSize - 1, colSize - 1);
        if(row >= rowSize) cout << "�Y���G���[\n";
        if(col >= colSize) cout << "�Y���G���[\n";
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

    //===== �s���ϊ����� =====//
    // �]�u�s��𓾂�
    Matrix<Type> transposed() {
        Matrix<Type> valMat(colSize, rowSize);
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                valMat[j][i] = (*mat)[i][j];
            }
        }
        return valMat;
    }

    //==== �s��̐ڍ� ====//
    // ������ɍs���ڍ�����
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

    // �������ɍs���ڍ�����
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
    // �������ɍs���ڍ�����
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
    // �E�����ɍs���ڍ�����
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

    //==== �s�܂��͗�̒u��(����) ====//
    // col1���col2����u��
    void colSubstitution(const unsigned long col1, const unsigned long col2) {
        vector<Type> tmpVec(rowSize);
        for(unsigned long i = 0; i < rowSize; i++) {
            tmpVec[i] = (*mat)[i][col1];
            (*mat)[i][col1] = (*mat)[i][col2];
            (*mat)[i][col2] = tmpVec[i];
        }
    }
    // row1�s��row2�s���s�u��
    void rowSubstitution(const unsigned long row1, const unsigned long row2) {
        vector<Type> tmpVec(colSize);
        for(unsigned long j = 0; j < colSize; j++) {
            tmpVec[j] = (*mat)[row1][j];
            (*mat)[row1][j] = (*mat)[row2][j];
            (*mat)[row2][j] = tmpVec[j];
        }
    }
    // �����_����u����1�񂾂��s��
    void randomColSubstitution() {
        RandomValue rv;
        unsigned long r1 = rv.getRand(colSize);
        unsigned long r2 = rv.getRand(colSize);
        while(r1 == r2) {r2 = rv.getRand(colSize);}
        //cout << "r1 = " << r1 << ", r2 = " << r2 << "\n";
        colSubstitution(r1, r2);
    }

    //==== �s��{�ό` ====//
    // A�s�ڂ�B�s�ڂɑ���
    void AddRowFromAToB(const unsigned long& A, const unsigned long& B) {
        for(unsigned long j = 0; j < colSize; j++) {
            (*mat)[B][j] += (*mat)[A][j];
        }
    }
    // A�s�ڂ�B�s�ڂ������
    void SubRowFromAToB(const unsigned long& A, const unsigned long& B) {
        for(unsigned long j = 0; j < colSize; j++) {
            (*mat)[B][j] -= (*mat)[A][j];
        }
    }

    //==== �s��̕W���`���擾(�����s���) ====//
    // ���s��{�ό`�݂̂ŕW���`�Ɏ���Ȃ��ꍇ�́A��u�����s��
    Matrix<Type> getRENF_H() {
        Matrix<Type> valMat = (*this);    // �����s��
        unsigned long hr = rowSize;
        unsigned long hc = colSize;
        bool deleteCheck = false;
        for(unsigned long s = 0; s < hr; s++) {
            unsigned long n = hc - hr;
            //cout << s << "�s" << s+n << "��ڂ̏����F";
            if((int)valMat[s][s + n] == 0) {
                int checker = 0;
                for(unsigned long i = s + 1; i < hr; i++) {
                    //cout << valMat[i][s + n] << "\n";
                    if((int)valMat[i][s + n] != 0) {
                        valMat.rowSubstitution(s, i);    // [s][s+n]�v�f��1�ƂȂ�悤�ɍs�u��
                        checker = 1;
                        break;
                    }
                }
                if(checker == 0) {    // ��u�����K�v
                    //cout << "[*info*] " << s+n << "��ڂɂ�1�̂���s������܂���B";
                    //cout << "\n";
                    for(unsigned long j = 0; j < n; j++) {
                        if((int)valMat[s][j] != 0) {
                            valMat.colSubstitution(j, s + n);    // [s][s+n]�v�f��1�ƂȂ�悤�ɗ�u��
                            (*this).colSubstitution(j, s + n);
                            //cout << "[*info*] " << j << "��ڂ�" << s + n << "��ڂ�u�����܂����B�����s��̗�u�����s���܂����B\n";
                            checker = 2;
                            break;
                        }
                    }
                    if(checker == 0) {    // ��u�����ł��Ȃ�
                        //cout << "[*info*] " << s << "�s�ڂ��S��ł��B�W���`�ɕό`�ł��܂���B\n";
                        //valMat[s][s+n] = 1;
                        //valMat.exportMatrix("tmp1.txt");
                        valMat = valMat.deleteRowCol(s, s+n);
                        hr--;
                        hc--;
                        n = hc - hr;
                        //valMat.exportMatrix("tmp2.txt");
                        deleteCheck = true;
                        //cout << "[*info*] " << s << "�s�A�����" << s+n << "��ڂ��폜���܂����B\n";
                    }
                } else {}//cout << "\n";
            }
            if(deleteCheck == false && (int)valMat[s][s + n] != 0) {
                // s+1�s�`rowSize�s��s+n��ڂ�1�������Ă���s��s�s�ڂ𑫂�
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

    //===== �I�u�W�F�N�g����p���\�b�h =====//
    void regionSet(const unsigned long r, const unsigned long c) {    // r�~c�s���ݒu(�����̂��̂͏���)
        try {
            mat = new vector<vector<Type> >(r, vector<Type>(c));
        } catch(std::bad_alloc ba) {    // �̈�m�ێ��s
            //throw MatrixBadAlloc(&ba, r, c);
            throw this;
        }
    }
    void matCopy(const Matrix<Type>& M) {    // M���R�s�[(�����̂��̂͏���)
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
    void reserve(unsigned long r, unsigned long c) {    // r�~c�s��̗̈���m��(�����̂��̂��ێ�)
        vector<vector<Type> >* tmpMat;    // �R�s�[�����p�̔z��
        Type* a;
        tmpMat = getMatrixRegion(r, c, a);
        for(unsigned long i = 0; i < rowSize || i < r; i++) {
            for(unsigned long j = 0; j < colSize || j < c; j++) {
                (*tmpMat)[i][j] = (*mat)[i][j];    // �s��R�s�[
            }
        }
        rowSize = r;
        colSize = c;
        delete mat;
        mat = NULL;
        mat = tmpMat;
    }
    // �P�ʍs��(�����̂��̂͏���)
    void setIdentity() {
        for(unsigned long i = 0; i < rowSize; i++) {
            for(unsigned long j = 0; j < colSize; j++) {
                if(i == j) (*mat)[i][j] = 1;
                else (*mat)[i][j] = 0;
            }
        }
    }

    // �t�@�C���ɕۑ�
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

    //===== ���Z�q�̓�d��`[1] =====//
    Matrix<Type>& operator=(const Matrix<Type>& M) {    // ������Z�q
        if(&M != this) {    // �������g�ł͖���
            rowSize = M.row();
            colSize = M.col();
            delete mat;
            mat = NULL;
            regionSet(rowSize, colSize);
            matCopy(M);
        }
        return *this;
    }
    vector<Type>& operator[](const int i) const {return (*mat)[i];}    // �Y�����Z�q
    vector<Type>& operator()(const int i) const {return (*mat)[i];}    // ��Ɠ���
    Type& operator()(const int i, const int j) const {return (*mat)[i][j];}    // �A�N�Z�X
};

//===== �N���X�Ŏg���֐� =====//
// rowSize�~colSize�̍s��̈���m�ۂ���
// �^Type�Ƀf�t�H���g�R���X�g���N�^�������ꍇ�͖���`
// rowSize : �s��
// colSize : ��
// a : �^�����p�ϐ�(�����l�͎g�p���Ȃ�)
template <class Type>
vector<vector<Type> >* getMatrixRegion(const unsigned long rowSize, const unsigned long colSize, const Type* a) {
    vector<vector<Type> >* mat;
    try{
        mat = new vector<vector<Type> >(rowSize, vector<Type>(colSize));
    } catch(std::bad_alloc ba) {    // �̈�m�ێ��s
        cout << "�������G���[���������܂����B";
        /* throw MatrixBadAlloc(&ba, rowSize, colSize); */
        throw ba;
    }
    return mat;
}

//===== ���Z�q�̓�d��`[2] =====//
// <<���Z�q
template <class Type> ostream& operator<<(ostream& ss, const Matrix<Type>& M) {return ss << M.toString();}
// ���Z
template <class Type> Matrix<Type> operator+(const Matrix<Type>& M, const Matrix<Type>& N) {
    unsigned long mr = M.row();
    unsigned long mc = M.col();
    unsigned long nr = N.row();
    unsigned long nc = N.col();
    if(mr != nr || mc != nc) throw MOperationErr<Type>(&M, &N);
    Matrix<Type> valMat(mr, mc);    // M�Ɠ����傫���̍s��
    for(unsigned long i = 0; i < mr; i++) {
        for(unsigned long j = 0; j < mc; j++) {
            valMat[i][j] = M[i][j] + N[i][j];
        }
    }
    return valMat;
}
// ���Z
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
// ���Z
template <class Type> Matrix<Type> operator-(const Matrix<Type>& M, const Matrix<Type>& N) {
    return M + (-1) * N;
}
// ���Z
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
// ��Z
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
// ��Z
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
// ����
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

//===== �s��v�Z�G���[ =====//
// �y�����v���z
// �E�T�C�Y���قȂ�s�񓯎m�̉��Z�E���Z���s����
// �E�񐔂ƍs������v���Ȃ�2�̍s������̏��Ԃŏ�Z���悤�Ƃ���
template < class Type >
class MOperationErr{
  const Matrix<Type>* identA;
  const Matrix<Type>* identB;
public:
  MOperationErr(const Matrix<Type>* a, const Matrix<Type>* b) : identA(a), identB(b) {}
};


#endif
