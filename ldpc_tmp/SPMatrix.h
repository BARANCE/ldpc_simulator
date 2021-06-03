/*
    alist��spmat�`���̑a�s��������N���X
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
    int rowSize;    // �s��
    int colSize;    // ��
    int maxWr;    // �s�d��
    int maxWc;    // ��d��
    vector<int> wr;    // �e�s�̍s�d�݂��L�^�����x�N�g��
    vector<int> wc;    // �e�s�̗�d�݂��L�^�����x�N�g��
    vector<vector<int> >* matCol;
    vector<vector<int> >* matRow;
    string fn;    // �t�@�C����


    static const int WDYM = 0;    // spmat���p���̃t���O
    //==== �R���X�g���N�^ ====//
    SPMatrix();
    explicit SPMatrix(string fileName);    // �t�@�C������ǂݍ���
    SPMatrix(string fileName, const int wdym);    // SPMat�`����ǂݍ���
    SPMatrix(int r, int c);
    SPMatrix(const Matrix<BinaryFiniteField>& M);    // �ϊ�
    SPMatrix(const SPMatrix& M);    // �R�s�[
    ~SPMatrix(){
        delete matCol;
        delete matRow;
    }    // �f�X�g���N�^

    //==== �擾���\�b�h ====//
    Matrix<BinaryFiniteField> getMatrix() const;
    string getFileName() const {return fn;}
    int row() const {return rowSize;}
    int col() const {return colSize;}
    int getMaxWr() const {return maxWr;}
    int getMaxWc() const {return maxWc;}
    vector<int> getWrVector() const {return wr;}
    vector<int> getWcVector() const {return wc;}
    vector<vector<int> >* getSPMatrixCol() const {return matCol;}
    vector<vector<int> >* getSPMatrixColCopy() const;    // �R�s�[���󂯎��o�[�W����
    vector<vector<int> >* getSPMatrixRow() const {return matRow;}
    vector<vector<int> >* getSPMatrixRowCopy() const;
    int isZero() const;    // �s�񂪗�s�񂩂ǂ���
    void printProgress(int n) const;

    //==== �s���؂�o�� ====//
    SPMatrix getRow(int index) const;
    SPMatrix getCol(int index) const;

    //==== �s���ϊ����� ====//
    SPMatrix transposed() const;    // �]�u�s��𓾂�

    //==== �s�f�[�^�E��f�[�^�̕ϊ� ====//
    void replaceMatRowFromMatCol();    // matCol����matRow�𐶐�
    void replaceMatColFromMatRow();    // matRow����matCol�𐶐�

    //==== �ݒ胁�\�b�h ====//
    void readAlist();    // �t�@�C������ǂݍ���Ńv���p�e�B�ɃZ�b�g
    void readAlist(string fileName);
    void readSpmat();    // �t�@�C������ǂݍ���Ńv���p�e�B�ɃZ�b�g
    void readSpmat(string fileName);
    void writeAlist();
    void writeAlist(string fileName);
    void writeSpmat();    // �v���p�e�B�����Ƀt�@�C���ɏ�������
    void writeSpmat(string fileName);
    void setFileName(string fileName) {fn = fileName;}
    void setMatrix(const Matrix<BinaryFiniteField>& M);
    void setMatrix(const SPMatrix& M);

    void ___setMaxWr(int n) {maxWr = n;}
    void ___setMaxWc(int n) {maxWc = n;}

    BinaryFiniteField getValue(int x, int y) const;    // x�sy��̒l���擾����B
    void setValue(int x, int y, const BinaryFiniteField value);    // x�sy���value��z�u����B

    void removeCycles(); // 4�T�C�N�����[�v�̏���

    //==== �����v�Z�p���\�b�h ====//
    int vecMatchCounter(const vector<int>& v1, const vector<int>& v2) const;    // v1��v2�ň�v���Ă���v�f�̐��𐔂���
    int checkValueInVec(const vector<int>& v, const int& value) const;    // v��value���܂܂�Ă��邩�𒲂ׂ�
    vector<int> vecExclusiveOr(const vector<int>& v1, const vector<int>& v2) const;    // v1��v2�̕Е������Ɋ܂܂��v�f�����o��
    vector<int> removeValueByVec(const vector<int>& v, const int& value) const;    // v����value����菜��

    //==== ���Z�q�̓�d��` ====//
    operator Matrix<BinaryFiniteField>() const;
    SPMatrix& operator=(const Matrix<BinaryFiniteField>& M);
    SPMatrix& operator=(const SPMatrix& M);
    friend SPMatrix operator+(const SPMatrix& M, const SPMatrix& N);    // ���Z
    friend SPMatrix operator-(const SPMatrix& M, const SPMatrix& N);
    friend SPMatrix operator*(const SPMatrix& M, const SPMatrix& N);    // ��Z
    friend SPMatrix operator*(const SPMatrix& M, const BinaryFiniteField& N);
    friend SPMatrix operator*(const BinaryFiniteField& M, const SPMatrix& N);
    friend int operator==(const SPMatrix& M, const SPMatrix& N);    // ��r

    static const int ZERO = 0;    // ��s��
    static const int E = 1;    // �P�ʍs��
    static const int ONE = 2;    // �S��s��
    friend int operator==(const SPMatrix& M, const int& identifier);    // ��s��E�P�ʍs��E�S��s�񂩂ǂ������`�F�b�N
    friend int operator!=(const SPMatrix& M, const int& identifier);

    //vector<int>& operator[](const int& i) {return (*matRow)[i];};


    //==== �G���[�N���X ====//
    class UndefErr{    // ����`�l�G���[
        SPMatrix* p;
    public:
        UndefErr(SPMatrix* ident) : p(ident) {}
    };

    class IOExceptions{    // �t�@�C���I�[�v���G���[
        SPMatrix* p;
    public:
        IOExceptions(SPMatrix* ident) : p(ident) {}
    };

    class EOFExceptions{    // �\�肵�Ă��Ȃ��Ƃ����EOF���o������
    };

    class FileFormatErr {    // �t�@�C���̃t�H�[�}�b�g(spmat�`��)���������Ȃ�
    };

    //==== ������ϊ� ====//
    stringstream& getStringstream() const {    // ���֐��Ŏg�p���镶����X�g���[���擾�p�֐�
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
