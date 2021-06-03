/*
    �񌳒l��v�f�Ƃ��Ď��a�s����\������N���X
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
    //==== �\���@��\���\���� ====//
    /*typedef enum ___CONS_MODE {
        GallagerRandomConstruction, ModifyedArrayCode, VMProtograph, LRLCode, ProtographCode, SpatiallyCoupledCode, ProgressiveEdgeGrowth, PEGSCC, RingedPEGSCC, ModifyedPEGSCC, RandomSCC, blockBasedSCC
    } ConsMode;*/

    //==== �a�s��ǂݍ��݃��\�b�h(�N���X���\�b�h) ====//
    static SPMatrix ReadSpmat(string fileName) {
        SPMatrix spm(fileName, SPMatrix::WDYM);
        return spm;
    }
    static SPMatrix ReadAlist(string fileName) {
        SPMatrix spm(fileName);
        return spm;
    }

    //==== �a�s��\�����\�b�h(��{�I�ɃN���X���\�b�h) ====//
    //==== �M�����K�[�̃����_���\���@ ====//
    // ��u���ɂ���Đ���LDPC�����̌����s��𐶐�����B��������s��̓t�������N�Ƃ͌���Ȃ��B
    // colSize : �s��̗�
    // wr : �s�d��
    // wc : ��d��
    static SPMatrix gallagerRandomConstruction(int colSize, int wr, int wc) {
        srand( (unsigned int)time( NULL ) );	// ����w��
		if(colSize % wc != 0) {
			cout << "[SMConstructor::gallagerRandomConstruction][warning]colSize��wc�Ŋ��肫��܂���B\n";
		}
		SPMatrix H(colSize * wc / wr, colSize);	// �s��{��
		
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
    static Matrix<BinaryFiniteField> gallagerRandomConstruction(int rowSize, int colSize) {    // ����������
        int tmp = gcd(rowSize, colSize);
        int wc = colSize / tmp;
        int wr = rowSize / tmp;
        return gallagerRandomConstruction(colSize, wr, wc);
    }

    // �ő����
    inline static int gcd( int m, int n ) {
        // �����ɂO������ꍇ�͂O��Ԃ�
        if ( ( 0 == m ) || ( 0 == n ) )
            return 0;

        // ���[�N���b�h�̕��@
        while( m != n ) {
            if ( m > n ) m = m - n;
            else         n = n - m;
        }
        return m;
    }//gcd

    //==== Modified Array Code(MAC) ====//
    // Array�����𔭓W�����������̌����s��𐶐�����B�����������̐��\���ǂ��B
    // j : �c�����̃u���b�N��
    // k : �������̃u���b�N��
    // p : �s����\������u���b�N(�����s��)�T�C�Y(p�~p)
    static Matrix<BinaryFiniteField> MAC(int paramJ, int paramK, int paramP) {
        Matrix<BinaryFiniteField> retMat(paramJ * paramP, paramK * paramP);    // jp�~kp�s��
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

    //==== �t�@���f�������h�s��Ɋ�Â��v���g�O���t ====//
    // �v���g�O���t�����ŗ��p���鏬�K�͂Ȍ����s��𐶐�����B���ꎩ�̂�Array�����ɂ�錟���s��Ɉ�v����B
    // paramJ : ��d��
    // paramK : �s�d��
    // paramP : �u���b�N�T�C�Y
    static Matrix<BinaryFiniteField> VMProtograph(int paramJ, int paramK, int paramP) {
        Matrix<BinaryFiniteField> retMat(paramJ * paramP, paramK * paramP);
        for(int i = 0; i < paramJ; i++) {
            for(int j = 0; j < paramK; j++) {
                // ��I�t�Z�b�g�̌���
                int c = i * paramP;    // c�s��
                int d = j * paramP;    // d���

                // �C���^���[�v���̌���
                int x = (i * j) % paramP;    // I�̊e�v�f��x�����E�ɏ���V�t�g

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

    //==== (l, r, L')��Ԍ��������̃v���g�O���t ====//
    // l : ��d��
    // r : �ő�s�d��
    // n : ������
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

    // ��d�݂ƍs���̊�l�A�񐔂���(l, r, L')��Ԍ��������̃v���g�O���t���쐬
    // m�͕K���s����菬�����Ȃ�
    static Matrix<BinaryFiniteField> LRL2(int m, int n, int wc) {
        int r = (int)((double)(n * wc) / (m - wc + 1));
        return LRL(wc, r, n);
    }
    // l,r,LHat�̃p�����[�^����쐬
    // l : ��d��
    // r : �ő�s�d��
    // LHat : ������
    static Matrix<BinaryFiniteField> LRLProtograph(unsigned long l, unsigned long r, unsigned long LHat) {
        unsigned long connectNum = LHat;
        unsigned long hr = l + (connectNum - 1);    // �s��
        unsigned long hc = r * connectNum / l;    // ��
        Matrix<BinaryFiniteField> retMat(hr, hc);
        for(unsigned long i = 0; i < connectNum; i++) {    // �����u���b�N����
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

    // LRL�v���g�O���t�̒[���𒆉����ɐڑ�����
    // modeValue : �ڑ�������ʒu��2�i���Ō���[1:���A2:�E�A4:��A8:��](��)13�Ȃ獶�㉺�ɐڑ�
    // centerCoverWidth : �����ɂȂ���ۂ̌�����
    // blockWidth : �����Ɍq����u���b�N�̑傫��(�傫���قǑ����̃m�[�h�Œ����Ɍq����)
    static Matrix<BinaryFiniteField> CCLRLP(unsigned long l, unsigned long r, unsigned long LHat, int modeValue, unsigned long centerCoverWidth, unsigned long blockWidth) {
        Matrix<BinaryFiniteField> retMat = LRLProtograph(l, r, LHat);
        // ���[�h�̉��
        vector<int> mode(4, 0);
        if(modeValue % 2 > 0) mode[0] = 1;    // ��
        if((modeValue % 4) / 2 > 0) mode[1] = 1;    // �E
        if((modeValue % 8) / 4 > 0) mode[2] = 1;    // ��
        if(modeValue / 8 > 0) mode[3] = 1;    // ��

        unsigned long hr = retMat.row();
        unsigned long hc = retMat.col();
        // ��`��`�F�b�N
        if(mode[0] == 1 || mode[1] == 1) {    // ���E�ڑ�
            if(centerCoverWidth > hr) centerCoverWidth = hr;
            if(blockWidth > hc) blockWidth = hc;
        }
        if(mode[2] == 1 || mode[3] == 1) {    // �㉺�ڑ�
            if(centerCoverWidth > hc) centerCoverWidth = hc;
            if(blockWidth > hr) blockWidth = hr;
        }

        // �ڑ��r�b�g�̔z�u
        if(mode[0] == 1) {    //����
            unsigned long x1 = hr / 2 - centerCoverWidth / 2;
            unsigned long x2 = x1 + centerCoverWidth;
            unsigned long y1 = 0;
            unsigned long y2 = blockWidth;
            placeBitsOnMatrix(&retMat, x1, x2, y1, y2);
        }
        if(mode[1] == 1) {    // �E��
            unsigned long x1 = hr / 2 - centerCoverWidth / 2;
            unsigned long x2 = x1 + centerCoverWidth;
            unsigned long y1 = hc - blockWidth;
            unsigned long y2 = hc;
            placeBitsOnMatrix(&retMat, x1, x2, y1, y2);
        }
        if(mode[2] == 1) {    // �㑤
            unsigned long x1 = 0;
            unsigned long x2 = blockWidth;
            unsigned long y1 = hc / 2 - centerCoverWidth / 2;
            unsigned long y2 = y1 + centerCoverWidth;
            placeBitsOnMatrix(&retMat, x1, x2, y1, y2);
        }
        if(mode[3] == 1) {    // ����
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


    //==== �v���g�O���t����(protograph code) ====//
    // �v���g�O���t�𕡐���R�s�[���A�m�[�h�Ԃ̒u�����s�����Ƃɂ���āA�����s��𐶐�����
    static SPMatrix ProtographCode(Matrix<BinaryFiniteField>* P, unsigned long T) {
        /*unsigned long pr = P->row();
        unsigned long pc = P->col();
        cout << "rowSize : " << pr * T << ", colSize = " << pc * T << "\n";

        Matrix<BinaryFiniteField> retMat(pr * T, pc * T);
        for(unsigned long i = 0; i < pr; i++) {
            for(unsigned long j = 0; j < pc; j++) {
                if((int)(*P)[i][j] == 1) {
                    // ��I�t�Z�b�g�̌���
                    unsigned long c = i * T;
                    unsigned long d = j * T;
                    vector<int> perm = getPerm(T);    // �u���x�N�g�����擾
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
                    // ��I�t�Z�b�g�̌���
                    unsigned long c = i * T;
                    unsigned long d = j * T;
                    vector<int> perm = getPerm(T);    // �u���x�N�g�����擾
                    for(unsigned long k = 0; k < T; k++) {
                        retMat.setValue(c + perm[k], d + k, 1);
                    }
                }
            }
        }
        return retMat;
    }

    //==== ��Ԍ�������(spatially-coupled code) ====//
    // �����Ȍ����s�����Ƃ��A�Ίp����ɑя�Ƀp���e�B�V���{�������Ԍ����s��𐶐�����B
    // P : ��ƂȂ鏬���ȃp���e�B�����s��
    // T : �R�s�[��
    static Matrix<BinaryFiniteField> SpatiallyCoupledCode(Matrix<BinaryFiniteField>* P, int T) {
        int pr = P->row();
        int pc = P->col();
        Matrix<BinaryFiniteField> retMat(pr * (T + 1), pc * T);
        for(int k = 0; k < T; k++) {
            // ��I�t�Z�b�g�̌���
            int c = k * pr;    // c�s��
            int d = k * pc;    // d���

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
    // �^�i�[�O���t�̓��a(girth)���\�Ȍ���傫������悤�ɁA�ӂ̐ڑ���I�����Ȃ���O���t���\�������@
    // rowSize : �����s��̍s��
    // colSize : �����s��̗�
    // valDim : �ϐ��m�[�h���Ƃ̎������i�[����x�N�g���B�傫����colSize�ƈ�v���邱��
    static Matrix<BinaryFiniteField> ProgressiveEdgeGrowth(int rowSize, int colSize, const vector<int>& valDim) {
        int hr = rowSize;
        int hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> checkDim(hr, 0);    // �`�F�b�N�m�[�h�̎������i�[����x�N�g��
        for(int i = 0; i < hc; i++) {
            // 1�{��
            if(valDim[i] > 0) {
                int checkSelect = 0;    // ���ƂȂ�`�F�b�N�m�[�h
                for(int k = 1; k < hr; k++) {
                    if(checkDim[checkSelect] > checkDim[k]) {
                        checkSelect = k;    // �ł������������̃m�[�h��I��
                    }
                }
                retMat[checkSelect][i] = 1;    // �ӂ�z�u
                checkDim[checkSelect]++;    // �`�F�b�N�m�[�h�̎����㏸
            }

            // 2�{�ڈȍ~
            for(int j = 1; j < valDim[i]; j++) {
                vector<int> prevCheckVec(0);    // ��O�̃`�F�b�N�m�[�h�̔ԍ����i�[����x�N�g��
                vector<int> nowCheckVec(0);    // ���݂̃`�F�b�N�m�[�h�̔ԍ����i�[����x�N�g��
                vector<int> valList(0);    // ���ڂ���ϐ��m�[�h�̃��X�g
                vector<int> nextValList(0);
                vector<int> valEnd(0);    // �`�F�b�N�ςݕϐ��m�[�h�̃��X�g
                valList.push_back(i);

                //==== �������烋�[�v���� ====//
                while(1) {
                    int valNode = valList[valList.size() - 1];
                    contains(valEnd, valNode);
                    valList.pop_back();

                    // �ڑ�����`�F�b�N�m�[�h�̃��X�g
                    for(int k = 0; k < hr; k++) {
                        if((int)retMat[k][valNode] == 1) {
                            //cout << "nowCheckVec.size = " << nowCheckVec.size() << "��";
                            contains(nowCheckVec, k);
                            //cout << nowCheckVec.size() << "\n";
                        }
                    }

                    // �`�F�b�N�m�[�h���X�g�ɐڑ�����ϐ��m�[�h�̃��X�g
                    for(int k = 0; k < nowCheckVec.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[nowCheckVec[k]][l] == 1 && l != valNode) {    // �������s��
                                if(!contains2(valEnd, l)) {
                                    contains(nextValList, l);
                                }
                            }
                        }
                    }

                    //cout << prevCheckVec.size() << " �� " << nowCheckVec.size() << "\n";

                    //cout << "nowCheckVec = " << nowCheckVec;

                    if(nowCheckVec.size() == hr) {
                        //cout << "���^��";
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
                            // �`�F�b�N�m�[�h���X�g���X�V
                            if(eqVec(prevCheckVec, nowCheckVec)) {
                                //cout << "����";
                                break;
                            }
                            //cout << "��";
                            break;
                        }
                    }

                }
                //==== ���[�v���������܂� ====//

                // ��W�����擾
                vector<int> comp = complement(nowCheckVec, hr);

                /*cout << "comp = [";
                for(int b = 0; b < comp.size(); b++) {
                    cout << comp[b] << " ";
                }
                cout << "]\n";*/

                int checkSelect = comp[0];
                for(int k = 1; k < comp.size(); k++) {
                    if(checkDim[checkSelect] > checkDim[comp[k]]) {
                        checkSelect = comp[k];    // �ł������������̃m�[�h��I��
                    }
                }
                retMat[checkSelect][i] = 1;    // �ӂ�z�u
                checkDim[checkSelect]++;    // �`�F�b�N�m�[�h�̎����㏸
            }
        }
        return retMat;
    }

    //==== PEG��Ԍ�������(PEG-SCC) ====//
    // PEG�ɕӂ�ڑ�����ۂ̐����t����������
    // bandWidth : �т̑傫��
    static Matrix<BinaryFiniteField> PEGSCC(int rowSize, int colSize, const vector<int>& valDim, int bandWidth) {
        int hr = rowSize;
        int hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> checkDim(hr, 0);    // �`�F�b�N�m�[�h�̎������i�[����x�N�g��
        for(int i = 0; i < hc; i++) {
            //cout << "i = " << i << "\n";
            // ��I�t�Z�b�g�̌���
            int c = ((double)(hr - bandWidth + 1) / hc) * i;
            //cout << "c = " << c << "\n";

            // 1�{��
            if(valDim[i] > 0) {
                int checkSelect = c;    // ���ƂȂ�`�F�b�N�m�[�h
                for(int k = c; k < c + bandWidth; k++) {
                    if(checkDim[checkSelect] > checkDim[k]) {
                        checkSelect = k;    // �ł������������̃m�[�h��I��
                    }
                }
                //cout << "checkSekect : " << checkSelect << "\n";
                retMat[checkSelect][i] = 1;    // �ӂ�z�u
                checkDim[checkSelect]++;    // �`�F�b�N�m�[�h�̎����㏸
            }

            // 2�{�ڈȍ~
            for(int j = 1; j < valDim[i]; j++) {
                vector<int> prevCheckVec(0);    // ��O�̃`�F�b�N�m�[�h�̔ԍ����i�[����x�N�g��
                vector<int> nowCheckVec(0);    // ���݂̃`�F�b�N�m�[�h�̔ԍ����i�[����x�N�g��
                vector<int> valList(0);    // ���ڂ���ϐ��m�[�h�̃��X�g
                vector<int> nextValList(0);
                vector<int> valEnd(0);    // �`�F�b�N�ςݕϐ��m�[�h�̃��X�g
                valList.push_back(i);

                //==== �������烋�[�v���� ====//
                while(1) {
                    int valNode = valList[valList.size() - 1];
                    contains(valEnd, valNode);
                    valList.pop_back();

                    // �ڑ�����`�F�b�N�m�[�h�̃��X�g
                    for(int k = 0; k < hr; k++) {
                        if((int)retMat[k][valNode] == 1) {
                            //cout << "nowCheckVec.size = " << nowCheckVec.size() << "��";
                            contains(nowCheckVec, k);
                            //cout << nowCheckVec.size() << "\n";
                        }
                    }

                    // �`�F�b�N�m�[�h���X�g�ɐڑ�����ϐ��m�[�h�̃��X�g
                    for(int k = 0; k < nowCheckVec.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[nowCheckVec[k]][l] == 1 && l != valNode) {    // �������s��
                                if(!contains2(valEnd, l)) {
                                    contains(nextValList, l);
                                }
                            }
                        }
                    }

                    //cout << prevCheckVec.size() << " �� " << nowCheckVec.size() << "\n";

                    //cout << "nowCheckVec = " << nowCheckVec;

                    if(containMinToMax(nowCheckVec, c, c + bandWidth)) {
                    //if(nowCheckVec.size() == hr) {    // ������ς���̂��ȁH
                        //cout << "���^��";
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
                            // �`�F�b�N�m�[�h���X�g���X�V
                            if(eqVec(prevCheckVec, nowCheckVec)) {
                                //cout << "����";
                                break;
                            }
                            //cout << "��";
                            break;
                        }
                    }

                }
                //==== ���[�v���������܂� ====//

                // ��W�����擾
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
                        checkSelect = comp[k];    // �ł������������̃m�[�h��I��
                    }
                }
                retMat[checkSelect][i] = 1;    // �ӂ�z�u
                checkDim[checkSelect]++;    // �`�F�b�N�m�[�h�̎����㏸
            }
        }
        return retMat;
    }

    //==== Ringed-PEG��Ԍ�������(R-PEGSCC) ====//
    static Matrix<BinaryFiniteField> RPEGSCC(int rowSize, int colSize, const vector<int>& valDim, int bandWidth) {
        int hr = rowSize;
        int hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> checkDim(hr, 0);    // �`�F�b�N�m�[�h�̎������i�[����x�N�g��
        for(int i = 0; i < hc; i++) {
            //cout << "i = " << i << "\n";
            // ��I�t�Z�b�g�̌���
            int c = ((double)hr / hc) * i;
            //cout << "c = " << c << "\n";

            // 1�{��
            if(valDim[i] > 0) {
                int checkSelect = c;    // ���ƂȂ�`�F�b�N�m�[�h
                for(int k = c; k < c + bandWidth; k++) {
                    int sel = k;
                    if(sel >= hr) sel -= hr;
                    if(checkDim[checkSelect] > checkDim[sel]) {
                        checkSelect = sel;    // �ł������������̃m�[�h��I��
                    }
                }
                //cout << "checkSekect : " << checkSelect << "\n";
                retMat[checkSelect][i] = 1;    // �ӂ�z�u
                checkDim[checkSelect]++;    // �`�F�b�N�m�[�h�̎����㏸
            }

            // 2�{�ڈȍ~
            for(int j = 1; j < valDim[i]; j++) {
                vector<int> prevCheckVec(0);    // ��O�̃`�F�b�N�m�[�h�̔ԍ����i�[����x�N�g��
                vector<int> nowCheckVec(0);    // ���݂̃`�F�b�N�m�[�h�̔ԍ����i�[����x�N�g��
                vector<int> valList(0);    // ���ڂ���ϐ��m�[�h�̃��X�g
                vector<int> nextValList(0);
                vector<int> valEnd(0);    // �`�F�b�N�ςݕϐ��m�[�h�̃��X�g
                valList.push_back(i);

                //==== �������烋�[�v���� ====//
                while(1) {
                    int valNode = valList[valList.size() - 1];
                    contains(valEnd, valNode);
                    valList.pop_back();

                    // �ڑ�����`�F�b�N�m�[�h�̃��X�g
                    for(int k = 0; k < hr; k++) {
                        if((int)retMat[k][valNode] == 1) {
                            //cout << "nowCheckVec.size = " << nowCheckVec.size() << "��";
                            contains(nowCheckVec, k);
                            //cout << nowCheckVec.size() << "\n";
                        }
                    }

                    // �`�F�b�N�m�[�h���X�g�ɐڑ�����ϐ��m�[�h�̃��X�g
                    for(int k = 0; k < nowCheckVec.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[nowCheckVec[k]][l] == 1 && l != valNode) {    // �������s��
                                if(!contains2(valEnd, l)) {
                                    contains(nextValList, l);
                                }
                            }
                        }
                    }

                    //cout << prevCheckVec.size() << " �� " << nowCheckVec.size() << "\n";

                    //cout << "nowCheckVec = " << nowCheckVec;
                    if(c + bandWidth < hr) {
                        if(containMinToMax(nowCheckVec, c, c + bandWidth)) {
                        //if(nowCheckVec.size() == hr) {    // ������ς���̂��ȁH
                            //cout << "���^��";
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
                            // �`�F�b�N�m�[�h���X�g���X�V
                            if(eqVec(prevCheckVec, nowCheckVec)) {
                                //cout << "����";
                                break;
                            }
                            //cout << "��";
                            break;
                        }
                    }

                }
                //==== ���[�v���������܂� ====//

                // ��W�����擾
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
                        checkSelect = comp[k];    // �ł������������̃m�[�h��I��
                    }
                }
                retMat[checkSelect][i] = 1;    // �ӂ�z�u
                checkDim[checkSelect]++;    // �`�F�b�N�m�[�h�̎����㏸
            }
        }
        return retMat;
    }

    //==== �Ӕz�u������ύX����PEG��Ԍ�������(PEG-SCC) ====//
    // PEG�ɕӂ�ڑ�����ۂ̐����t����������
    // bandWidth : �т̑傫��
    static Matrix<BinaryFiniteField> MPEGSCC(int rowSize, int colSize, const vector<int>& valDim, int bandWidth) {
        int hr = rowSize;
        int hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> checkDim(hr, 0);    // �`�F�b�N�m�[�h�̎������i�[����x�N�g��
        for(int a = 0; a < hc; a++) {
            int i;    // �z�u�ʒu��O����끨�O���c�Ƃ���
            if(a % 2 == 0) {
                i = a / 2;
            } else {
                i = hc - (a / 2) - 1;
            }

            //cout << "i = " << i << "\n";
            // ��I�t�Z�b�g�̌���
            int c = ((double)(hr - bandWidth + 1) / hc) * i;
            //cout << "c = " << c << "\n";

            // 1�{��
            if(valDim[i] > 0) {
                int checkSelect = c;    // ���ƂȂ�`�F�b�N�m�[�h
                for(int k = c; k < c + bandWidth; k++) {
                    if(checkDim[checkSelect] > checkDim[k]) {
                        checkSelect = k;    // �ł������������̃m�[�h��I��
                    }
                }
                //cout << "checkSekect : " << checkSelect << "\n";
                retMat[checkSelect][i] = 1;    // �ӂ�z�u
                checkDim[checkSelect]++;    // �`�F�b�N�m�[�h�̎����㏸
            }

            // 2�{�ڈȍ~
            for(int j = 1; j < valDim[i]; j++) {
                vector<int> prevCheckVec(0);    // ��O�̃`�F�b�N�m�[�h�̔ԍ����i�[����x�N�g��
                vector<int> nowCheckVec(0);    // ���݂̃`�F�b�N�m�[�h�̔ԍ����i�[����x�N�g��
                vector<int> valList(0);    // ���ڂ���ϐ��m�[�h�̃��X�g
                vector<int> nextValList(0);
                vector<int> valEnd(0);    // �`�F�b�N�ςݕϐ��m�[�h�̃��X�g
                valList.push_back(i);

                //==== �������烋�[�v���� ====//
                while(1) {
                    int valNode = valList[valList.size() - 1];
                    contains(valEnd, valNode);
                    valList.pop_back();

                    // �ڑ�����`�F�b�N�m�[�h�̃��X�g
                    for(int k = 0; k < hr; k++) {
                        if((int)retMat[k][valNode] == 1) {
                            //cout << "nowCheckVec.size = " << nowCheckVec.size() << "��";
                            contains(nowCheckVec, k);
                            //cout << nowCheckVec.size() << "\n";
                        }
                    }

                    // �`�F�b�N�m�[�h���X�g�ɐڑ�����ϐ��m�[�h�̃��X�g
                    for(int k = 0; k < nowCheckVec.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[nowCheckVec[k]][l] == 1 && l != valNode) {    // �������s��
                                if(!contains2(valEnd, l)) {
                                    contains(nextValList, l);
                                }
                            }
                        }
                    }

                    //cout << prevCheckVec.size() << " �� " << nowCheckVec.size() << "\n";

                    //cout << "nowCheckVec = " << nowCheckVec;

                    if(containMinToMax(nowCheckVec, c, c + bandWidth)) {
                    //if(nowCheckVec.size() == hr) {    // ������ς���̂��ȁH
                        //cout << "���^��";
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
                            // �`�F�b�N�m�[�h���X�g���X�V
                            if(eqVec(prevCheckVec, nowCheckVec)) {
                                //cout << "����";
                                break;
                            }
                            //cout << "��";
                            break;
                        }
                    }

                }
                //==== ���[�v���������܂� ====//

                // ��W�����擾
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
                        checkSelect = comp[k];    // �ł������������̃m�[�h��I��
                    }
                }
                retMat[checkSelect][i] = 1;    // �ӂ�z�u
                checkDim[checkSelect]++;    // �`�F�b�N�m�[�h�̎����㏸
            }
        }
        return retMat;
    }

    //==== �����_����Ԍ������� ====//
    // bandWidth : �т̑傫��
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

    //==== �u���b�N�z�u��Ԍ������� ====//
    // p�~q�s����΂߂ɁA���������`�F�b�N�m�[�h���d�Ȃ�悤�ɔz�u���Č����s����\������
    // �e�s�񂲂Ƃ̃T�C�Y���ςƂ���
    static Matrix<BinaryFiniteField> BlockSCC(vector<Matrix<BinaryFiniteField> > blocks, unsigned long conSize) {
        int blockNum = blocks.size();    // �z�u����u���b�N�̐�
        unsigned long hr = 0;
        unsigned long hc = 0;
        for(int i = 0; i < blockNum; i++) {
            hr += blocks[i].row();
            hc += blocks[i].col();
            if(blocks[i].row() < conSize) {
                cout << "ERROR : �u���b�N�s��̍s�����A��������菬�������ߌ����s����\���ł��܂���B\n";
            }
        }
        hr -= conSize * (blockNum - 1); // ����������������������
        cout << "hr = " << hr << ", hc = " << hc << "\n";
        Matrix<BinaryFiniteField> retMat(hr, hc);
        unsigned long hrSum = 0;
        unsigned long hcSum = 0;
        for(int i = 0; i < blockNum; i++) {
            unsigned long rbase = hrSum - (i * blockNum);    // ��ʒu(�s)
            unsigned long cbase = hcSum;    // ��ʒu(��)

            unsigned long shr = blocks[i].row();
            unsigned long shc = blocks[i].col();
            for(int j = 0; j < shr; j++) {
                for(int k = 0; k < shc; k++) {
                    retMat[rbase + j][cbase + k] = blocks[i][j][k];
                }
            }

            //==== ��ʒu���X�V ====//
            hrSum += shr;
            hcSum += shc;
        }
        return retMat;
    }

    //==== Progressive Edge Growth(PEG) ������ ====//
    // �^�i�[�O���t�̓��a(girth)���\�Ȍ���傫������悤�ɁA�ӂ̐ڑ���I�����Ȃ���O���t���\�������@
    // rowSize : �����s��̍s��
    // colSize : �����s��̗�
    // valDim : �ϐ��m�[�h���Ƃ̎������i�[����x�N�g���B�傫����colSize�ƈ�v���邱��
    static Matrix<BinaryFiniteField> PEG2(int rowSize, int colSize, const vector<int>& valDim) {
        int hr = rowSize;
        int long hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> nowCheckDim(hr);    // �\�����̊e�`�F�b�N�m�[�h�̌��݂̎�����\���x�N�g��
        vector<int> distance(hr);    // �C�ӂ̕ϐ��m�[�h�̕Ӑڑ��X�e�b�v�ɂ����āA�e�`�F�b�N�m�[�h�ւ̋������ꎞ�i�[����x�N�g���B

        //==== �`�F�b�N�m�[�h�̎�����0�ɏ����� ====//
        for(int i = 0; i < hr; i++) nowCheckDim[i] = 0;

        for(int i = 0; i < hc; i++) {
            int selectedVid = i;    // �I�����ꂽ�ϐ��m�[�h�B��ɕύX���邱�Ƃ�z��

            //==== i�Ԗڂ̕ϐ��m�[�h����e�`�F�b�N�m�[�h�ւ̋����𖳌���(-1)�ɏ����� ====//
            for(int j = 0; j < hr; j++) distance[j] = -1;

            //==== �`�F�b�N�m�[�h�̑I���\�͈�(PEG��Ԍ��������Ȃǂŗ��p) ====//
            int startPos = 0;
            int endPos = hr;

            //==== 1�{�ڂ͎����̈�ԒႢ�`�F�b�N�m�[�h�ɐڑ� ====//
            int minDimId = startPos;
            int minDim = nowCheckDim[startPos];
            for(int j = startPos + 1; j < endPos; j++) {
                if(nowCheckDim[j] < minDim) {
                    minDimId = j;
                    minDim = nowCheckDim[j];
                }
            }
            retMat[minDimId][i] = 1;    // 1�{�ڂ̕ӂ�ݒu
            nowCheckDim[minDimId]++;

            //==== �������X�g�쐬(�T��)�̂��߂̏����m�[�h ====//
            int depth = 0;    // ���݂̐[��
            int maxDepth = 10;    // �ő�[��
            vector<int> ovList(0);
            ovList.push_back(i);    // �����m�[�h�͒����ϐ��m�[�h

            //==== �ϐ��m�[�h�̎����𖞂����܂ŌJ��Ԃ� ====//
            for(int j = 1; j < valDim[i]; j++) {
                //==== �ł������`�F�b�N�m�[�h�𔭌����邽�߂̃��[�v ====//
                for(depth = 0; depth < maxDepth; depth++) {
                    //==== �`�F�b�N�m�[�h�ւ̋������X�V ====//
                    vector<int> children(0);
                    for(int l = 0; l < ovList.size(); l++) {
                        for(int k = 0; k < hr; k++) {
                            if((int)retMat[k][ovList[l]] == 1) {    // �ڑ����Ă���`�F�b�N�m�[�h�ɑ΂���
                                if(distance[k] == -1) {    // ���������ݒ�(������)���X�V
                                    children.push_back(k);
                                    distance[k] = depth;    // �������X�V
                                }else if(distance[k] != -1 && depth < distance[k]) {    // ���݂̋���(depth)���傫���������ݒ肳��Ă��遨�X�V
                                    children.push_back(k);
                                    distance[k] = depth;    // �������X�V
                                }
                            }
                        }
                    }
                    //==== �`�F�b�N�m�[�h���X�gchildren�̏d������菜�� ====//
                    vector<int> childrenRev(0);    // �d������菜����children
                    for(int k = 0; k < children.size(); k++) {
                        bool frag = false;    // �d���������true
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

                    //==== �`�F�b�N�m�[�h���X�gchildren�ɐڑ�����ϐ��m�[�h���Z�o ====//
                    vector<int> connectedValueList(0);
                    for(int k = 0; k < childrenRev.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[childrenRev[k]][l] == 1) {
                                connectedValueList.push_back(l);
                            }
                        }
                    }

                    //==== connectedValueList�̏d������菜��(ovList�Ƃ̏d��������) ====//
                    vector<int> connectedRev(0);
                    for(int k = 0; k < connectedValueList.size(); k++) {
                        bool frag = false;    // �d���������true
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
                    if(ovList.empty()) break;    // ���ɒ��ׂ�ϐ��m�[�h���X�g����

                }
                //==== �������ł��傫���m�[�h���X�g���Z�o ====//
                int maxValue = 1;    // �ő勗��(������Ȃ��-1)
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

                //==== ���̒��ōł��������������`�F�b�N�m�[�h���Z�o ====//
                int cid = highGirthList[0];
                int cdim = nowCheckDim[highGirthList[0]];
                for(int j = 1; j < highGirthList.size(); j++) {
                    if(nowCheckDim[highGirthList[j]] < cdim) {
                        cid = highGirthList[j];    // �X�V
                        cdim = nowCheckDim[highGirthList[j]];
                    }
                }

                //==== �ӂ�ݒu ====//
                if((int)retMat[cid][i] == 1) {
                    cout << "distance[" << cid << "] = " << distance[cid] << "\n";
                    cout << "retMat[" << cid << "][" << i << "]�ɑ��d��\n";

                }
                retMat[cid][i] = 1;
                nowCheckDim[cid]++;    // �����㏸

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

    //==== PEG��Ԍ�������2 ====//
    static Matrix<BinaryFiniteField> PEGSCC2(int rowSize, int colSize, const vector<int>& valDim, int bandWidth) {
        int hr = rowSize;
        int long hc = colSize;
        Matrix<BinaryFiniteField> retMat(hr, hc);
        vector<int> nowCheckDim(hr);    // �\�����̊e�`�F�b�N�m�[�h�̌��݂̎�����\���x�N�g��
        vector<int> distance(hr);    // �C�ӂ̕ϐ��m�[�h�̕Ӑڑ��X�e�b�v�ɂ����āA�e�`�F�b�N�m�[�h�ւ̋������ꎞ�i�[����x�N�g���B

        //==== �`�F�b�N�m�[�h�̎�����0�ɏ����� ====//
        for(int i = 0; i < hr; i++) nowCheckDim[i] = 0;

        for(int i = 0; i < hc; i++) {
            int selectedVid = i;    // �I�����ꂽ�ϐ��m�[�h�B��ɕύX���邱�Ƃ�z��

            //==== i�Ԗڂ̕ϐ��m�[�h����e�`�F�b�N�m�[�h�ւ̋����𖳌���(-1)�ɏ����� ====//
            for(int j = 0; j < hr; j++) distance[j] = -1;

            //==== �`�F�b�N�m�[�h�̑I���\�͈�(PEG��Ԍ��������Ȃǂŗ��p) ====//
            int c = ((double)(hr - bandWidth + 1) / hc) * i;
            int startPos = c;
            int endPos = c + bandWidth;

            //==== 1�{�ڂ͎����̈�ԒႢ�`�F�b�N�m�[�h�ɐڑ� ====//
            int minDimId = startPos;
            int minDim = nowCheckDim[startPos];
            for(int j = startPos + 1; j < endPos; j++) {
                if(nowCheckDim[j] < minDim) {
                    minDimId = j;
                    minDim = nowCheckDim[j];
                }
            }
            retMat[minDimId][i] = 1;    // 1�{�ڂ̕ӂ�ݒu
            nowCheckDim[minDimId]++;

            //==== �������X�g�쐬(�T��)�̂��߂̏����m�[�h ====//
            int depth = 0;    // ���݂̐[��
            int maxDepth = 10;    // �ő�[��
            vector<int> ovList(0);
            ovList.push_back(i);    // �����m�[�h�͒����ϐ��m�[�h

            //==== �ϐ��m�[�h�̎����𖞂����܂ŌJ��Ԃ� ====//
            for(int j = 1; j < valDim[i]; j++) {
                //==== �ł������`�F�b�N�m�[�h�𔭌����邽�߂̃��[�v ====//
                for(depth = 0; depth < maxDepth; depth++) {
                    //==== �`�F�b�N�m�[�h�ւ̋������X�V ====//
                    vector<int> children(0);
                    for(int l = 0; l < ovList.size(); l++) {
                        for(int k = 0; k < hr; k++) {
                            if((int)retMat[k][ovList[l]] == 1) {    // �ڑ����Ă���`�F�b�N�m�[�h�ɑ΂���
                                if(distance[k] == -1) {    // ���������ݒ�(������)���X�V
                                    children.push_back(k);
                                    distance[k] = depth;    // �������X�V
                                }else if(distance[k] != -1 && depth < distance[k]) {    // ���݂̋���(depth)���傫���������ݒ肳��Ă��遨�X�V
                                    children.push_back(k);
                                    distance[k] = depth;    // �������X�V
                                }
                            }
                        }
                    }
                    //==== �`�F�b�N�m�[�h���X�gchildren�̏d������菜�� ====//
                    vector<int> childrenRev(0);    // �d������菜����children
                    for(int k = 0; k < children.size(); k++) {
                        bool frag = false;    // �d���������true
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

                    //==== �`�F�b�N�m�[�h���X�gchildren�ɐڑ�����ϐ��m�[�h���Z�o ====//
                    vector<int> connectedValueList(0);
                    for(int k = 0; k < childrenRev.size(); k++) {
                        for(int l = 0; l < hc; l++) {
                            if((int)retMat[childrenRev[k]][l] == 1) {
                                connectedValueList.push_back(l);
                            }
                        }
                    }

                    //==== connectedValueList�̏d������菜��(ovList�Ƃ̏d��������) ====//
                    vector<int> connectedRev(0);
                    for(int k = 0; k < connectedValueList.size(); k++) {
                        bool frag = false;    // �d���������true
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
                    if(ovList.empty()) break;    // ���ɒ��ׂ�ϐ��m�[�h���X�g����

                }
                //==== �������ł��傫���m�[�h���X�g���Z�o ====//
                int maxValue = 1;    // �ő勗��(������Ȃ��-1)
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

                //==== ���̒��ōł��������������`�F�b�N�m�[�h���Z�o ====//
                int cid = highGirthList[0];
                int cdim = nowCheckDim[highGirthList[0]];
                for(int j = 1; j < highGirthList.size(); j++) {
                    if(nowCheckDim[highGirthList[j]] < cdim) {
                        cid = highGirthList[j];    // �X�V
                        cdim = nowCheckDim[highGirthList[j]];
                    }
                }

                //==== �ӂ�ݒu ====//
                if((int)retMat[cid][i] == 1) {
                    cout << "distance[" << cid << "] = " << distance[cid] << "\n";
                    cout << "retMat[" << cid << "][" << i << "]�ɑ��d��\n";

                }
                retMat[cid][i] = 1;
                nowCheckDim[cid]++;    // �����㏸

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

    //==== LRL���������ĐV�����v���g�O���t���\�� ====//
    // lrl1 : �O����LRL��Ԍ��������v���g�O���t���w��
    // lrl2 : ������LRL��Ԍ��������v���g�O���t���w��
    // splitId : �����ӏ��̕ϐ��m�[�hid
    static Matrix<BinaryFiniteField> connectLRL(Matrix<BinaryFiniteField>* lrl1, Matrix<BinaryFiniteField>* lrl2, unsigned long splitId) {
        unsigned long hr1 = lrl1->row();
        unsigned long hr2 = lrl2->row();
        unsigned long hc1 = lrl1->col();
        unsigned long hc2 = lrl2->col();
        if(splitId < 0 || splitId > hc1) throw ConstructExeption();

        //==== LRL1�̐؂���̍��E�̃f�[�^���擾 ====//
        // ����
        BandData left1 = {0, 0, 0};
        if(splitId != 0) {
            left1 = getBandData(lrl1, splitId - 1);
        }
        // �E��
        BandData right1 = {0, 0, 0};
        if(splitId != hc1) {
            right1 = getBandData(lrl1, splitId);
        }

        //==== �ڑ�����LRL2�̐؂���̃f�[�^���擾 ====//
        // ����
        BandData left2 = getBandData(lrl2, 0);
        // �E��
        BandData right2 = getBandData(lrl2, hc2 - 1);

    }

    //==== Loop��Ԍ������� ====//
    static Matrix<BinaryFiniteField> loopLDPC(const vector<LoopLDPC::LBandData>& bandData) {
        LoopLDPC lldpc(bandData);
        return lldpc.generate();
    }

    //==== �����s��t�@�C����ǂݍ��� ====//
    static Matrix<BinaryFiniteField> loadMatrixFile(const string fileName) {
        SPMatrix H(fileName);
        Matrix<BinaryFiniteField> retMat = H;

        //==== �����Ԍ��������̃e�X�g(�f�o�b�O�p) ====//
        // ��������
        /*unsigned long hr = 31;
        unsigned long hc = 58;
        Matrix<BinaryFiniteField> PosMat(hr, hc);*/

        // �߂��Ƃ���͋߂��Ɛڑ�
        /*unsigned long d = 2;
        for(unsigned long i = 0; i < d; i++) {
            PosMat[30 - i][56 - i * 2] = 1;
            PosMat[30 - i][57 - i * 2] = 1;
        }*/

        // �����Ƃ���قǋ߂��Ɛڑ�
        /*unsigned long e = 2;
        for(unsigned long i = 0; i < e; i++) {
            PosMat[30 - i][58 - e * 2 + i * 2] = 1;
            PosMat[30 - i][59 - e * 2 + i * 2] = 1;
        }*/

        // �O�p�`�h��Ԃ�
        /*unsigned long f = 10;
        for(unsigned long i = 0; i < f; i++) {
            for(unsigned long j = 0; j < f * 2 - i * 2; j++) {
                PosMat[30 - i][58 - f * 2 + i * 2 + j] = 1;
            }
        }*/

        // EX�`�F�b�N�m�[�h��L�΂��ڑ�
        /*unsigned long g = 29;
        for(unsigned long i = 0; i < g; i++) {
            PosMat[30][58 - g * 2 + i * 2] = 1;
            PosMat[30][59 - g * 2 + i * 2] = 1;
        }*/

        // EX�ϐ��m�[�h��L�΂��ڑ�
        /*unsigned long h = 29;
        for(unsigned long i = 0; i < h; i++) {
            PosMat[30 - i][56] = 1;
            PosMat[30 - i][57] = 1;
        }*/

        // ����������
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
        // �����܂�

        return retMat;
    }

    //==== �v���g�O���t�����Ԍ��������̍쐬 ====//
    // �Ɨ��P�[�X�v���g�O���t�̍s��t�@�C����ǂݍ��݁A�ڑ�����ǉ�
    // fileName : [string]�ǂݍ��ރt�@�C���̃p�X�{���O
    // l : [int]��d��
    // r : [int]�ő�s�d��
    // L : [int]������
    // copyNum : [int]���̖{��
    // connectType : [int]�ڑ����@
    // connectSize : [int]�ڑ����̑傫��
    static Matrix<BinaryFiniteField> ConnectIndependentLRLProtograph (int l, int r, int L, int copyNum, int connectSize) {
        //Matrix<BinaryFiniteField> retMat = SMConstructor::loadMatrixFile(fileName);
        unsigned long protHr = l + L - 1;    // lrL�v���g�O���t�̍s��
        unsigned long protHc = r * L / l;    // lrL�v���g�O���t�̗�
        unsigned long hr = protHr * copyNum;    // retMat�̍s��
        unsigned long hc = protHc * copyNum;    // retMat�̗�
        Matrix<BinaryFiniteField> retMat(hr, hc);

        // ��ѐ��G���[�`�F�b�N
        /*if(hr != retMat.row() || hc != retMat.col()) {    // �w�肵�����e�ƃv���g�O���t���قȂ�
            cout << "[SMConstructor::loadAndConnectIndependentLRLProtograph][error]�ǂݍ��񂾃v���g�O���t�̍s���E�񐔂ƁA�w�肳�ꂽ�p�����[�^(l,r,L,copyNum)����Z�o�����s���E�񐔂���v���܂���B�v���g�O���t�̎w�肪����Ă��邩�A�p�����[�^�𐳂����w�肵�Ă��������B\n";
            Matrix<BinaryFiniteField> nu(0,0);
            return nu;
        }
        if(connectSize > L) {
            cout << "[SMConstructor::loadAndConnectIndependentLRLProtograph][error]�ڑ����̑傫��connectSize�́A������L�ȉ��łȂ���΂Ȃ�܂���B\n";
            Matrix<BinaryFiniteField> nu(0,0);
            return nu;
        }*/

        // �e�u���b�N��ʒu�v�Z
        // �����ӁF�ڑ���C�̊�ʒu�͈ȉ��̒ʂ�
        // a�Ԗڂ̃v���g�O���t��(a + 1 mod copyNum)�Ԗڂ̃v���g�O���t�Ɛڑ�����ꍇ�A�ڑ���C��
        // ��s�FrowStartPow[a] ... rowStartPow[a + 1] - 1
        // ���FcolStartPow[(a + 1) % copyNum] ... colStartPow[(a + 1) % copyNum + 1] - 1
        // �ƂȂ�_�ɒ��ӁB
        vector<unsigned long> rowStartPos(copyNum + 1, 0);
        vector<unsigned long> colStartPos(copyNum + 1, 0);
        for(int i = 0 ; i < copyNum + 1; i++) {
            rowStartPos[i] = protHr * i;
            colStartPos[i] = protHc * i;
        }

        // lrL��Ԍ��������̃v���g�O���t��z�u
        for(int i = 0; i < copyNum; i++) {
            int baseX = rowStartPos[i];
            int baseY = colStartPos[i];
            int bHeight = l;
            int bWidth = r / l;    // �e�u���b�N�̉���
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

        // �ڑ����̍쐬
        // i�Ԗڂ̃u���b�N�ƁA(i + 1) % copyNum �Ԗڂ̃u���b�N��ڑ�
        for(int i = 0; i < copyNum; i++) {
            // i�Ԗڂ̐ڑ����̔z�u�ʒu����
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

    //==== �v���g�O���t�����Ԍ��������̍쐬(�������ɉ����Đڑ����������I�ɑ傫������Version) ====//
    // �����ӁFLoopLDPC���\�b�h�Ƃ͈قȂ�A�t�@�C������ǂݍ��ޓ_�ɒ��ӂ��邱��
    // rate : [double]�ڑ����̑傫����L�ɑ΂��銄��(L=8��rate=0.5�Ȃ�ڑ����̑傫����4)
    static Matrix<BinaryFiniteField> ConnectIndependentLRLProtographWithProportionalConnection (int l, int r, int L, int copyNum, double rate) {
        return ConnectIndependentLRLProtograph(l, r, L, copyNum, (int)(L * rate));
    }


private :
    typedef struct ___BANDDATA {
        unsigned long top;
        unsigned long bottom;
        unsigned long bandWidth;
    } BandData;

    // �w�肵���ʒu�̑т̍ŏ㕔�E�ŉ����E�c�����擾����֐�
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

    // b-a�̒�����n�����_���Ɏ�����x�N�g���𐶐�
    static vector<int> getRandomSequence(int a, int b, int n) {
        if(a > b) {
            int temp = a;
            a = b;
            b = temp;
        }
        int c = b - a;
        if(n > c) n = c;    // n�͍ő�b-c�܂�
        vector<int> retVec(n);
        RandomValue rv;
        //cout << "bandWidth = " << c << "\n";
        for(int i = 0; i < n; i++) {
            retVec[i] = rv.getRand(c - i) + a;
            //if(retVec[i] < a) cout << "��������\n";
        }
        for(int i = n - 1; i >= 0; i--) {
            for(int j = i - 1; j >= 0; j--) {
                if(retVec[i] >= retVec[j]) {
                    retVec[i]++;
                    //if(retVec[i] >= b) cout << "��������\n";
                }
            }
        }
        std::sort(retVec.begin(), retVec.end());
        return retVec;
    }

    static void _matrixIdentBlockPut(int row, int col, int paramP, Matrix<BinaryFiniteField>* M) {    // �P�ʍs���ݒu
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

    // �u���x�N�g�����擾
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

    // vector�̒l��val���܂܂�Ă��邩�ǂ����𒲂ׁA�Ȃ���Βǉ����\�[�g����
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
    // vector�̒���min�ȏ�max�����̒l���S�Ċ܂܂�Ă���ꍇ��1��Ԃ�
    inline static int containMinToMax(const vector<int>& vec, int min, int max) {
        int counter = 0;
        for(int i = min; i < max; i++) {
            if(contains2(vec, i)) counter++;
        }
        if(counter == max - min) {
            return 1;
        } else return 0;
    }

    // ��W�����擾
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
    class ConstructExeption{};    // ��O�N���X

};

#endif
