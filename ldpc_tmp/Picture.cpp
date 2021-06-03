#include "Picture.h"

//==== �R���X�g���N�^ ====//
Picture::Picture(const Matrix<BinaryFiniteField>& H) : x(H.col()), y(H.row()), r(x * y), g(x * y), b(x * y), a(x * y) {
    //static const unsigned long int BITCOLOR = 0xff0000;
    for(int i = 0; i < y; i++) {
        for(int j = 0; j < x; j++) {
            if((int)H[i][j] == 1) {
                r[i * x + j] = 255;
                g[i * x + j] = 0;
                b[i * x + j] = 0;
                a[i * x + j] = 255;
            } else {
                r[i * x + j] = 255;
                g[i * x + j] = 255;
                b[i * x + j] = 255;
                a[i * x + j] = 255;
            }
        }
    }
}

//==== row�scol��ڂ̐F��ύX ====//
void Picture::setColor(int row, int col, unsigned char red, unsigned char green, unsigned char blue) {
    r[row * x + col] = red;
    g[row * x + col] = green;
    b[row * x + col] = blue;
}
void Picture::setColor(int row, int col, unsigned long rgb) {
    unsigned char red = rgb / (0x10000);
    unsigned char green = (rgb % (0x10000)) / 0x100;
    unsigned char blue = rgb % 0x100;
    setColor(row, col, red, green, blue);
}

//==== row�scol��ڂ̐F���擾 ====//
unsigned long Picture::getColor(int row, int col) {
    unsigned long retValue = r[row * x + col];
    retValue *= 0xff;
    retValue += g[row * x + col];
    retValue *= 0xff;
    retValue += b[row * x + col];
    return retValue;
}

//==== bmp��ǂݍ��� ====//
void Picture::readBMP(string fileName) {
    int i, j;
    int p;    // �J���[�p���b�g�ԍ�
    int start_y, stop_y, sign_y;    // ���������̈Ⴂ���L�^
    int pad;    // 32bit���E�̂��߂̃p�f�B���O
    int mask, shift;    // 8bit�ȉ��̒l�����o�����ɗ��p
    int palsize;    // �J���[�p���b�g�̃T�C�Y
    int imgsize;    // �C���[�W�̃T�C�Y
    int colbit;    // 1��f������̃r�b�g��
    unsigned char* buf;
    unsigned char* buf_top;
    BMPFILEHEADER bf;    // �t�@�C���w�b�_
    BMPINFOHEADER bi;    // �t�@�C���w�b�_
    RGBQUAD* rgbq;    // �J���[�p���b�g

    ifstream fin;
    fin.open(fileName.c_str(), ios::in|ios::binary);
    if(!fin) {
        cout << "�t�@�C�����J���܂���\n";
        return;
    }

    if(fin.eof()) return;
    fin.read((char*)&bf, sizeof(BMPFILEHEADER));

    fin.close();
}

//==== bmp�ŏ����o�� ====//
void Picture::writeBMP(string fileName) {
    int i, j;
    int pad;    // 32bit���E�̂��߂̃p�f�B���O
    int imgsize;    // �C���[�W�̃T�C�Y
    unsigned char* buf;
    unsigned char* buf_top;
    BMPFILEHEADER bf;    // �t�@�C���w�b�_
    BMPINFOHEADER bi;    // �t�@�C���w�b�_

    //==== �w�b�_�쐬 ====//
    pad = (x * 3 + 3) / 4 * 4 - x * 3;    // 32bit���E�����ɂ��p�f�B���O
    imgsize = (x * 3 + pad) * y;    // �o�͂����f�[�^�T�C�Y

    bf.bfType = *(WORD*)"BM";
    bf.bfSize = imgsize + sizeof(BMPFILEHEADER) + sizeof(BMPINFOHEADER);
    bf.bfReserved1 = 0;
    bf.bfReserved2 = 0;
    bf.bfOffBits = sizeof(BMPFILEHEADER) + sizeof(BMPINFOHEADER);

    bi.biSize = sizeof(BMPINFOHEADER);
    bi.biWidth = x;
    bi.biHeight = y;
    bi.biPlanes = 1;
    bi.biBitCount = 24;
    bi.biCompression = 0;
    bi.biSizeImage = imgsize;
    bi.biXPelsPerMeter = 0;
    bi.biYPelsPerMeter = 0;
    bi.biClrUsed = 0;
    bi.biClrImportant = 0;

    //==== �o�̓o�b�t�@�m�� ====//
    try {
        buf_top = buf = new unsigned char[imgsize];
    } catch(std::bad_alloc) {
        cout << "�������m�ۂɎ��s���܂����B\n";
    }

    //==== �f�[�^���� ====//
    char tmp;
    for(int i = y - 1; i >= 0; i--) {
        for(j = 0; j < x; j++) {
            //cout << j + i * x << ", " << b[j + i * x] << "\n";
            *(buf++) = b[j + i * x];
            //tmp = b[j + i * x];
            *(buf++) = g[j + i * x];
            //tmp = g[j + i * x];
            *(buf++) = r[j + i * x];
            //tmp = r[j + i * x];
        }
        for(j = 0; j < pad; j++) {
            *(buf++) = 0;
            //tmp = 0;
        }
        //cout << tmp << " ";
    }

    //==== �w�b�_�o�� ====//
    ofstream fout;
    fout.open(fileName.c_str(), ios::out|ios::binary|ios::trunc);
    if(!fout) {
        cout << "�t�@�C�����J���܂���";
        return;
    }
    fout.write((char*)&bf, sizeof(BMPFILEHEADER));
    fout.write((char*)&bi, sizeof(BMPINFOHEADER));

    //==== �摜�f�[�^�o�� ====//
    fout.write((char*)buf_top, imgsize);
    fout.close();

    free(buf_top);
}
