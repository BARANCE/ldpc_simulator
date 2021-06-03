/*
�摜�f�[�^��\���N���X
�Q�l�Fhttp://www.mm2d.net/c/
*/

#if !defined(___Class_Picture)
#define ___Class_Picture

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "Matrix.h"
#include "BinaryFiniteField.h"
using namespace std;

//==== BMP�`���t�@�C���ŗ��p����L�q ====//
typedef unsigned char       BYTE;            /* 1byte�����Ȃ����� */
typedef unsigned short      WORD;            /* 2byte�����Ȃ����� */
typedef unsigned long       DWORD;           /* 4byte�����Ȃ����� */
typedef long                LONG;            /* 4byte����         */

#pragma pack(1)
typedef struct BMPFILEHEADER{
  WORD  bfType;
  DWORD bfSize;
  WORD  bfReserved1;
  WORD  bfReserved2;
  DWORD bfOffBits;
}BMPFILEHEADER;

typedef struct BMPINFOHEADER{
  DWORD biSize;
  LONG  biWidth;
  LONG  biHeight;
  WORD  biPlanes;
  WORD  biBitCount;
  DWORD biCompression;
  DWORD biSizeImage;
  LONG  biXPelsPerMeter;
  LONG  biYPelsPerMeter;
  DWORD biClrUsed;
  DWORD biClrImportant;
}BMPINFOHEADER;

typedef struct RGBQUAD{
  BYTE  rgbBlue;
  BYTE  rgbGreen;
  BYTE  rgbRed;
  BYTE  rgbReserved;
}RGBQUAD;

#pragma pack()

class Picture {
    int x;             /* x�����̃s�N�Z���� */
    int y;             /* y�����̃s�N�Z���� */
public:
    //==== �F�萔 ====//
    static const unsigned long BLACK = 0x000000;
    static const unsigned long WHITE = 0xffffff;
    static const unsigned long RED = 0xff0000;
    static const unsigned long GREEN = 0x00ff00;
    static const unsigned long BLUE = 0x0000ff;
    static const unsigned long YELLOW = 0xffff00;
    static const unsigned long CYAN = 0x00ffff;
    static const unsigned long MAGENTA = 0xff00ff;

    //==== �摜�̎��f�[�^ ====//
    vector<unsigned char> r;  /* R�v�f�̋P�x(0�`255) */
    vector<unsigned char> g;  /* G�v�f�̋P�x(0�`255) */
    vector<unsigned char> b;  /* B�v�f�̋P�x(0�`255) */
    vector<unsigned char> a;  /* ���l(�����x)(0�`255) */
    /* �s�N�Z���̃f�[�^�͍�����E�C�ォ�牺�֑����������Ŋi�[ */

    Picture() : x(0), y(0), r(x * y), g(x * y), b(x * y), a(x * y){}
    Picture(int xsize, int ysize) : x(xsize), y(ysize), r(x * y), g(x * y), b(x * y), a(x * y){
        for(int i = 0; i < y; i++) {
            for(int j = 0; j < x; j++) {
                r[i * x + j] = 255;
                g[i * x + j] = 255;
                b[i * x + j] = 255;
                a[i * x + j] = 255;
            }
        }
        //cout << "x * y = " << x * y << "\n";
    }
    Picture(int xsize, int ysize, const vector<unsigned char>& red, const vector<unsigned char>& green, const vector<unsigned char>& blue);
    Picture(const Matrix<BinaryFiniteField>& H);

    //==== row�scol��ڂ̐F���擾 ====//
    unsigned long getColor(int row, int col);

    //==== row�scol��ڂ̐F��ύX ====//
    void setColor(int row, int col, unsigned char red, unsigned char green, unsigned char blue);
    void setColor(int row, int col, unsigned long rgb);

    //==== bmp�œǂݏ��� ====//
    void readBMP(string fileName);
    void writeBMP(string fileName);
};

#endif
