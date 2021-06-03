#include "Picture.h"

//==== コンストラクタ ====//
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

//==== row行col列目の色を変更 ====//
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

//==== row行col列目の色を取得 ====//
unsigned long Picture::getColor(int row, int col) {
    unsigned long retValue = r[row * x + col];
    retValue *= 0xff;
    retValue += g[row * x + col];
    retValue *= 0xff;
    retValue += b[row * x + col];
    return retValue;
}

//==== bmpを読み込み ====//
void Picture::readBMP(string fileName) {
    int i, j;
    int p;    // カラーパレット番号
    int start_y, stop_y, sign_y;    // 走査方向の違いを記録
    int pad;    // 32bit境界のためのパディング
    int mask, shift;    // 8bit以下の値を取り出す時に利用
    int palsize;    // カラーパレットのサイズ
    int imgsize;    // イメージのサイズ
    int colbit;    // 1画素あたりのビット数
    unsigned char* buf;
    unsigned char* buf_top;
    BMPFILEHEADER bf;    // ファイルヘッダ
    BMPINFOHEADER bi;    // ファイルヘッダ
    RGBQUAD* rgbq;    // カラーパレット

    ifstream fin;
    fin.open(fileName.c_str(), ios::in|ios::binary);
    if(!fin) {
        cout << "ファイルが開けません\n";
        return;
    }

    if(fin.eof()) return;
    fin.read((char*)&bf, sizeof(BMPFILEHEADER));

    fin.close();
}

//==== bmpで書き出し ====//
void Picture::writeBMP(string fileName) {
    int i, j;
    int pad;    // 32bit境界のためのパディング
    int imgsize;    // イメージのサイズ
    unsigned char* buf;
    unsigned char* buf_top;
    BMPFILEHEADER bf;    // ファイルヘッダ
    BMPINFOHEADER bi;    // ファイルヘッダ

    //==== ヘッダ作成 ====//
    pad = (x * 3 + 3) / 4 * 4 - x * 3;    // 32bit境界条件によるパディング
    imgsize = (x * 3 + pad) * y;    // 出力されるデータサイズ

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

    //==== 出力バッファ確保 ====//
    try {
        buf_top = buf = new unsigned char[imgsize];
    } catch(std::bad_alloc) {
        cout << "メモリ確保に失敗しました。\n";
    }

    //==== データ整列 ====//
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

    //==== ヘッダ出力 ====//
    ofstream fout;
    fout.open(fileName.c_str(), ios::out|ios::binary|ios::trunc);
    if(!fout) {
        cout << "ファイルが開けません";
        return;
    }
    fout.write((char*)&bf, sizeof(BMPFILEHEADER));
    fout.write((char*)&bi, sizeof(BMPINFOHEADER));

    //==== 画像データ出力 ====//
    fout.write((char*)buf_top, imgsize);
    fout.close();

    free(buf_top);
}
