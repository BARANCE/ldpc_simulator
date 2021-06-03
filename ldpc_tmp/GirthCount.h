#if !defined (___Class_GirthCount)
#define ___Class_GirthCount

#include "Matrix.h"
#include "BinaryFiniteField.h"

/*
    タナーグラフのループの長さを求めるクラス
*/
class GirthCount {
public:
    //static unsigned long loopNum(const Matrix<BinaryFiniteField>& H) {
    static unsigned long loopNum(const SPMatrix& H) {
        unsigned long counter = 0;
        int hr = H.row();
        int hc = H.col();

        for(int i = 0; i < hr; i++) {
            for(int j = 0; j < hc; j++) {
                if((int)H.getValue(i,j) == 1) {
                    for(int k = i + 1; k < hr; k++) {
                        if((int)H.getValue(k,j) == 1) {
                            for(int l = j + 1; l < hc; l++) {
                                if((int)H.getValue(i,l) == 1) {
                                    if((int)H.getValue(k,l) == 1) {
                                        counter++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return counter;
    }
};

#endif
