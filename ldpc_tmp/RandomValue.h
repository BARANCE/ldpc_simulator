/*
    �����������N���X
*/

#if !defined(___Class_RandomValue)
#define ___Class_RandomValue

#include <cstdlib>
#include <ctime>
#include "mt19937ar-cok.h"

class RandomValue {
public :
    //===== �R���X�g���N�^ =====//
    RandomValue(){
        //srand((unsigned)time(NULL));    // ������
    }

    //===== �����������\�b�h =====//
    // min�ȏ�max�����̗����𐶐�����
    int getRand(int min, int max) {
        double r = genrand_real2();
        while (r >= (double)1) r = genrand_real2();
        r = r * (max - min);
        r += min;

        int tmp = (int)r;
        return tmp;
    }
    unsigned long getRand(unsigned long min, unsigned long max) {
        double r = genrand_real2();
        while (r >= (double)1) r = genrand_real2();
        r = r * (max - min);
        r += min;

        unsigned long tmp = (unsigned long)r;
        return tmp;
    }
    // 0�ȏ�max�����̗����𐶐�����
    int getRand(int max) {
        return getRand(0, max);
    }
    unsigned long getRand(unsigned long max) {
        return getRand(0UL, max);
    }
};

#endif
