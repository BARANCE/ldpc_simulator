/*
    集合演算を行うクラステンプレート
*/

#if !defined(___Class_SetCalc)
#define ___Class_SetCalc

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <functional>

template <class Type> class SetCalc {
public :
    // AまたはB
    static vector<Type> AorB (const vector<Type>* A, const vector<Type>* B) {
        vector<Type> retVec;
        int asize = A->size();
        int bsize = B->size();
        for(int i = 0; i < asize; i++) {
            retVec.push_back((*A)[i]);
        }
        for(int i = 0; i < bsize; i++) {
            int check = 0;
            for(int j = 0; j < asize; j++) {
                if((*B)[i] == (*A)[j]) {
                    check = 1;
                    break;
                }
            }
            if(check == 0) {
                retVec.push_back((*B)[i]);
            }
        }
        return retVec;
    }

    // AかつB
    static vector<Type> AandB (const vector<Type>* A, const vector<Type>* B) {
        vector<Type> retVec;
        int asize = A->size();
        int bsize = B->size();
        for(int i = 0; i < asize; i++) {
            for(int j = 0; j < bsize; j++) {
                if((*A)[i] == (*B)[j]) {
                    retVec.push_back((*A)[i]);
                    break;
                }
            }
        }
        return retVec;
    }

    // 差集合
    static vector<Type> AminusB (const vector<Type>* A, const vector<Type>* B) {
        vector<Type> retVec;
        int asize = A->size();
        int bsize = B->size();
        for(int i = 0; i < asize; i++) {
            int check = 0;
            for(int j = 0; j < bsize; j++) {
                if((*A)[i] == (*B)[j]) {
                    check = 1;
                    break;
                }
            }
            if (check == 0) {
                retVec.push_back((*A)[i]);
            }
        }
        return retVec;
    }

    // 排他的論理和
    static vector<Type> AxorB (const vector<Type>* A, const vector<Type>* B) {
        vector<Type> retVec;
        int asize = A->size();
        int bsize = B->size();
        for(int i = 0; i < asize; i++) {
            int check = 0;
            for(int j = 0; j < bsize; j++) {
                if((*A)[i] == (*B)[j]) {
                    check = 1;
                    break;
                }
            }
            if(check == 0) {
                retVec.push_back((*A)[i]);
            }
        }
        for(int i = 0; i < bsize; i++) {
            int check = 0;
            for(int j = 0; j < asize; j++) {
                if((*A)[j] == (*B)[i]) {
                    check = 1;
                    break;
                }
            }
            if(check == 0) {
                retVec.push_back((*B)[i]);
            }
        }
        return retVec;
    }

    // 特別用途
    static vector<Type> spor (const vector<Type>* A, const vector<Type>* B) {
        vector<Type> retVec;
        vector<Type> tmp;
        retVec = SetCalc::AminusB(B, A);
        tmp = SetCalc::AandB(A, B);
        if(tmp.size() > 0) retVec.push_back(tmp[0]);
        std::sort(retVec.begin(), retVec.end());
        return retVec;
    }
};

#endif
