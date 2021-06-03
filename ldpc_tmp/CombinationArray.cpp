#include "CombinationArray.h"

vector<int> CombinationArray::next() {
    vector<int> retVec(0);

    vector<int> decisionArray = reverse(countToArray());    // �ǉ�����v�f�̈ʒu�����肷��x�N�g��
    unsigned long dsize = decisionArray.size();
    if(dsize <= size) {
        for(int i = 0; i < dsize; i++) {
            if(decisionArray[i] == 1) {
                retVec.push_back((*vec)[i]);
            }
        }
    } else {
        // FULL
    }
    count++;

    return retVec;
}

//==== �x�N�^�̋t�������߂� ====//
vector<int> CombinationArray::reverse(const vector<int>& v) const {
    unsigned long rsize = v.size();
    vector<int> retVec(rsize);
    for(int i = 0; i <= rsize / 2; i++) {
        retVec[i] = v[rsize - i - 1];
        retVec[rsize - i - 1] = v[i];
    }
    return retVec;
}

//==== �ϐ�count�̒��g���x�N�^�ɕϊ� ====//
vector<int> CombinationArray::countToArray() const {
    unsigned long ketaSize = keta2(count);
    vector<int> retVec(ketaSize);
    for(int i = 0; i < ketaSize; i++) {
        unsigned long val1 = count % pow2(i + 1);
        unsigned long val2 = val1 / pow2(i);
        if(val2 == 1) {
            retVec[ketaSize - i - 1] = 1;
        } else {
            retVec[ketaSize - i - 1] = 0;
        }
    }

    return retVec;
}

//==== a��2�i���ł̌��������߂� ====//
unsigned long CombinationArray::keta2(unsigned long val) const {
    unsigned long c = 0;
    if(val == 0) c++;
    while(val != 0) {
        val /= 2;
        c++;
    }
    return c;
}

//==== 2��val�悵���l���Z�o���� ====//
unsigned long CombinationArray::pow2(unsigned long val) const {
    unsigned long retVal = 1;
    for(int i = 0; i < val; i++) {
        retVal *= 2;
    }
    return retVal;
}
