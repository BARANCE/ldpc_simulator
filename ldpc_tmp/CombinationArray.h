/*
    �z��̗v�f�̑g�ݍ��킹��\���N���X
*/

#if !defined (___Class_CombinationArray)
#define ___Class_CombinationArray

#include <vector>
using namespace std;

class CombinationArray {
    const vector<int>* vec;    // �g�ݍ��킹���擾����Ώ̂̃x�N�^
    const unsigned long size;    // �x�N�^�̃T�C�Y
    unsigned long count;
public:
    //==== �R���X�g���N�^ ====//
    CombinationArray(const vector<int>* v) : vec(v), size((*v).size()), count(0) {}

    //==== ���̑g�ݍ��킹�����o�� ====//
    vector<int> next();

    int end();    // �擾�\�ȑg�ݍ��킹���Ō�܂ŒB���Ă����1�A�����łȂ����0��Ԃ�
    vector<int> countToArray() const;
    unsigned long keta2(unsigned long val) const;    // 2�i���ł̌��������߂�
    unsigned long pow2(unsigned long val) const;    // 2��val������߂�
    vector<int> reverse(const vector<int>& v) const;

    unsigned long getCount() const {return count;}
};

#endif
