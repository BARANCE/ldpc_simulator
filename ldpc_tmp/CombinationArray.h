/*
    配列の要素の組み合わせを表すクラス
*/

#if !defined (___Class_CombinationArray)
#define ___Class_CombinationArray

#include <vector>
using namespace std;

class CombinationArray {
    const vector<int>* vec;    // 組み合わせを取得する対称のベクタ
    const unsigned long size;    // ベクタのサイズ
    unsigned long count;
public:
    //==== コンストラクタ ====//
    CombinationArray(const vector<int>* v) : vec(v), size((*v).size()), count(0) {}

    //==== 次の組み合わせを取り出す ====//
    vector<int> next();

    int end();    // 取得可能な組み合わせが最後まで達していれば1、そうでなければ0を返す
    vector<int> countToArray() const;
    unsigned long keta2(unsigned long val) const;    // 2進数での桁数を求める
    unsigned long pow2(unsigned long val) const;    // 2のval乗を求める
    vector<int> reverse(const vector<int>& v) const;

    unsigned long getCount() const {return count;}
};

#endif
