/*
    有限体を定義するクラス
*/
#if !defined (___Class_FiniteField)
#define ___Class_FiniteField

#include <string>
#include <iostream>
#include <vector>
using namespace std;

enum OperationType{OP_ADD, OP_SUB, OP_MUL, OP_DIV};    // 演算の種類を表す列挙体

//===== 有限体クラス =====//
class FiniteField{

    friend class BinaryFiniteField;
    int value;    // Fp上の値
    int mod;    // 上限値
    bool valueExist;    // valueが存在するかどうか
    bool modExist;    // modが存在するかどうか
public :
    explicit FiniteField();
    //FiniteField(const int& p);
    FiniteField(const int& p, const int& val);    // Fp上の値valが格納されたオブジェクトを定義する
    FiniteField(const FiniteField& fa);    // コピーコンストラクタ

    //===== エラークラス =====//
    class UndefErr {    // 未定義値エラー
        const FiniteField* ident;
    public:
        UndefErr(const FiniteField* p) : ident(p) {};
    };
    class OperationErr {    // 計算エラー
        const FiniteField* ident;
    public:
        OperationErr(const FiniteField* p) : ident(p) {};
    };

    //===== データ変更 =====//
    FiniteField& setModular(const int& m);
    FiniteField& setValue(const int& val);    // modが定義されていないとエラーとなる

    //===== データアクセス =====//
    int getModular() const;
    int getValue() const;
    bool getModExist() const;
    bool getValueExist() const;
    bool getOperationPropriety() const;
    vector<FiniteField> getAllElements() const;    // mod元有限体の全ての元を取り出す
    FiniteField getInverse(OperationType op) const;    // 演算opに関する逆元を取り出す
    string toString() const;    // 文字列表現に変更

    //===== 演算子の多重定義 =====//
    FiniteField operator=(const FiniteField& a);
    FiniteField operator=(const int& a);    // modが定義されていないとエラーとなる
    FiniteField& operator++();
    FiniteField operator++(int);
    FiniteField& operator--();
    FiniteField operator--(int);
    friend FiniteField operator+(const FiniteField& x, const FiniteField& y);
    friend FiniteField operator-(const FiniteField& x, const FiniteField& y);
    friend FiniteField operator*(const FiniteField& x, const FiniteField& y);
    friend FiniteField operator/(const FiniteField& x, const FiniteField& y);
    friend FiniteField& operator+=(FiniteField& x, const FiniteField& y);
    friend FiniteField& operator-=(FiniteField& x, const FiniteField& y);
    friend FiniteField& operator*=(FiniteField& x, const FiniteField& y);
    friend FiniteField& operator/=(FiniteField& x, const FiniteField& y);
    friend FiniteField operator-(FiniteField& x);    // 単項-演算子

    friend FiniteField operator+(const FiniteField& x, const int& y);
    friend FiniteField operator-(const FiniteField& x, const int& y);
    friend FiniteField operator*(const FiniteField& x, const int& y);
    friend FiniteField operator/(const FiniteField& x, const int& y);

    friend FiniteField operator+(const int& x, const FiniteField& y);
    friend FiniteField operator-(const int& x, const FiniteField& y);
    friend FiniteField operator*(const int& x, const FiniteField& y);
    friend FiniteField operator/(const int& x, const FiniteField& y);

    friend FiniteField& operator+=(FiniteField& x, const int& y);
    friend FiniteField& operator-=(FiniteField& x, const int& y);
    friend FiniteField& operator*=(FiniteField& x, const int& y);
    friend FiniteField& operator/=(FiniteField& x, const int& y);

    friend int& operator+=(int& x, const FiniteField& y);
    friend int& operator-=(int& x, const FiniteField& y);
    friend int& operator*=(int& x, const FiniteField& y);
    friend int& operator/=(int& x, const FiniteField& y);

    friend ostream& operator<<(ostream& ss,const FiniteField& x);
    friend fstream& operator>>(fstream& fs,FiniteField& x);

    operator int();
    operator string();
};

#endif
