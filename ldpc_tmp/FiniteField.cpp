#include "FiniteField.h"
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

//===== コンストラクタ =====//
FiniteField::FiniteField() : mod(1) , value(0) , modExist(false), valueExist(false){}
/*FiniteField::FiniteField(const int& p) : mod(p), modExist(true), valueExist(false) {
    if (p < 0) mod = p * (-1);
    if (p == 0) throw OperationErr(this);
}*/
FiniteField::FiniteField(const int& p, const int& val) : mod(p), value(val % p), modExist(true), valueExist(true) {
    if (p < 0) mod = p * (-1);
    if (p == 0) throw OperationErr(this);
    if (val < 0) {
        int tmp = val * (-1);
        FiniteField nfa(p, tmp);
        nfa = nfa.getInverse(OP_ADD);
        value = nfa.getValue();
    }
}
FiniteField::FiniteField(const FiniteField& fa) : modExist(fa.getModExist()), valueExist(fa.getValueExist()) {
    if (modExist) {
        mod = fa.getModular();
        if (valueExist) value = fa.getValue();
    } else if(valueExist) throw UndefErr(this);
}

//===== データ変更 =====//
FiniteField& FiniteField::setModular(const int& m) {
    mod = m;
    modExist = true;
    return (*this);
}
FiniteField& FiniteField::setValue(const int& val) {
    if(!modExist) throw UndefErr(this);
    value = val % mod;
    valueExist = true;
    return (*this);
}

//===== データアクセス =====//
int FiniteField::getModular() const {
    if(!modExist) throw UndefErr(this);
    return mod;
}
int FiniteField::getValue() const {
    //cout << "this\n";
    if(!valueExist) throw UndefErr(this);
    return value;
}
bool FiniteField::getModExist() const {return modExist;}
bool FiniteField::getValueExist() const {return valueExist;}
bool FiniteField::getOperationPropriety() const {
    return (modExist && valueExist);
}
vector<FiniteField> FiniteField::getAllElements() const {
    vector<FiniteField> vec(mod);
    for(int i = 0; i < mod; i++) {
        FiniteField tmp(mod, i);
        vec[i] = tmp;
    }
    return vec;
}
FiniteField FiniteField::getInverse(OperationType op) const {
    if(!getOperationPropriety()) throw UndefErr(this);
    if(op == OP_ADD || op == OP_SUB) {    // 加法逆元
        FiniteField nfa(mod, mod-value);
        return nfa;
    } else if (op == OP_MUL || op == OP_DIV) {    // 乗法逆元(本当は拡張ユークリッドで求める)
        FiniteField nfa;
        nfa.setModular(mod);
        int i;
        for(i = 0; i < mod; i++) {
            int tmp = value * i % mod;
            if(tmp == 1) {
                nfa.setValue(i);
                break;
            }
        }
        if(i == mod) {
            throw OperationErr(this);
        }

        return nfa;
    } else {
        throw UndefErr(this);
        FiniteField nfa(0,0);
        return nfa;
    }
}

string FiniteField::toString() const {
    stringstream ss;
    ss << value;
    return ss.str();
}

//===== 演算子の多重定義 =====//
FiniteField FiniteField::operator=(const FiniteField& a) {
    modExist = a.getModExist();
    valueExist = a.getValueExist();
    if(modExist) mod = a.getModular();
    if(valueExist) value = a.getValue();
    return (*this);
}
FiniteField FiniteField::operator=(const int& a) {
    if(!modExist) throw UndefErr(this);
    int tmp = a;
    bool minus = (tmp < 0) ? true : false;
    if (minus) tmp = tmp * (-1);
    value = tmp % mod;
    valueExist = true;
    if (minus) {
        FiniteField nfa = getInverse(OP_ADD);
        value = nfa.getValue();
    }
    return (*this);
}
FiniteField& FiniteField::operator++() {
    if(!getOperationPropriety()) throw UndefErr(this);
    value = (value + 1) % mod;
    return (*this);
}
FiniteField FiniteField::operator++(int) {
    if(!getOperationPropriety()) throw UndefErr(this);
    FiniteField tmp = (*this);
    value = (value + 1) % mod;
    return tmp;
}
FiniteField& FiniteField::operator--() {
    if(!getOperationPropriety()) throw UndefErr(this);
    FiniteField nfa(mod, 1);
    value = (value + nfa.getInverse(OP_ADD)) % mod;
    return (*this);
}
FiniteField FiniteField::operator--(int) {
    if(!getOperationPropriety()) throw UndefErr(this);
    FiniteField tmp = (*this);
    FiniteField nfa(mod, 1);
    value = (value + nfa.getInverse(OP_ADD)) % mod;
    return tmp;
}

FiniteField operator+(const FiniteField& x, const FiniteField& y){
    if (!x.getOperationPropriety()) throw FiniteField::UndefErr(&x);
    if (!y.getOperationPropriety()) throw FiniteField::UndefErr(&y);
    FiniteField nfa(x.getModular(), x.getValue() + y.getValue());
    return nfa;
}
FiniteField operator-(const FiniteField& x, const FiniteField& y){
    return x + y.getInverse(OP_ADD);
}
FiniteField operator*(const FiniteField& x, const FiniteField& y){
    if (!x.getOperationPropriety()) throw FiniteField::UndefErr(&x);
    if (!y.getOperationPropriety()) throw FiniteField::UndefErr(&y);
    FiniteField nfa(x.getModular(), x.getValue() * y.getValue());
    return nfa;
}
FiniteField operator/(const FiniteField& x, const FiniteField& y){
    return x + y.getInverse(OP_MUL);
}
FiniteField& operator+=(FiniteField& x, const FiniteField& y){
    x = x + y;
    return x;
}
FiniteField& operator-=(FiniteField& x, const FiniteField& y){
    x = x - y;
    return x;
}
FiniteField& operator*=(FiniteField& x, const FiniteField& y){
    x = x * y;
    return x;
}
FiniteField& operator/=(FiniteField& x, const FiniteField& y){
    x = x / y;
    return x;
}
FiniteField operator-(FiniteField& x){    // 単項-演算子
    return x.getInverse(OP_ADD);
}

FiniteField operator+(const FiniteField& x, const int& y){
    if (!x.getOperationPropriety()) throw FiniteField::UndefErr(&x);
    bool minus = (y < 0) ? true : false;    // yが負数かどうか
    int opy = y;
    if(minus) opy *= (-1);    // 正数に直す
    FiniteField nfa(x.getModular(), opy);
    if(minus) nfa = nfa.getInverse(OP_ADD);
    //cout << "nfa = " << nfa << "\n";
    return x + nfa;
}
FiniteField operator-(const FiniteField& x, const int& y){
    return x + (-y);
}
FiniteField operator*(const FiniteField& x, const int& y){
    if (!x.getOperationPropriety()) throw FiniteField::UndefErr(&x);
    bool minus = (y < 0) ? true : false;    // yが負数かどうか
    int opy = y;
    if(minus) opy *= (-1);    // 正数に直す
    FiniteField nfa(x.getModular(), opy);
    if(minus) nfa = nfa.getInverse(OP_ADD);
    //cout << "nfa = " << nfa << "\n";
    return x * nfa;
}
FiniteField operator/(const FiniteField& x, const int& y){
    if (!x.getOperationPropriety()) throw FiniteField::UndefErr(&x);
    bool minus = (y < 0) ? true : false;    // yが負数かどうか
    int opy = y;
    if(minus) opy *= (-1);    // 正数に直す
    FiniteField nfa(x.getModular(), opy);
    if(minus) nfa = nfa.getInverse(OP_ADD);
    //cout << "nfa = " << nfa << "\n";
    return x / nfa;
}

FiniteField operator+(const int& x, const FiniteField& y){
    return y + x;
}
FiniteField operator-(const int& x, const FiniteField& y){
    return y.getInverse(OP_ADD) + x;
}
FiniteField operator*(const int& x, const FiniteField& y){
    return y * x;
}
FiniteField operator/(const int& x, const FiniteField& y){
    if (!y.getOperationPropriety()) throw FiniteField::UndefErr(&y);
    FiniteField nfa(y.getModular(), x);
    return nfa / y;
}

FiniteField& operator+=(FiniteField& x, const int& y){
    x = x + y;
    return x;
}
FiniteField& operator-=(FiniteField& x, const int& y){
    x = x - y;
    return x;
}
FiniteField& operator*=(FiniteField& x, const int& y){
    x = x * y;
    return x;
}
FiniteField& operator/=(FiniteField& x, const int& y){
    x = x / y;
    return x;
}

int& operator+=(int& x, const FiniteField& y){
    x = x + y;
    return x;
}
int& operator-=(int& x, const FiniteField& y){
    x = x - y;
    return x;
}
int& operator*=(int& x, const FiniteField& y){
    x = x * y;
    return x;
}
int& operator/=(int& x, const FiniteField& y){
    x = x / y;
    return x;
}

ostream& operator<<(ostream& ss,const FiniteField& x) {
    return ss << x.toString();
}
fstream& operator>>(fstream& fs, FiniteField& x) {
    if(fs.fail()) throw FiniteField::UndefErr(&x);
    string str;
    fs >> str;
    int val = atoi(str.c_str());
    //cout << "val = " << val << "\n";
    x.setValue(val);
    return fs;
}

//===== 変換コンストラクタ =====//
FiniteField::operator int() { return value;}
FiniteField::operator string() { return toString();}
