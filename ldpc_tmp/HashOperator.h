/*
    ハッシュ構造を定義するクラス
*/

#if !defined (___Class_HashOperator)
#define ___Class_HashOperator

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

class HashOperator {
    vector<string> keyList;    // キー
    vector<string> valueList;    // 値
public :
    //==== コンストラクタ ====//
    HashOperator();
    HashOperator(const vector<string>& keys, const vector<string>& vals);
    HashOperator(const vector<string>& keys, const vector<int>& vals);
    HashOperator(const vector<string>& keys, const vector<double>& vals);
    HashOperator(string fileName);

    //==== 値が存在するかを調べる ====//
    bool tryKey(string keyName) const {
        unsigned long size = keyList.size();
        unsigned long i = 0;
        for(i = 0; i < size; i++) {
            if(keyList[i] == keyName) return true;
        }
        return false;
    }

    //==== 値の読み取り ====//
    template <class Type> Type getValueByKey(string keyName) const {
        unsigned long size = keyList.size();
        unsigned long i = 0;
        for(i = 0; i < size; i++) {
            if(keyList[i] == keyName) break;
        }
        if(i == size) {
            cout << "[HashOperator.cpp->getValueByKey] WARNING:指定されたキー\"" << keyName << "\"がハッシュに存在しません。処理を中断します。\n";
            Type t = 0;
            return t;
        }
        stringstream ss;
        Type t;
        ss << valueList[i];
        ss >> t;
        return t;
    }

    //==== 値の設定 ====//
    template <class Type> void setValueByKey(string keyName, const Type& value) {
        unsigned long size = keyList.size();
        unsigned long i = 0;
        for(i = 0; i < size; i++) {
            if(keyList[i] == keyName) break;
        }

        stringstream ss;
        string tmp;
        ss << value;
        ss >> tmp;

        if(i == size) {
            keyList.push_back(keyName);
            valueList.push_back(tmp);
        } else {
            valueList[i] = tmp;
        }
    }

	//==== template関連のdebug用 ====//
	template <> void HashOperator::setValueByKey(string keyName, const string& value) {
    unsigned long size = keyList.size();
    unsigned long i = 0;
    for(i = 0; i < size; i++) {
        if(keyList[i] == keyName) break;
    }

    if(i == size) {
        keyList.push_back(keyName);
        valueList.push_back(value);
    } else {
        valueList[i] = value;
    }
}

template <> string HashOperator::getValueByKey<string>(string keyName) const {
    unsigned long size = keyList.size();
    unsigned long i = 0;
    for(i = 0; i < size; i++) {
        if(keyList[i] == keyName) break;
    }
    if(i == size) {
        cout << "[HashOperator.cpp->getValueByKey] WARNING:指定されたキー\"" << keyName << "\"がハッシュに存在しません。処理を中断します。\n";
        string t;
        return t;
    }
    return valueList[i];
}
	//==== debugここまで ====//

    //==== ファイルから読み込み ====//
    void readHashFile(string fileName);    // 現在の値は上書きされることに注意

    //==== ファイルへ書き出し ====//
    void writeHashFile(string fileName) const;

    //==== 文字列から分割 ====//
    // ライブラリ：http://goodjob.boy.jp/chirashinoura/id/100.html
    static vector<string> split(const string& str, const string& delim) {
        vector<string> result;    // 返却ベクトル
        int cutAt;
        string ostr = str;    // コピー
        while( (cutAt = ostr.find_first_of(delim)) != ostr.npos ) {    // 文字列delimが見つからなくなるまで繰り返す。delimの位置をcutAtに格納
            if(cutAt > 0) {
                result.push_back(ostr.substr(0, cutAt));
            }
            ostr = ostr.substr(cutAt + 1);
        }
        if(ostr.size() > 0) {
            result.push_back(ostr);
        }
        return result;
    }

    //==== 文字列 <-> 型変換テンプレート ====//
    template <class Type> string typeToString(const Type* f) const {    // 特定の型→文字列
        string retVal;
        stringstream ss;
        ss << (*f);
        ss >> retVal;
        return retVal;
    }
    template <class Type> Type stringToType(const string* s, const Type* ident) const {    // 文字列→特定の型(identは型識別用)
        Type retVal;
        stringstream ss;
        ss << (*s);
        ss >> retVal;
        return retVal;
    }

    //==== 文字列変換 ====//
    stringstream& getStringstream() const {    // 他関数で使用する文字列ストリーム取得用関数
        stringstream ss;
        for(unsigned long i = 0; i < keyList.size(); i++) {
            ss << keyList[i] << ":" << valueList[i] << "\n";
        }
        return ss;
    }
    string toString() const {return getStringstream().str();}
    void print() const {cout << toString();}
    friend ostream& operator<<(ostream& ss, const HashOperator& v) {return ss << v.toString();}
};

#endif
