#include "HashOperator.h"

//==== コンストラクタ ====//
HashOperator::HashOperator() : keyList(), valueList() {}
HashOperator::HashOperator(const vector<string>& keys, const vector<string>& vals) : keyList(keys.size()), valueList(keys.size()) {
    unsigned long size = keys.size();
    if(size != vals.size()) {
        cout << "[HashOperator.cpp->HashOperator] WARNING:キーと値のベクトルサイズが異なります。ハッシュの長さは入力値のkeysベクトルの長さになるように構成されます。\n";
        cout << "keys.size() :" << keys.size() << ", vals.size() : " << vals.size() << "\n";
    }
    for(unsigned long i = 0; i < size; i++) {    // コピー
        keyList[i] = keys[i];
        valueList[i] = vals[i];
    }
}
HashOperator::HashOperator(const vector<string>& keys, const vector<int>& vals) : keyList(keys.size()), valueList(keys.size()) {
    unsigned long size = keys.size();
    if(size != vals.size()) {
        cout << "[HashOperator.cpp->HashOperator] WARNING:キーと値のベクトルサイズが異なります。ハッシュの長さは入力値のkeysベクトルの長さになるように構成されます。\n";
        cout << "keys.size() :" << keys.size() << ", vals.size() : " << vals.size() << "\n";
    }
    for(unsigned long i = 0; i < size; i++) {    // コピー
        keyList[i] = keys[i];
        stringstream ss;
        string st;
        ss << vals[i];
        ss >> st;
        valueList[i] = st;
    }
}
HashOperator::HashOperator(const vector<string>& keys, const vector<double>& vals) : keyList(keys.size()), valueList(keys.size()) {
    unsigned long size = keys.size();
    if(size != vals.size()) {
        cout << "[HashOperator.cpp->HashOperator] WARNING:キーと値のベクトルサイズが異なります。ハッシュの長さは入力値のkeysベクトルの長さになるように構成されます。\n";
        cout << "keys.size() :" << keys.size() << ", vals.size() : " << vals.size() << "\n";
    }
    for(unsigned long i = 0; i < size; i++) {    // コピー
        keyList[i] = keys[i];
        stringstream ss;
        string st;
        ss << vals[i];
        ss >> st;
        valueList[i] = st;
    }
}
HashOperator::HashOperator(string fileName) {
    readHashFile(fileName);
}

//==== ファイルから読み込み ====//
void HashOperator::readHashFile(string fileName) {
    // ハッシュ初期化
    vector<string> v1, v2;
    keyList = v1, valueList = v2;

    ifstream ifs;
    ifs.open(fileName.c_str(), fstream::in);
    if(ifs.fail()) {
        cout << "[HashOperator.cpp->readHashFile] WARNING:ファイル\"" << fileName << "\"を開けません。ファイルが存在しないか、ファイルまたはディレクトリへの読み取り権限がありません。読み込み処理を中断します。\n";
        return;
    }
    if(ifs.eof()) {
        cout << "[HashOperator.cpp->readHashFile] WARNING:ファイル\"" << fileName << "\"の内容が空です。取得したデータはありません。\n";
        return;
    }
    while(!ifs.eof()) {
        string tmp1;
        ifs >> tmp1;
        vector<string> tmp2 = split(tmp1, ":");
        if(tmp2.size() == 1) {
            cout << "[HashOperator.cpp->readHashFile] WARNING:デリミタが無い行があります。キーと値の組み合わせ以外のデータが記載されています。該当データはキーとして認識され、値は空文字列として保持されます。\n";
            keyList.push_back(tmp2[0]);
            valueList.push_back("");
        } else if(tmp2.size() >= 2) {
            if(tmp2.size() >= 3) {
                cout << "[HashOperator.cpp->readHashFile] WARNING:デリミタが同じ行に2つ以上記述されています。キーと値の組み合わせ以外のデータが記載されています。2つ目までのデータを取得します。\n";
            }
            keyList.push_back(tmp2[0]);
            valueList.push_back(tmp2[1]);
        }
    }
    ifs.close();
}

//==== ファイルへ書き出し ====//
void HashOperator::writeHashFile(string fileName) const {
    ofstream ofs;
    ofs.open(fileName.c_str(), fstream::out);
    if(ofs.fail()) {
        cout << "[HashOperator.cpp->writeHashFile] WARNING:ファイル\"" << fileName << "\"を開けません。使用中であるか、ファイルまたはディレクトリへの書き込み権限がありません。書き込み処理を中断します。\n";
        return;
    }
    unsigned long size = keyList.size();
    for(unsigned long i = 0; i < size; i++) {
        ofs << keyList[i] << ":" << valueList[i] << "\n";
    }
    ofs.close();
}

/*template <> void HashOperator::setValueByKey(string keyName, const string& value) {
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
}*/
