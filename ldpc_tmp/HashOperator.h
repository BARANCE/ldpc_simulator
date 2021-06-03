/*
    �n�b�V���\�����`����N���X
*/

#if !defined (___Class_HashOperator)
#define ___Class_HashOperator

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

class HashOperator {
    vector<string> keyList;    // �L�[
    vector<string> valueList;    // �l
public :
    //==== �R���X�g���N�^ ====//
    HashOperator();
    HashOperator(const vector<string>& keys, const vector<string>& vals);
    HashOperator(const vector<string>& keys, const vector<int>& vals);
    HashOperator(const vector<string>& keys, const vector<double>& vals);
    HashOperator(string fileName);

    //==== �l�����݂��邩�𒲂ׂ� ====//
    bool tryKey(string keyName) const {
        unsigned long size = keyList.size();
        unsigned long i = 0;
        for(i = 0; i < size; i++) {
            if(keyList[i] == keyName) return true;
        }
        return false;
    }

    //==== �l�̓ǂݎ�� ====//
    template <class Type> Type getValueByKey(string keyName) const {
        unsigned long size = keyList.size();
        unsigned long i = 0;
        for(i = 0; i < size; i++) {
            if(keyList[i] == keyName) break;
        }
        if(i == size) {
            cout << "[HashOperator.cpp->getValueByKey] WARNING:�w�肳�ꂽ�L�[\"" << keyName << "\"���n�b�V���ɑ��݂��܂���B�����𒆒f���܂��B\n";
            Type t = 0;
            return t;
        }
        stringstream ss;
        Type t;
        ss << valueList[i];
        ss >> t;
        return t;
    }

    //==== �l�̐ݒ� ====//
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

	//==== template�֘A��debug�p ====//
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
        cout << "[HashOperator.cpp->getValueByKey] WARNING:�w�肳�ꂽ�L�[\"" << keyName << "\"���n�b�V���ɑ��݂��܂���B�����𒆒f���܂��B\n";
        string t;
        return t;
    }
    return valueList[i];
}
	//==== debug�����܂� ====//

    //==== �t�@�C������ǂݍ��� ====//
    void readHashFile(string fileName);    // ���݂̒l�͏㏑������邱�Ƃɒ���

    //==== �t�@�C���֏����o�� ====//
    void writeHashFile(string fileName) const;

    //==== �����񂩂番�� ====//
    // ���C�u�����Fhttp://goodjob.boy.jp/chirashinoura/id/100.html
    static vector<string> split(const string& str, const string& delim) {
        vector<string> result;    // �ԋp�x�N�g��
        int cutAt;
        string ostr = str;    // �R�s�[
        while( (cutAt = ostr.find_first_of(delim)) != ostr.npos ) {    // ������delim��������Ȃ��Ȃ�܂ŌJ��Ԃ��Bdelim�̈ʒu��cutAt�Ɋi�[
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

    //==== ������ <-> �^�ϊ��e���v���[�g ====//
    template <class Type> string typeToString(const Type* f) const {    // ����̌^��������
        string retVal;
        stringstream ss;
        ss << (*f);
        ss >> retVal;
        return retVal;
    }
    template <class Type> Type stringToType(const string* s, const Type* ident) const {    // �����񁨓���̌^(ident�͌^���ʗp)
        Type retVal;
        stringstream ss;
        ss << (*s);
        ss >> retVal;
        return retVal;
    }

    //==== ������ϊ� ====//
    stringstream& getStringstream() const {    // ���֐��Ŏg�p���镶����X�g���[���擾�p�֐�
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
