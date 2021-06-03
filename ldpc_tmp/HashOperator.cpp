#include "HashOperator.h"

//==== �R���X�g���N�^ ====//
HashOperator::HashOperator() : keyList(), valueList() {}
HashOperator::HashOperator(const vector<string>& keys, const vector<string>& vals) : keyList(keys.size()), valueList(keys.size()) {
    unsigned long size = keys.size();
    if(size != vals.size()) {
        cout << "[HashOperator.cpp->HashOperator] WARNING:�L�[�ƒl�̃x�N�g���T�C�Y���قȂ�܂��B�n�b�V���̒����͓��͒l��keys�x�N�g���̒����ɂȂ�悤�ɍ\������܂��B\n";
        cout << "keys.size() :" << keys.size() << ", vals.size() : " << vals.size() << "\n";
    }
    for(unsigned long i = 0; i < size; i++) {    // �R�s�[
        keyList[i] = keys[i];
        valueList[i] = vals[i];
    }
}
HashOperator::HashOperator(const vector<string>& keys, const vector<int>& vals) : keyList(keys.size()), valueList(keys.size()) {
    unsigned long size = keys.size();
    if(size != vals.size()) {
        cout << "[HashOperator.cpp->HashOperator] WARNING:�L�[�ƒl�̃x�N�g���T�C�Y���قȂ�܂��B�n�b�V���̒����͓��͒l��keys�x�N�g���̒����ɂȂ�悤�ɍ\������܂��B\n";
        cout << "keys.size() :" << keys.size() << ", vals.size() : " << vals.size() << "\n";
    }
    for(unsigned long i = 0; i < size; i++) {    // �R�s�[
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
        cout << "[HashOperator.cpp->HashOperator] WARNING:�L�[�ƒl�̃x�N�g���T�C�Y���قȂ�܂��B�n�b�V���̒����͓��͒l��keys�x�N�g���̒����ɂȂ�悤�ɍ\������܂��B\n";
        cout << "keys.size() :" << keys.size() << ", vals.size() : " << vals.size() << "\n";
    }
    for(unsigned long i = 0; i < size; i++) {    // �R�s�[
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

//==== �t�@�C������ǂݍ��� ====//
void HashOperator::readHashFile(string fileName) {
    // �n�b�V��������
    vector<string> v1, v2;
    keyList = v1, valueList = v2;

    ifstream ifs;
    ifs.open(fileName.c_str(), fstream::in);
    if(ifs.fail()) {
        cout << "[HashOperator.cpp->readHashFile] WARNING:�t�@�C��\"" << fileName << "\"���J���܂���B�t�@�C�������݂��Ȃ����A�t�@�C���܂��̓f�B���N�g���ւ̓ǂݎ�茠��������܂���B�ǂݍ��ݏ����𒆒f���܂��B\n";
        return;
    }
    if(ifs.eof()) {
        cout << "[HashOperator.cpp->readHashFile] WARNING:�t�@�C��\"" << fileName << "\"�̓��e����ł��B�擾�����f�[�^�͂���܂���B\n";
        return;
    }
    while(!ifs.eof()) {
        string tmp1;
        ifs >> tmp1;
        vector<string> tmp2 = split(tmp1, ":");
        if(tmp2.size() == 1) {
            cout << "[HashOperator.cpp->readHashFile] WARNING:�f���~�^�������s������܂��B�L�[�ƒl�̑g�ݍ��킹�ȊO�̃f�[�^���L�ڂ���Ă��܂��B�Y���f�[�^�̓L�[�Ƃ��ĔF������A�l�͋󕶎���Ƃ��ĕێ�����܂��B\n";
            keyList.push_back(tmp2[0]);
            valueList.push_back("");
        } else if(tmp2.size() >= 2) {
            if(tmp2.size() >= 3) {
                cout << "[HashOperator.cpp->readHashFile] WARNING:�f���~�^�������s��2�ȏ�L�q����Ă��܂��B�L�[�ƒl�̑g�ݍ��킹�ȊO�̃f�[�^���L�ڂ���Ă��܂��B2�ڂ܂ł̃f�[�^���擾���܂��B\n";
            }
            keyList.push_back(tmp2[0]);
            valueList.push_back(tmp2[1]);
        }
    }
    ifs.close();
}

//==== �t�@�C���֏����o�� ====//
void HashOperator::writeHashFile(string fileName) const {
    ofstream ofs;
    ofs.open(fileName.c_str(), fstream::out);
    if(ofs.fail()) {
        cout << "[HashOperator.cpp->writeHashFile] WARNING:�t�@�C��\"" << fileName << "\"���J���܂���B�g�p���ł��邩�A�t�@�C���܂��̓f�B���N�g���ւ̏������݌���������܂���B�������ݏ����𒆒f���܂��B\n";
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
        cout << "[HashOperator.cpp->getValueByKey] WARNING:�w�肳�ꂽ�L�[\"" << keyName << "\"���n�b�V���ɑ��݂��܂���B�����𒆒f���܂��B\n";
        string t;
        return t;
    }
    return valueList[i];
}*/
