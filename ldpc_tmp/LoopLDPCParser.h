/*
    LoopLDPC�����̃p�����[�^��ݒ�t�@�C������ǂݎ�邽�߂̃p�[�T
*/

#if !defined(___Class_LoopLDPCParser)
#define ___Class_LoopLDPCParser

#include "LoopLDPC.h"
#include "HashOperator.h"

class LoopLDPCParser {
    vector<LoopLDPC::LBandData> bandData;
public:
    //==== �R���X�g���N�^ ====//
    LoopLDPCParser(const HashOperator& hash);

    //==== ���J���\�b�h ====//
    // Hash����ݒ��ǂݏo���Avector�f�[�^�ɕϊ�����
    vector<LoopLDPC::LBandData> translate();
};

#endif
