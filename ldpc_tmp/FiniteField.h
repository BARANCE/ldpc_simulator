/*
    �L���̂��`����N���X
*/
#if !defined (___Class_FiniteField)
#define ___Class_FiniteField

#include <string>
#include <iostream>
#include <vector>
using namespace std;

enum OperationType{OP_ADD, OP_SUB, OP_MUL, OP_DIV};    // ���Z�̎�ނ�\���񋓑�

//===== �L���̃N���X =====//
class FiniteField{

    friend class BinaryFiniteField;
    int value;    // Fp��̒l
    int mod;    // ����l
    bool valueExist;    // value�����݂��邩�ǂ���
    bool modExist;    // mod�����݂��邩�ǂ���
public :
    explicit FiniteField();
    //FiniteField(const int& p);
    FiniteField(const int& p, const int& val);    // Fp��̒lval���i�[���ꂽ�I�u�W�F�N�g���`����
    FiniteField(const FiniteField& fa);    // �R�s�[�R���X�g���N�^

    //===== �G���[�N���X =====//
    class UndefErr {    // ����`�l�G���[
        const FiniteField* ident;
    public:
        UndefErr(const FiniteField* p) : ident(p) {};
    };
    class OperationErr {    // �v�Z�G���[
        const FiniteField* ident;
    public:
        OperationErr(const FiniteField* p) : ident(p) {};
    };

    //===== �f�[�^�ύX =====//
    FiniteField& setModular(const int& m);
    FiniteField& setValue(const int& val);    // mod����`����Ă��Ȃ��ƃG���[�ƂȂ�

    //===== �f�[�^�A�N�Z�X =====//
    int getModular() const;
    int getValue() const;
    bool getModExist() const;
    bool getValueExist() const;
    bool getOperationPropriety() const;
    vector<FiniteField> getAllElements() const;    // mod���L���̂̑S�Ă̌������o��
    FiniteField getInverse(OperationType op) const;    // ���Zop�Ɋւ���t�������o��
    string toString() const;    // ������\���ɕύX

    //===== ���Z�q�̑��d��` =====//
    FiniteField operator=(const FiniteField& a);
    FiniteField operator=(const int& a);    // mod����`����Ă��Ȃ��ƃG���[�ƂȂ�
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
    friend FiniteField operator-(FiniteField& x);    // �P��-���Z�q

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
