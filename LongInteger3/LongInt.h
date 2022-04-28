#pragma once

#include<vector>
#include"fourier_trans.h"
#include<string>
#include<random>
#include<iostream>
#include<string.h>
#include<sstream>
#include<iomanip>
#include<limits>
#include<thread>

#define LINT_BASE 100000
#define LINT_BASE_DNUM 5
class LongInt
{

private:

	template<class T>
	void swap(T& a, T& b)
	{
		T tmp = std::move(a);
		a = std::move(b);
		b = std::move(tmp);
	}

	std::vector<long long> data;
	bool sign = 0;//���̂Ƃ�1

	void del_zero() & noexcept
	{
		if (data.empty())
		{
			data.push_back(0);
			sign = 0;
			return;
		}
		while (data.size() > 1 && !data.back())
		{
			data.pop_back();
		}
		if (data.size() == 1 && data[0] == 0)
		{
			sign = 0;
		}
	}

	void fix_carry()&
	{
		long long carry = 0;
		for (std::size_t i = 0; i < data.size(); ++i)
		{
			data[i] += carry;
			carry = 0;
			if (data[i] > 0)
			{
				carry = data[i] / LINT_BASE;
				data[i] %= LINT_BASE;
			}
			else if (data[i] < 0)
			{
				carry = (data[i] + 1) / LINT_BASE - 1;
				data[i] -= carry * LINT_BASE;
			}
		}

		while (carry > 0)
		{
			data.push_back(carry % LINT_BASE);
			carry /= LINT_BASE;
		}
		while (carry <= -LINT_BASE)
		{
			data.push_back((carry + 1) % LINT_BASE + LINT_BASE - 1);
			carry = (carry + 1) / LINT_BASE - 1;
		}
		if (carry < 0)data.push_back(carry);
		if (data.back() < 0)
		{
			sign = !sign;
			carry = 1;
			for (std::size_t i = 0; i < data.size(); ++i)
			{
				data[i] = LINT_BASE - 1 - data[i] + carry;
				carry = data[i] / LINT_BASE;
				data[i] %= LINT_BASE;
			}
			data.back() %= LINT_BASE;
		}
		del_zero();
	}
	static void multiply_nr(const std::vector<long long>& num1, const std::vector<long long>& num2, std::vector<long long>& result, std::size_t truncate);
	static void multiply(const std::vector<long long>& num1, const std::vector<long long>& num2, std::vector<long long>& result, std::size_t truncate);
	static void fix_carry_p(std::vector<long long>& num, std::size_t truncate);
public:

	LongInt();//�f�t�H���g�R���X�g���N�^
	LongInt(long long);//�R���X�g���N�^int
	LongInt(std::size_t, bool);
	LongInt(bool, bool);//infinity��Ԃ�

	LongInt(const char*);//�R���X�g���N�^char*
	LongInt(const LongInt&) = default;//�R�s�[�R���X�g���N�^
	LongInt(LongInt&&) = default;//���[�u�R���X�g���N�^

	virtual ~LongInt()noexcept = default;

	LongInt& operator=(const LongInt&) = default;
	LongInt& operator=(LongInt&&) = default;

	bool operator<(const LongInt&)const& noexcept;
	bool operator>(const LongInt&)const& noexcept;
	bool operator<=(const LongInt&)const& noexcept;
	bool operator>=(const LongInt&)const& noexcept;
	bool operator==(const LongInt&)const& noexcept;
	bool operator!=(const LongInt&)const& noexcept;

	LongInt& operator++();
	LongInt operator++(int);
	LongInt& operator--();
	LongInt operator--(int);

	LongInt operator+()const&;
	LongInt operator-()const&;

	LongInt operator+(LongInt)const&;//���Z
	LongInt operator-(LongInt)const&;//���Z
	LongInt operator*(const LongInt&)const&;//��Z
	LongInt operator/(const LongInt&)const&;//���Z
	LongInt operator%(const LongInt&)const&;//���W�������Z

	LongInt& operator+=(const LongInt&);
	LongInt& operator-=(const LongInt&);
	LongInt& operator*=(const LongInt&);
	LongInt& operator/=(const LongInt&);
	LongInt& operator%=(const LongInt&);

	long long operator[](std::size_t)const& noexcept;

	//size���ȉ��̗����𐶐�
	void random(std::size_t size)
	{
		std::size_t i, s = size / LINT_BASE_DNUM + 1;
		static std::random_device seed;
		static std::mt19937_64 mt(seed());
		static std::uniform_int_distribution<> rand(0, LINT_BASE - 1);
		data.resize(s);
		for (i = 0; i < s - 1; ++i)
		{
			data[i] = rand(mt);
		}
		data[s - 1] = rand(mt) % (long long)(std::pow)(10, size % LINT_BASE_DNUM);
		sign = rand(mt) <= LINT_BASE / 2 - 1;
		del_zero();
	}

	static LongInt zero()
	{
		return LongInt();
	}

	//���ɂ���
	void absolute()
	{
		sign = 0;
	}

	bool is_inf()const&
	{
		return (data.size() == 1) && (data[0] == -1);
	}

	bool get_sign()const&
	{
		return sign;
	}
	void tostring(std::string&)const&;
};

namespace std
{
	//�ݏ�v�Z
	LongInt pow(const LongInt&, unsigned long long);

	template<>
	class numeric_limits<LongInt> : public _Num_base
	{
	public:
		static constexpr bool has_infinity = true;
		static constexpr bool is_exact = true;
		static constexpr bool is_integer = true;
		static constexpr bool is_signed = true;
		static constexpr bool is_specialized = true;
		static constexpr bool traps = true;
		static constexpr int digits = 67108864;
		static constexpr int digits10 = 335544320;
		static constexpr int radix = 100000;
		_NODISCARD static LongInt(min)()
		{
			return LongInt(true, true);
		}
		_NODISCARD static LongInt(max)()
		{
			return LongInt(false, true);
		}
		_NODISCARD static LongInt lowest()
		{
			return LongInt(true, false);
		}
		_NODISCARD static LongInt epsilon()
		{
			return LongInt();
		}
		_NODISCARD static LongInt round_error()
		{
			return LongInt();
		}
		_NODISCARD static LongInt denorm_min()
		{
			return LongInt();
		}
		_NODISCARD static LongInt infinity()
		{
			return LongInt(false, true);
		}
		_NODISCARD static LongInt quiet_NaN()
		{
			return LongInt();
		}
		_NODISCARD static LongInt signaling_NaN()
		{
			return LongInt();
		}
	};
}
void product(std::vector<LongInt>&);
//�K��v�Z
LongInt factorial(unsigned int);

std::ostream& operator<<(std::ostream&, const LongInt&);//�W���o��
std::istream& operator>>(std::istream&, LongInt&);//�W������