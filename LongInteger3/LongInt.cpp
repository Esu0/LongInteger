//#define _MEMO
#include"LongInt.h"

LongInt::LongInt() :data(1), sign(0)
{}

LongInt::LongInt(long long num) : data(1, num), sign(0)
{
	fix_carry();
}

LongInt::LongInt(const char* snum) : data(), sign(0)
{
	int tmp = 0, m = 1;
	std::size_t i = 0;
	std::size_t size = strlen(snum);

	if (*snum == '-')
	{
		++i;
		sign = 1;
	}
	if (*snum == '+')
	{
		++i;
	}
	for (std::size_t j = size; j > i;)
	{
		--j;
		tmp += m * (snum[j] - '0');
		m *= 10;
		if ((size - j) % LINT_BASE_DNUM == 0)
		{
			data.push_back(tmp);
			tmp = 0;
			m = 1;
		}
	}
	if (tmp != 0)data.push_back(tmp);
	del_zero();
}

LongInt::LongInt(std::size_t size, bool) :data(size)
{}

LongInt::LongInt(bool _sig, bool): data(1,-1), sign(_sig)
{}

long long LongInt::operator[](std::size_t index)const& noexcept
{
	return index < data.size() ? data[index] : 0;
}

void LongInt::tostring(std::string& str)const&
{
	if (is_inf())
	{
		str = sign ? "-inf" : "inf";
		return;
	}
	if (sign)str.push_back('-');
	std::ostringstream oss;
	oss << data.back();
	for (unsigned long long i = data.size() - 1; i > 0;)
	{
		--i;
		oss << std::setfill('0') << std::setw(5) << data[i];
	}
	str += oss.str();
}

bool LongInt::operator<(const LongInt& num2)const& noexcept
{
	if (sign == num2.sign)
	{
		bool _inf1 = is_inf(), _inf2 = num2.is_inf();
		if (_inf1 && _inf2)return false;
		else if (_inf1)return sign;
		else if (_inf2)return !sign;
		else if (data.size() == num2.data.size())
		{
			for (std::size_t i = data.size(); i > 0;)
			{
				--i;
				if (data[i] < num2.data[i])return !sign;
				if (data[i] > num2.data[i])return sign;
			}
			return false;
		}
		else if (data.size() < num2.data.size())return !sign;
		else return sign;
	}
	else return sign;
}

bool LongInt::operator>(const LongInt& num2)const& noexcept
{
	if (sign == num2.sign)
	{
		bool _inf1 = is_inf(), _inf2 = num2.is_inf();
		if (_inf1 && _inf2)return false;
		else if (_inf1)return !sign;
		else if (_inf2)return sign;
		else if (data.size() == num2.data.size())
		{
			for (std::size_t i = data.size(); i > 0;)
			{
				--i;
				if (data[i] > num2.data[i])return !sign;
				if (data[i] < num2.data[i])return sign;
			}
			return false;
		}
		else if (data.size() > num2.data.size())return !sign;
		else return sign;
	}
	else return num2.sign;
}

bool LongInt::operator>=(const LongInt& num2)const& noexcept
{
	return !operator<(num2);
}

bool LongInt::operator<=(const LongInt& num2)const& noexcept
{
	return !operator>(num2);
}

bool LongInt::operator==(const LongInt& num2)const& noexcept
{
	if (sign == num2.sign && data.size() == num2.data.size())
	{
		for (std::size_t i = 0; i < data.size(); ++i)
		{
			if (data[i] != num2.data[i])return false;
		}
		return true;
	}
	else return false;
}

bool LongInt::operator!=(const LongInt& num2)const& noexcept
{
	return !operator==(num2);
}

LongInt& LongInt::operator++()&
{
	if (is_inf())return *this;
	if (sign)
	{
		if (data.size() == 1 && data[0] == 1)
		{
			data[0] = 0;
			sign = 0;
			return *this;
		}
		for (std::size_t i = 0; i < data.size(); ++i)
		{
			if (!(data[i]))
			{
				data[i] = LINT_BASE - 1;
			}
			else
			{
				--(data[i]);
				return *this;
			}
		}
	}
	else
	{
		for (std::size_t i = 0; i < data.size(); ++i)
		{
			++(data[i]);
			if (data[i] < LINT_BASE)return *this;
			else data[i] %= LINT_BASE;
		}
		data.push_back(1);
	}
	return *this;
}

LongInt LongInt::operator++(int)&
{
	LongInt tmp = *this;
	this->operator++();
	return tmp;
}

LongInt& LongInt::operator--()&
{
	if (is_inf())return *this;
	if (sign)
	{
		for (std::size_t i = 0; i < data.size(); ++i)
		{
			++(data[i]);
			if (data[i] < LINT_BASE)return *this;
			else data[i] %= LINT_BASE;
		}
		data.push_back(1);
	}
	else
	{
		if (data.size() == 1 && data[0] == 0)
		{
			data[0] = 1;
			sign = 1;
			return *this;
		}
		for (std::size_t i = 0; i < data.size(); ++i)
		{
			if (!(data[i]))
			{
				data[i] = LINT_BASE - 1;
			}
			else
			{
				--(data[i]);
				return *this;
			}
		}
	}
	return *this;
}

LongInt LongInt::operator--(int)&
{
	LongInt tmp = *this;
	this->operator--();
	return tmp;
}

LongInt LongInt::operator+()const&
{
	return *this;
}

LongInt LongInt::operator-()const&
{
	LongInt tmp = *this;
	if (!(tmp.data.size() == 1 && tmp.data[0] == 0))tmp.sign = !tmp.sign;
	return tmp;
}

LongInt LongInt::operator+(LongInt num)const&
{
	{
		bool _inf1 = is_inf(), _inf2 = num.is_inf();
		if (_inf1 || _inf2)
		{
			if (sign != num.sign && _inf1 && _inf2)return LongInt();
			else return _inf1 ? *this : num;
		}
	}
	std::size_t i;
	if (num.sign != sign)
	{
		for (i = 0; i < num.data.size(); ++i)num.data[i] -= (*this)[i];
		for (; i < data.size(); ++i)num.data.push_back(-(this->data[i]));
	}
	else
	{
		for (i = 0; i < num.data.size(); ++i)num.data[i] += (*this)[i];
		for (; i < data.size(); ++i)num.data.push_back(this->data[i]);
	}
	num.fix_carry();
	return num;
}

LongInt& LongInt::operator+=(const LongInt& num)&
{
	{
		bool _inf1 = is_inf(), _inf2 = num.is_inf();
		if (_inf1 || _inf2)
		{
			if (sign == num.sign && _inf1 && _inf2)*this = LongInt();
			else if(!_inf1 && _inf2)*this = num;
			return *this;
		}
	}
	unsigned long long i;
	if (sign == num.sign)
	{
		for (i = 0; i < (std::min)(num.data.size(), data.size()); ++i)data[i] += num.data[i];
		for (; i < num.data.size(); ++i)data.push_back(num.data[i]);
	}
	else
	{
		for (i = 0; i < (std::min)(data.size(), num.data.size()); ++i)data[i] -= num.data[i];
		for (; i < num.data.size(); ++i)data.push_back(-num.data[i]);
	}
	fix_carry();
	return *this;
}

LongInt LongInt::operator-(LongInt num)const&
{
	{
		bool _inf1 = is_inf(), _inf2 = num.is_inf();
		if (_inf1 || _inf2)
		{
			if (sign == num.sign && _inf1 && _inf2)return LongInt();
			else return _inf1 ? *this : -num;
		}
	}
	std::size_t i;
	num.sign = !num.sign;
	if (sign != num.sign)
	{
		for (i = 0; i < num.data.size(); ++i)num.data[i] -= (*this)[i];
		for (; i < data.size(); ++i)num.data.push_back(-(this->data[i]));
	}
	else
	{
		for (i = 0; i < num.data.size(); ++i)num.data[i] += (*this)[i];
		for (; i < data.size(); ++i)num.data.push_back(this->data[i]);
	}
	num.fix_carry();
	return num;
}

LongInt& LongInt::operator-=(const LongInt& num)&
{
	{
		bool _inf1 = is_inf(), _inf2 = num.is_inf();
		if (_inf1 || _inf2)
		{
			if (sign != num.sign && _inf1 && _inf2)*this = LongInt();
			else if (!_inf1 && _inf2)*this = -num;
			return *this;
		}
	}
	unsigned long long i;
	if (sign != num.sign)
	{
		for (i = 0; i < (std::min)(num.data.size(), data.size()); ++i)data[i] += num.data[i];
		for (; i < num.data.size(); ++i)data.push_back(num.data[i]);
	}
	else
	{
		for (i = 0; i < (std::min)(num.data.size(), data.size()); ++i)data[i] -= num.data[i];
		for (; i < num.data.size(); ++i)data.push_back(-num.data[i]);
	}
	fix_carry();
	return *this;
}



#define NTT_PRIME1 2013265921
#define NTT_ROOT1 137
#define NTT_PRIME2 1811939329
#define NTT_ROOT2 136

//オーバーフローあり
void LongInt::fix_carry_p(std::vector<long long>& num, std::size_t truncate = 0)
{
	long long carry = 0;
	std::size_t i;
	for (i = 0; i < truncate; i++)
	{
		num[i] += carry;
		carry = num[i] / LINT_BASE;
	}
	for (; i < num.size(); ++i)
	{
		num[i] += carry;
		carry = num[i] / LINT_BASE;
		num[i - truncate] = num[i] % LINT_BASE;
	}
	for (i = 0; i < truncate; ++i)num.pop_back();
}

//resultはメモリ確保の必要無し
void LongInt::multiply_nr(const std::vector<long long>& num1, const std::vector<long long>& num2, std::vector<long long>& result, std::size_t truncate = 0)
{
	result.resize(num1.size() + num2.size());
	for (unsigned long long i = 0; i < num1.size(); ++i)
	{
		for (unsigned long long j = 0; j < num2.size(); ++j)
		{
			result[i + j] += num1[i] * num2[j];
		}
	}
	LongInt::fix_carry_p(result, truncate);
}

//resultはメモリ確保の必要無し
void LongInt::multiply(const std::vector<long long>& num1, const std::vector<long long>& num2, std::vector<long long>& result, std::size_t truncate = 0)
{
	unsigned long long size = 1;
	int cnt = 0;
	unsigned long long i;
	{
		std::size_t _Size = (std::max)(num1.size(), num2.size());
		while (size < _Size)
		{
			size <<= 1;
			cnt += 1;
		}
	}
	if ((unsigned long long)(cnt * 20 + 2) > num2.size() || (unsigned long long)(cnt * 20 + 2) > num1.size())
	{
		LongInt::multiply_nr(num1, num2, result);
	}
	else
	{
		static discrete::Number_Theoretic_Transform ntt1(NTT_PRIME1, NTT_ROOT1), ntt2(NTT_PRIME2, NTT_ROOT2);
		std::vector<std::vector<long long>> tmp1(1, num1);
		std::vector<std::vector<long long>> tmp2(1, num2);
		tmp1[0].resize(size);
		tmp2[0].resize(size);
		tmp1.resize(4, tmp1[0]);
		tmp2.resize(4, tmp2[0]);
		long long _w1 = ntt1.getRoot(cnt), _w2 = ntt2.getRoot(cnt);
		long long _mul = 1;
		for (i = 0; i < size; ++i)
		{
			tmp1[2][i] *= _mul;
			tmp1[2][i] %= NTT_PRIME1;
			tmp2[2][i] *= _mul;
			tmp2[2][i] %= NTT_PRIME1;
			_mul *= _w1;
			_mul %= NTT_PRIME1;
		}
		_mul = 1;
		for (i = 0; i < size; ++i)
		{
			tmp1[3][i] *= _mul;
			tmp1[3][i] %= NTT_PRIME2;
			tmp2[3][i] *= _mul;
			tmp2[3][i] %= NTT_PRIME2;
			_mul *= _w2;
			_mul %= NTT_PRIME2;
		}
		long long size_inv1 = number::minv((long long)(size * 2), NTT_PRIME1);
		long long size_inv2 = number::minv((long long)(size * 2), NTT_PRIME2);
		ntt1.trans_f(tmp1[0], tmp1[2], tmp2[0], tmp2[2], cnt - 1, 0, size - 1);
		ntt2.trans_f(tmp1[1], tmp1[3], tmp2[1], tmp2[3], cnt - 1, 0, size - 1);
		for (i = 0; i < size; ++i)
		{
			tmp1[0][i] *= tmp2[0][i];
			tmp1[0][i] %= NTT_PRIME1;
			tmp1[2][i] *= tmp2[2][i];
			tmp1[2][i] %= NTT_PRIME1;
			tmp1[1][i] *= tmp2[1][i];
			tmp1[1][i] %= NTT_PRIME2;
			tmp1[3][i] *= tmp2[3][i];
			tmp1[3][i] %= NTT_PRIME2;
		}
		ntt1.itrans_t(tmp1[0], tmp1[2], cnt - 1, 0, size - 1);
		ntt2.itrans_t(tmp1[1], tmp1[3], cnt - 1, 0, size - 1);
		result.resize((std::max)(num1.size() + num2.size(), size));
		_w1 = ntt1.getiRoot(cnt);
		_w2 = ntt2.getiRoot(cnt);
		_mul = 1;
		long long _mul2 = 1;
		for (i = 0; i < size; ++i)
		{
			tmp1[2][i] *= _mul; tmp1[2][i] %= NTT_PRIME1;
			tmp1[3][i] *= _mul2; tmp1[3][i] %= NTT_PRIME2;
			long long _t1 = tmp1[0][i] + tmp1[2][i], _t2 = tmp1[1][i] + tmp1[3][i];
			result[i] = number::garner((_t1 * size_inv1) % NTT_PRIME1, (_t2 * size_inv2) % NTT_PRIME2, NTT_PRIME1, NTT_PRIME2);
			_mul *= _w1; _mul %= NTT_PRIME1;
			_mul2 *= _w2; _mul2 %= NTT_PRIME2;
		}
		for (; i < result.size(); ++i)
		{
			unsigned long long _ti = i - size;
			long long _t1 = tmp1[0][_ti] - tmp1[2][_ti], _t2 = tmp1[1][_ti] - tmp1[3][_ti];
			if (_t1 < 0)_t1 += NTT_PRIME1;
			if (_t2 < 0)_t2 += NTT_PRIME2;
			result[i] = number::garner((_t1 * size_inv1) % NTT_PRIME1, (_t2 * size_inv2) % NTT_PRIME2, NTT_PRIME1, NTT_PRIME2);
		}
	}
	LongInt::fix_carry_p(result, truncate);
}

LongInt LongInt::operator*(const LongInt& num)const&
{
	if (is_inf() || num.is_inf())return LongInt((bool)(sign ^ num.sign), true);
	LongInt tmp;
	multiply(data, num.data, tmp.data);
	while (tmp.data.size() > 1 && !tmp.data.back())
	{
		tmp.data.pop_back();
	}
	tmp.sign = sign ^ num.sign;
	return tmp;
}

LongInt& LongInt::operator*=(const LongInt& num)&
{
	if (is_inf())
	{
		sign ^= num.sign;
		return *this;
	}
	else if (num.is_inf())
	{
		sign ^= num.sign;
		data = num.data;
		return *this;
	}
	std::vector<long long> tmp;
	multiply(data, num.data, tmp);
	while (tmp.size() > 1 && !tmp.back())
	{
		tmp.pop_back();
	}
	data = std::move(tmp);
	sign ^= num.sign;
	return *this;
}

LongInt LongInt::operator/(const LongInt& num)const&
{
	LongInt tmp = *this;
	tmp /= num;
	return tmp;
}

LongInt LongInt::operator%(const LongInt& num)const&
{
	return *this - (*this / num) * num;
}

void divide(std::vector<long long>& num1, long long num2)
{
	long long prev = num1.back() % num2;
	num1.back() /= num2;
	for (std::size_t i = num1.size() - 1; i > 0;)
	{
		--i;
		num1[i] += prev * LINT_BASE;
		prev = num1[i] % num2;
		num1[i] /= num2;
	}
}

LongInt& LongInt::operator/=(const LongInt& num)&
{
	{
		bool _inf1 = is_inf(), _inf2 = num.is_inf();
		if (_inf1 && _inf2)
		{
			*this = LongInt();//未定義
			return *this;
		}
		else if (_inf1)
		{
			sign ^= num.sign;//符号が変転
			return *this;
		}
		else if (_inf2)
		{
			*this = LongInt();//0を返す
			return *this;
		}
	}
	if (data.size() < num.data.size())
	{
		*this = zero();
		return *this;
	}
	if (data.size() == 1)
	{
		data[0] /= num.data[0];
		if (data[0] != 0)sign ^= num.sign;
		else sign = false;
		return *this;
	}
	if (num.data.size() == 1)
	{
		divide(data, num.data[0]);
		del_zero();
		sign ^= num.sign;
		return *this;
	}
	std::vector<long long> mul = { num.data[num.data.size() - 2],num.data.back() };//numを現在の有効桁数で切り捨てたもの
	std::vector<long long> inv = { 0,1 }, tmp1, tmp2;
	std::size_t sagnif = 2;//有効桁数
	std::size_t next = 4;//次のループの有効桁数
	std::size_t i;
	std::size_t size = data.size() + 2 - num.data.size();//逆数の必要精度(有効桁数)
	long long carry = 0;
	while (true)
	{
		for (auto& element : tmp1)element = 0;
		multiply_nr(inv, mul, tmp1);
		carry = 0;
		//2から引く
		--next;
		for (i = 0; i < next; ++i)
		{
			tmp1[i] = LINT_BASE + carry - tmp1[i];
			carry = tmp1[i] / LINT_BASE - 1;
			tmp1[i] %= LINT_BASE;
		}
		tmp1[next] = 2 + carry - tmp1[next];
		++next;
		for (auto& element : tmp2)element = 0;
		multiply_nr(tmp1, inv, tmp2, 3);
		if (inv[1] == tmp2[1] && inv[0] == tmp2[0])break;
		inv[1] = tmp2[1]; inv[0] = tmp2[0];
	}
	inv.push_back(0);
	//inv.size() == sagnif + 1
	//Bを除数とする
	//A = A * (2 - A * B)を繰り返すとAが逆数に近づく
	while (true)
	{
		for (auto& element : tmp1)element = 0;
		multiply(mul, inv, tmp1);//tmp1 = A * B
		//tmp1.size() == next + 1
		tmp1.pop_back();//最上位桁は0
		//tmp1 = 2 - tmp1
		carry = 0;
		--next;
		for (i = 0; i < next; ++i)
		{
			tmp1[i] = LINT_BASE + carry - tmp1[i];
			carry = tmp1[i] / LINT_BASE - 1;
			tmp1[i] %= LINT_BASE;
		}
		tmp1[next] = 2 + carry - tmp1[next];
		++next;
		//tmp1.size() == next
		//下位sagnif - 1桁切り捨て
		multiply(tmp1, inv, tmp2, sagnif - 1);//tmp2 = A * tmp1 == A * (2 - A * B)
		//tmp2.size() == next + (sagnif + 1) - (sagnif - 1) == next + 2 
		tmp2.pop_back();
		//tmp2.size() == next + 1
		inv = std::move(tmp2);
		if (sagnif > size)
		{
			break;
		}
		mul.resize(next);
		if (next <= num.data.size())
		{
			for (i = 0; i < next; ++i)mul[i] = num.data[num.data.size() - next + i];
		}
		else
		{
			for (i = 0; i < num.data.size(); ++i)
			{
				mul[next - num.data.size() + i] = num.data[i];
			}
			for (i = next - num.data.size(); i > 0;)
			{
				--i;
				mul[i] = 0;
			}
		}
		tmp2 = std::vector<long long>();
		sagnif <<= 1;
		next <<= 1;
	}
	++size;
	if (inv[inv.size() - size - 1] >= LINT_BASE / 2)carry = 1;
	else carry = 0;
	for (i = 0; i < size; ++i)
	{
		inv[i] = inv[inv.size() - size + i];
	}
	inv.resize(size);
	for (auto& element : inv)
	{
		element += carry;
		if (element < LINT_BASE)break;
		else element %= LINT_BASE;
	}
	LongInt p;
	p.data = std::move(inv);
	std::size_t s = data.size() + 1;
	*this *= p;
	//サイズが小さいときは0
	if (data.size() <= s)return *this = zero();
	//桁合わせ
	for (i = s; i < data.size(); ++i)data[i - s] = data[i];
	data.resize(data.size() - s);
	sign ^= num.sign;
	return *this;
}

LongInt& LongInt::operator%=(const LongInt& num)&
{
	*this -= (*this / num) * num;
	return *this;
}

LongInt(std::pow)(const LongInt& base, unsigned long long expo)
{
	if (base.is_inf())return LongInt(base.get_sign() && (bool)(expo % 2), true);
	if (expo == 0)return LongInt(1);
	if (expo == 1)return LongInt(base);
	std::vector<LongInt> vec;
	constexpr unsigned long long mask = 0x8000000000000000;
	unsigned char bit_remaining = 64;
	while (!(expo & mask))
	{
		expo <<= 1;
		--bit_remaining;
	}
	LongInt tmp = base;
	expo <<= 1;
	--bit_remaining;
	while (bit_remaining != 0)
	{
		tmp *= tmp;
		if (expo & mask)tmp *= base;
		expo <<= 1;
		--bit_remaining;
	}
	return tmp;
}

void product(std::vector<LongInt>& integers)
{
	std::size_t to = integers.size();
	while (to != 1)
	{
		std::size_t mid = to / 2;
		for (std::size_t i = 0; i < mid; ++i)
		{
			integers[i] *= integers.back();
			integers.pop_back();
		}
		to = integers.size();
	}
}

LongInt factorial(unsigned int num)
{
	if (num < 2)return LongInt(1);
	std::list<unsigned int> prime_num = number::get_prime_number(num);
	std::vector<unsigned int> vexp(prime_num.size()), vprime_num(prime_num.begin(), prime_num.end());
	prime_num = std::list<unsigned int>();
	std::size_t i;
	for (i = 0; i < vprime_num.size(); ++i)
	{
		unsigned int tmp = 1;
		while (tmp <= num / vprime_num[i])
		{
			tmp *= vprime_num[i];
			vexp[i] += num / tmp;
		}
	}
	unsigned int vexp_max = vexp[0], prod_num = 1;
	std::vector<LongInt> prod_vec1, prod_vec2;
	while (vexp_max >= prod_num)
	{
		while (!vexp.empty() && vexp.back()==1)
		{
			prod_vec1.push_back(std::pow(LongInt(vprime_num.back()), prod_num));
			vprime_num.pop_back();
			vexp.pop_back();
		}
		for (i = vprime_num.size(); i > 0;)
		{
			--i;
			if (vexp[i] & (unsigned int)1)prod_vec1.push_back((std::pow)(LongInt(vprime_num[i]), prod_num));
			vexp[i] >>= 1;
		}
		product(prod_vec1);
		prod_vec2.push_back(std::move(prod_vec1[0]));
		prod_vec1 = std::vector<LongInt>();
		prod_num <<= 1;
	}
	product(prod_vec2);
	return prod_vec2[0];
}

std::ostream& operator<<(std::ostream& output, const LongInt& numout)
{
	std::string tmp;
	numout.tostring(tmp);
	output << tmp;
	return output;
}


std::istream& operator>>(std::istream& input, LongInt& numin)
{
	std::string tmp;
	input >> tmp;
	numin = tmp.c_str();
	return input;
}