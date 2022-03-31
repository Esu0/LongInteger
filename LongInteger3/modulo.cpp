#include "modulo.h"

unsigned long long number::mpow(unsigned long long base, unsigned long long exp, unsigned long long div)
{
	unsigned long long tmp;
	if (exp <= 0)return 1;
	if (exp == 1)return base;
	tmp = number::mpow(base, exp / 2, div);
	if (exp % 2)
		return (((tmp * tmp) % div) * base) % div;
	else 
		return (tmp * tmp) % div;
}

void number::gcd_extend(long long a, long long b, long long& x, long long& y)
{
	if (b == 1)
	{
		x = 1;
		y = 1 - a;
		return;
	}
	long long x2, y2;
	number::gcd_extend(b, a % b, x2, y2);
	x = y2;
	y = x2 - (a / b) * y2;
}

//aとmodは互いに素
long long number::minv(long long a, long long mod)
{
	long long b = mod, u = 1, v = 0;
	while (b)
	{
		long long tmp;

		tmp = u;
		u = v;
		v = tmp - (a / b) * v;

		tmp = a;
		a = b;
		b = tmp % b;
	}
	u %= mod;
	if (u < 0)u += mod;
	return u;
}

long long number::garner(long long r1, long long r2, long long m1, long long m2)
{
	r1 %= m1;
	r2 = ((r2 - r1) * number::minv(m1, m2)) % m2;
	if (r2 < 0)r2 += m2;
	r1 += r2 * m1;
	return r1;
}

unsigned int number::get_prim_root(unsigned int mod)
{
	std::vector<unsigned int> prime_factor = std::vector<unsigned int>();
	--mod;
	unsigned int m = mod;
	if (m % 2 == 0)
	{
		prime_factor.push_back(2);
		m /= 2;
	}
	while (m % 2 == 0)m /= 2;
	for (unsigned int i = 3; i * i < m; i += 2)
	{
		if (m % i == 0)
		{
			prime_factor.push_back(i);
			m /= i;
		}
		while (m % i == 0)m /= i;
	}
	if (m != 1)prime_factor.push_back(m);

	static std::random_device seed;
	static std::mt19937_64 mt(seed());
	std::uniform_int_distribution<> random(2, mod - 1);
	unsigned int candidate, prod;
	bool loopend;
	while (1)
	{
		candidate = random(mt);
		loopend = true;
		for (auto tmp : prime_factor)
		{
			prod = mod / tmp;
			if (mpow(candidate, prod, mod + 1) == 1)
			{
				loopend = false;
				break;
			}
		}
		if (loopend)break;
	}
	return candidate;
}

std::list<unsigned int> number::get_prime_number(unsigned int num)
{
	if (num < 2)return std::list<unsigned int>();
	std::list<unsigned int> num_list;
	unsigned int rt = (unsigned int)(std::sqrt)(num);
	for (unsigned long long i = 3; i <= (unsigned long long)num; i += 2)
	{
		num_list.push_back((unsigned int)i);
	}
	for (auto itr1 = num_list.begin(); itr1 != num_list.end() && *itr1 <= rt; ++itr1)
	{
		for (auto itr2 = std::next(itr1); itr2 != num_list.end();)
		{
			if (*itr2 % *itr1 == 0)
			{
				itr2 = num_list.erase(itr2);
			}
			else ++itr2;
		}
	}
	num_list.push_front(2);
	return num_list;
}