#pragma once

#include<vector>
#include<list>
#include<random>

namespace number
{
	unsigned long long mpow(unsigned long long base, unsigned long long exp, unsigned long long div);
	void gcd_extend(long long a, long long b, long long& x, long long& y);
	long long minv(long long a, long long mod);
	long long garner(long long r1, long long r2, long long m1, long long m2);
	unsigned int get_prim_root(unsigned int mod);
	std::list<unsigned int> get_prime_number(unsigned int num);
}