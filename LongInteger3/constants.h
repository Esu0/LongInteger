#pragma once
#include"LongInt.h"
#include<fstream>

struct frac_lint_result
{
	LongInt P, Q, B, T;
};

//a/b
struct frac_lint
{
	LongInt a, b;
};
void calc_pi();
frac_lint calc_root(unsigned int);
frac_lint calc_root2(unsigned int);
void bsa_pi(frac_lint_result&, size_t, size_t);