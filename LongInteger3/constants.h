#pragma once
#include"LongInt.h"
#include<fstream>

struct frac_lint_result
{
	LongInt P, Q, B, T;
};

void calc_pi();
void bsa_pi(frac_lint_result&, size_t, size_t);