#include"constants.h"
#include<chrono>

void calc_pi()
{
	size_t n = 4096;
	frac_lint_result result;
	frac_lint root10005 = calc_root(11);

	std::chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	bsa_pi(result, 0, n);
	LongInt pi = (result.Q * root10005.a * 426880).shift(n * 14) / (result.T * root10005.b);

	end = std::chrono::system_clock::now();
	double exetime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1000;

	std::cout << exetime << "ms" << std::endl;

	std::ofstream ofs;
	ofs.open("result.txt", std::ios::out);
	ofs << pi << std::endl;
	ofs.close();
}

frac_lint calc_root(unsigned int N)
{
	static constexpr unsigned int X = 10005;
	unsigned int n = 100;
	LongInt a = n * n + X, b = 2 * n;
	for (size_t i = 0; i < N; ++i)
	{
		LongInt tmp = std::move(a);
		a = tmp * tmp + b * X * b;
		b *= tmp * 2;
	}
	frac_lint res;
	res.a = std::move(a); res.b = std::move(b);
	return res;
}

frac_lint calc_root2(unsigned int N)
{
	frac_lint res;
	res.a = 99, res.b = 70;
	for (size_t i = 0; i < N; ++i)
	{
		LongInt tmp = std::move(res.a);
		res.a = tmp * tmp + res.b * 2 * res.b;
		res.b *= tmp * 2;
	}
	return res;
}
//[l,r)‚Ì•”•ª”—ñŒvŽZ
void bsa_pi(frac_lint_result& res, size_t l, size_t r)
{
	const LongInt A = LongInt(13591409);
	const LongInt B = LongInt(545140134);
	//const LongInt C = LongInt(262537412640768000);
	switch (r - l)
	{
	case 0:
		res.P = 1;
		res.Q = 1;
		res.B = 1;
		res.T = LongInt();
		break;

	case 1:
		if (l == 0)
		{
			res.P = 1;
			res.Q = 1;
			res.B = 1;
			res.T = A;
		}
		else
		{
			res.P = LongInt(-(long long)(6 * l - 5) * (2 * l - 1) * (6 * l - 1));
			res.Q = (LongInt(10939058860032000) * l) * (LongInt(l) * l);
			res.B = 1;
			res.T = (A + B * l) * res.P;
		}
		break;
		/*
	case 2:

		break;
	case 3:
		break;
		*/
	default:
		frac_lint_result _Left, _Right;
		size_t mid = (l + r) / 2;
		bsa_pi(_Left, l, mid);
		bsa_pi(_Right, mid, r);

		res.P = _Left.P * _Right.P;
		res.Q = _Left.Q * _Right.Q;
		res.B = _Left.B * _Right.B;
		res.T = _Left.B * _Left.P * _Right.T + _Right.B * _Right.Q * _Left.T;
		break;
	}
}