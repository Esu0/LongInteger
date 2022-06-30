#include"constants.h"

void calc_pi()
{
	size_t n = 32;
	frac_lint_result result;
	bsa_pi(result, 0, n);
	LongInt pi = result.Q * ("80019997500624804755833750301086") * "1000000000000000000" * 640320 / (result.T * 12);
	std::ofstream ofs;
	ofs.open("result.txt", std::ios::out);
	ofs << pi << std::endl;
	ofs.close();
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