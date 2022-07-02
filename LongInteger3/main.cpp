#include"LongInt.h"
#include<iostream>
#include<chrono>
#include"constants.h"
#include<fstream>

#define MODE_DEBUG

#ifndef MODE_DEBUG
#define MODE_TIME
#endif

using namespace std;
int main()
{
#ifdef MODE_DEBUG
	//LongInt a = -2000300000, b = -2000200;
	
	//a.random(10);
	//b.random(10);
	//a.absolute(); b.absolute();
	//cout << a << " * " << b << " = " << a * b << endl;
	//calc_pi();
	constexpr unsigned int N = 20;
	frac_lint res = calc_root2(N);
	LongInt root2 = res.a.shift((long long)(1.5 * (1 << N))) / res.b;
	std::ofstream ofs;
	ofs.open("root2.txt", std::ios::out);
	ofs << root2 << endl;
	cout << "digit size:" << root2.size() << endl;
	cout << "result.a size:" << res.a.size() << endl;
	cout << "result.b size:" << res.b.size() << endl;
	//ofs << root2 * root2 << endl;
#endif
#ifdef MODE_TIME
	chrono::system_clock::time_point start, end;
	LongInt a = "-12345678";
	LongInt b = "13572468";
	a.random(1000000);
	b.random(1000000);
	//cout << a << " * " << b << " = " << a * b << endl;
	start = chrono::system_clock::now();
	//factorial(1000000);
	a *= b;
	end = chrono::system_clock::now();
	double exetime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1000.0;
	cout << "execution time : " << exetime << " ms\n";
	//cout << a << endl;
#endif
}