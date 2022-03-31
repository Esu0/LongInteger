#pragma once

#include<vector>
#include"modulo.h"

namespace discrete
{
#pragma warning(push)
#pragma warning(disable : 4244)
	using namespace number;
	class Number_Theoretic_Transform
	{
	private:
		using ull = unsigned long long;

		std::vector<unsigned int> Root;//2^(m + 1)æª = Root[m]
		std::vector<unsigned int> inv_Root;//“¯—l
		unsigned int Div = 998244353;
		unsigned int Primitive_root = 3;
		unsigned int Times = 1;

	public:
		Number_Theoretic_Transform() : Root(), inv_Root()
		{
			unsigned int root_n = Div - 1, inv_root_n;
			while (!(root_n % 2))
			{
				root_n >>= 1;
				Times <<= 1;
				Root.push_back(0);
				inv_Root.push_back(0);
			}
			root_n = mpow(Primitive_root, root_n, Div);
			inv_root_n = mpow(root_n, Times - 1, Div);
			for (unsigned int i = 0; i < Root.size(); ++i)
			{
				Root[i] = mpow(root_n, Times >> (i + 1), Div);
				inv_Root[i] = mpow(inv_root_n, Times >> (i + 1), Div);
			}
		}

		Number_Theoretic_Transform(unsigned int mod) : Div{mod}, Root(), inv_Root()
		{
			Primitive_root = get_prim_root(mod);
			unsigned int root_n = Div - 1, inv_root_n;
			while (!(root_n % 2))
			{
				root_n >>= 1;
				Times <<= 1;
				Root.push_back(0);
				inv_Root.push_back(0);
			}
			root_n = mpow(Primitive_root, root_n, Div);
			inv_root_n = mpow(root_n, Times - 1, Div);
			for (unsigned int i = 0; i < Root.size(); ++i)
			{
				Root[i] = mpow(root_n, Times >> (i + 1), Div);
				inv_Root[i] = mpow(inv_root_n, Times >> (i + 1), Div);
			}
		}

		Number_Theoretic_Transform(unsigned int mod, unsigned int primitive) : Div{ mod }, Primitive_root{ primitive }, Root(), inv_Root()
		{
			unsigned int root_n = Div - 1, inv_root_n;
			while (!(root_n % 2))
			{
				root_n >>= 1;
				Times <<= 1;
				Root.push_back(0);
				inv_Root.push_back(0);
			}
			root_n = mpow(Primitive_root, root_n, Div);
			inv_root_n = mpow(root_n, Times - 1, Div);
			for (unsigned int i = 0; i < Root.size(); ++i)
			{
				Root[i] = mpow(root_n, Times >> (i + 1), Div);
				inv_Root[i] = mpow(inv_root_n, Times >> (i + 1), Div);
			}
		}

		unsigned int getRoot(unsigned int index)const&
		{
			return Root[index];
		}

		unsigned int getiRoot(unsigned int index)const&
		{
			return inv_Root[index];
		}

		//vec‚Ì—v‘f”‚ª2^n‚Ì‚Æ‚«m = n - 1
		void operator()(std::vector<ull>& vec, int m, bool inv, ull* copy, ull begin = 0, ull e = 1)const&
		{
			static ull i, j, k;
			if (m < 0)return;
			operator()(vec, m - 1, inv, copy, begin, e << 1);
			operator()(vec, m - 1, inv, copy, begin + e, e << 1);
			ull l = 1;
			for (i = begin; i < vec.size(); i += e)
			{
				copy[i] = vec[i];
			}
			j = begin;
			for (i = begin; j < vec.size(); i += e)
			{
				k = j + e;
				vec[i] = (copy[j] + (l * copy[k])) % Div;
				if (inv)l *= inv_Root[m];
				else l *= Root[m];
				l %= Div;
				j += e << 1;
			}
			for (j = begin; i < vec.size(); i += e)
			{
				k = j + e;
				vec[i] = (copy[j] + (l * copy[k])) % Div;
				if (inv)l *= inv_Root[m];
				else l *= Root[m];
				l %= Div;
				j += e << 1;
			}
		}

		//“Y‚¦š‚ÍƒrƒbƒgƒŠƒo[ƒX‚³‚ê‚Ä•Ô‚³‚ê‚é
		//ü”g”ŠÔˆø‚«Œ^FMT
		template<typename T>
		void trans_f(std::vector<T>& vec, int m, ull from, ull to)const&
		{
			T tmp;
			if (m == -1)return;
			else if (m <= 0)
			{
				tmp = vec[from];
				vec[from] = (vec[from] + vec[to]) % (T)Div;
				vec[to] = (tmp + (T)Div - vec[to]) % (T)Div;
				return;
			}
			ull mid = (to - from + 1) / 2;
			T w = 1;
			for (ull i = from + mid; i <= to; ++i)
			{
				tmp = vec[i - mid];
				vec[i - mid] = (tmp + vec[i]) % (T)Div;
				vec[i] = ((tmp + (T)Div - vec[i]) * w) % (T)Div;
				w *= (T)Root[m];
				w %= (T)Div;
			}
			trans_f<T>(vec, m - 1, from, from + mid - 1);
			trans_f<T>(vec, m - 1, from + mid, to);
		}

		//4‚Â“¯‚Éü”g”ŠÔˆø‚«Œ^FMT
		template<typename T>
		void trans_f(std::vector<T>& vec1, std::vector<T>& vec2, std::vector<T>& vec3, std::vector<T>& vec4, int m, ull from, ull to)const&
		{
			T tmp;
			if (m == -1)return;
			else if (m == 0)
			{
				tmp = vec1[from];
				vec1[from] = (vec1[from] + vec1[to]) % (T)Div;
				vec1[to] = (tmp + (T)Div - vec1[to]) % (T)Div;
				tmp = vec2[from];
				vec2[from] = (vec2[from] + vec2[to]) % (T)Div;
				vec2[to] = (tmp + (T)Div - vec2[to]) % (T)Div;
				tmp = vec3[from];
				vec3[from] = (vec3[from] + vec3[to]) % (T)Div;
				vec3[to] = (tmp + (T)Div - vec3[to]) % (T)Div;
				tmp = vec4[from];
				vec4[from] = (vec4[from] + vec4[to]) % (T)Div;
				vec4[to] = (tmp + (T)Div - vec4[to]) % (T)Div;
				return;
			}
			ull mid = (to - from + 1) / 2;
			T w = 1;
			for (ull i = from + mid; i <= to; ++i)
			{
				tmp = vec1[i - mid];
				vec1[i - mid] = (tmp + vec1[i]) % (T)Div;
				vec1[i] = ((tmp + (T)Div - vec1[i]) * w) % (T)Div;
				tmp = vec2[i - mid];
				vec2[i - mid] = (tmp + vec2[i]) % (T)Div;
				vec2[i] = ((tmp + (T)Div - vec2[i]) * w) % (T)Div;
				tmp = vec3[i - mid];
				vec3[i - mid] = (tmp + vec3[i]) % (T)Div;
				vec3[i] = ((tmp + (T)Div - vec3[i]) * w) % (T)Div;
				tmp = vec4[i - mid];
				vec4[i - mid] = (tmp + vec4[i]) % (T)Div;
				vec4[i] = ((tmp + (T)Div - vec4[i]) * w) % (T)Div;
				w *= (T)Root[m];
				w %= (T)Div;
			}
			trans_f<T>(vec1, vec2, vec3, vec4, m - 1, from, from + mid - 1);
			trans_f<T>(vec1, vec2, vec3, vec4, m - 1, from + mid, to);
		}

		//ŠÔŠÔˆø‚«Œ^‹tFMT
		template<typename T>
		void itrans_t(std::vector<T>& vec, int m, ull from, ull to)const&
		{
			T tmp;
			if (m == -1)return;
			else if (m <= 0)
			{
				tmp = vec[from];
				vec[from] = (vec[from] + vec[to]) % (T)Div;
				vec[to] = (tmp + (T)Div - vec[to]) % (T)Div;
				return;
			}
			ull mid = (to - from + 1) / 2;
			itrans_t<T>(vec, m - 1, from, from + mid - 1);
			itrans_t<T>(vec, m - 1, from + mid, to);
			T w = 1;
			for (ull i = from + mid; i <= to; ++i)
			{
				tmp = vec[i - mid];
				vec[i] *= w;
				vec[i] %= (T)Div;
				vec[i - mid] = (tmp + vec[i]) % (T)Div;
				vec[i] = (tmp + (T)Div - vec[i]) % (T)Div;
				w *= (T)inv_Root[m];
				w %= (T)Div;
			}
		}

		template<typename T>
		void itrans_t(std::vector<T>& vec1, std::vector<T>& vec2, int m, ull from, ull to)const&
		{
			T tmp;
			if (m == -1)return;
			else if (m <= 0)
			{
				tmp = vec1[from];
				vec1[from] = (vec1[from] + vec1[to]) % (T)Div;
				vec1[to] = (tmp + (T)Div - vec1[to]) % (T)Div;
				tmp = vec2[from];
				vec2[from] = (vec2[from] + vec2[to]) % (T)Div;
				vec2[to] = (tmp + (T)Div - vec2[to]) % (T)Div;
				return;
			}
			ull mid = (to - from + 1) / 2;
			itrans_t<T>(vec1, vec2, m - 1, from, from + mid - 1);
			itrans_t<T>(vec1, vec2, m - 1, from + mid, to);
			T w = 1;
			for (ull i = from + mid; i <= to; ++i)
			{
				tmp = vec1[i - mid];
				vec1[i] *= w;
				vec1[i] %= (T)Div;
				vec1[i - mid] = (tmp + vec1[i]) % (T)Div;
				vec1[i] = (tmp + (T)Div - vec1[i]) % (T)Div;
				tmp = vec2[i - mid];
				vec2[i] *= w;
				vec2[i] %= (T)Div;
				vec2[i - mid] = (tmp + vec2[i]) % (T)Div;
				vec2[i] = (tmp + (T)Div - vec2[i]) % (T)Div;
				w *= (T)inv_Root[m];
				w %= (T)Div;
			}
		}
	};
#pragma warning(pop)
}