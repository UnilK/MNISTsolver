#ifndef RANDOM_FUNCS_HPP_
#define RANDOM_FUNCS_HPP_

#include <algorithm>

namespace randoms{

	template <class T>  T float24(){
		// returns a floating point nmber in range ]-1, 1[
		T ret = (T)(std::rand()&((1<<24)-1))/(1<<24);
		if(rand()&1) ret = -ret;
		return ret;
	}

}

#endif
