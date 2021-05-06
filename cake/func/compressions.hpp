#ifndef COMPRESSIONS_HPP_
#define COMPRESSIONS_HPP_

#include <algorithm>
#include <math.h>
#include <vector>

using std::vector;

namespace compress{
	
	template<class T> void div_x(vector<T> &v, vector<T> &dv, vector<T> c){

		/*
		   values get compressed to range [0, 1].
		   Nice & continuous derivative. c[0] is
		   the x-axis squishification factor. The derivative
		   of the compression function is stored in dv.
		*/

		int32_t n = v.size();

		for(int32_t i=0; i<n; i++){
			v[i] *= c[0];
			if(v[i] > (T)0){
				v[i] += 1;
				dv[i] = ((T)0.5*c[0])/(v[i]*v[i]);
				v[i] = (T)1 - (T)0.5/v[i];
			} else {
				v[i] -= 1;
				dv[i] = ((T)0.5*c[0])/(v[i]*v[i]);
				v[i] = -(T)0.5/v[i];
			}
		}
	}
	
	template<class T> void div_xp2(vector<T> &v, vector<T> &dv, vector<T> c){

		// same idea as div_x, but the compression function is different

		int32_t n = v.size();

		for(int32_t i=0; i<n; i++){
			v[i] *= c[0];
			if(v[i] > (T)0){
				v[i] += 1;
				dv[i] = c[0]/(v[i]*v[i]*v[i]);
				v[i] = (T)1 - (T)0.5/(v[i]*v[i]);
			} else {
				v[i] -= 1;
				dv[i] = -c[0]/(v[i]*v[i]*v[i]);
				v[i] = (T)0.5/(v[i]*v[i]);
			}
		}
	}
	
	template<class T> void logistic(vector<T> &v, vector<T> &dv, vector<T> c){

		// Values get compressed to range [0, 1] with a logistic curve.

		int32_t n = v.size();

		for(int32_t i=0; i<n; i++){
			v[i] = std::exp(-c[0]*v[i]);
			if(v[i] < 1e20) dv[i] = (c[0]*v[i])/((v[i]+1)*(v[i]+1));
			else dv[i] = 0;
			v[i] = (T)1/(v[i]+1);
		}
	}

	

}

#endif
