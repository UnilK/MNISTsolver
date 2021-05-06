#ifndef FFT_HPP_
#define FFT_HPP_

#include <vector>
#include <algorithm>
#include <complex>

using std::vector;
using std::complex;

const long double PI_FFT = 3.14159265358979323;

template<class T> class FFT{

	/*

	   A class for calculating Fourier transforms
	   reasonably efficiently. Explaining FFT in
	   a comment is quite a tall order, so here
	   are some sources to learn how it works, if
	   you aren't familiar with it:

	   Source that I used to learn it:
       https://codeforces.com/blog/entry/43499

	   3Blue1Brown's works on related topics:
	   https://www.youtube.com/watch?v=v0YEaeIClKY
       https://www.youtube.com/watch?v=spUNpyF58BY
       https://www.youtube.com/watch?v=r6sGWTCMz2k
	   etc.

	   Reducible's video on FFT, the take is similar to the blog post's:
       https://www.youtube.com/watch?v=h7apO7q16V0&t=902s

	   The implementation I use is referred
	   to as the "bit reverse" technique, I believe.
	   Comments on the implementation assume that the reader
	   understands the basic recursive implementation of
	   fft well.

	*/

	protected:

		int32_t B = 0;
		vector<vector<complex<T> > > w;
		vector<int32_t > invbit;		

	public:

		FFT(){
			B = 0;
			resize_precalc_tables();
		}

		FFT(int32_t B_){
			B = B_;
			resize_precalc_tables();
		}

		void resize_precalc_tables(){
			
			while((int32_t)w.size() <= B){
				
				/*

				   precalculate w[b][x] = e^(i*pi*x/(2^b))

				   and the invbit array:
				   examples in base 2 of the invbit array:

				   invbit[1000111] = 1110001
				   invbit[0] = 0
				   invbit[10110111] = 11101101

				   it reverses the bit order.

				*/

				int32_t z = 1<<w.size();
				vector<complex<T> > prec(z);
				for(int32_t i=0; i<z; i++) prec[i] = std::polar((T)1, (T)(PI_FFT*i/z));
				w.push_back(prec);

				invbit.resize(z, 0);
				for(int i=0; i<z/2; i++){
					invbit[i] <<= 1;
					invbit[i+z/2] = invbit[i]+1;
				}
			}

		}

		void fft(vector<complex<T> > &v){
			
			/*
			   let's examine how the recursive FFT changes
			   the ordering of the input array:
			   
			   [000, 001, 010, 011, 100, 101, 110, 111]
			   [000, 010, 100, 110], [001, 011, 101, 111]
			   [000, 100], [010, 110], [001, 101], [011, 111]
			   [000], [100], [010], [110], [001], [101], [011], [111]

			   [000, 100, 010, 110, 001, 101, 011, 111]
			  
			   The index of element x in the resulting array
			   is it's index in the inupt array, but with the binary
			   representation flipped. This makes sense, since
			   on each layer we sort the array by the least
			   significant bit.

			   Just reorder the array, then run the basic
			   FFT algorithm on a loop.

			   O(1) extra memory used with O(N*log(N)) preclac,
			   O(N*log(N)) time with a good constant compared to
			   the recursive implementation.
			*/

			int32_t n = v.size(), b = 0;
			while((1<<b) < n) b++;

			// the invbit array for b-1 is the same as for b,
			// but each element is shifted back by 1 bit.
			int32_t shift = B-b;

			for(int32_t i=0; i<n; i++){
				if(i < invbit[i]>>shift) std::swap(v[i], v[invbit[i]>>shift]);
			}

			for(int32_t r=0; r<b; r++){
				int32_t rd = 1<<r;
				for(int32_t i=0; i<n; i+=2*rd){
					for(int32_t j=i; j<i+rd; j++){
						complex<T> tmp = w[r][j-i]*v[j+rd];
						v[j+rd] = v[j]-tmp;
						v[j] = v[j]+tmp;
					}
				}
			}
		}

		vector<T> convolution(
				vector<T> &x, vector<T> &y,
				int32_t n=0, bool inv1=0, bool inv2=0){
			
			int32_t b = 0, zx = x.size(), zy = y.size();
			if(!n) n = zx+zy-1;
			
			n = std::max(n, zx);
			n = std::max(n, zy);
			
			while(1<<b < n) b++;

			B = std::max(B, b);
			resize_precalc_tables();

			vector<complex<T> > cx(1<<b, {0, 0}), cy(1<<b, {0, 0});

			if(inv1) for(int32_t i=zx-1; i>=0; i--) cx[zx-1-i] = {x[i], 0};
			else for(int32_t i=0; i<zx; i++) cx[i] = {x[i], 0};
			
			if(inv2) for(int32_t i=zy-1; i>=0; i--) cy[zy-1-i] = {y[i], 0};
			else for(int32_t i=0; i<zy; i++) cy[i] = {y[i], 0};

			fft(cx);
			fft(cy);

			// clalculating convolution
			for(int32_t i=0; i<(1<<b); i++) cx[i] *= cy[i];

			// inverse
			fft(cx);
			std::reverse(cx.begin()+1, cx.end());
			vector<T> xy(n);
			for(int32_t i=0; i<n; i++) xy[i] = cx[i].real()/(T)(1<<b);

			return xy;
		}
};

#endif
