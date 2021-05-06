#ifndef CONVOLUTION_LAYER_HPP_
#define CONVOLUTION_LAYER_HPP_

#include <vector>
#include <algorithm>
#include <fstream>

#include "base.hpp"
#include "base-reversible.hpp"
#include "../func/fft.hpp"

using std::vector;
using std::ifstream;
using std::ifstream;

const int32_t CONVOLUTION_LAYER_ID = 0x0030;

template<class T> class ConvolutionLayer: public ReversibleLayer<T>{

	protected:

		FFT<T> *fft;

		vector<T> cn, cnC;
		
		// since the convolution operation is rather heavy,
		// the evaluation operations are cutting off if they
		// are not needed.
		bool is_first_layer = 0;

	public:

		/*
		   the data is shuffled with a convolution on to the
		   next layer. Here's a couple of ways to understand
		   what's going on:

		   Interpret the input array and convolution array
		   as polynomial coefficients. Multiply these
		   polynomials together and pass m of the most
		   significant coefficients of the result on to the next layer.

		   Think of the input and convolution arrays as
		   images. The result is an image, where each pixel tells:
		   how well does the convolution image centered (not an exact term)
		   at this point match the input image?

		   The brute force solution to the convolution is 
		   for(int i=0; i<n; i++){
				for(int j=0; i+j<n+m-1; j++) result[i+j] += cn[j]*v[i]
		   }
		   for(int i=n-1; i<n+m-1; i++) next->v[i-n-1] = result[i];

		   FFT is used to calculate the convolutions in O(n*log(n))

		   the inputs for fft->convolution(va, vb, n, fa, fb) are:
		   va -> input vector 1
		   vb -> input vector 2
		   n -> size of the output vector, n >= va.size()+vb.size()-1
		   fa -> flip vector 1 before convolution
		   fb -> flip vector 2 before convolution
		*/

		ConvolutionLayer(){ this->id = CONVOLUTION_LAYER_ID; }
		
		ConvolutionLayer(int32_t n_, int32_t m_, FFT<T> *fft_, T zero_, bool ifl_ = 0) : 
				ReversibleLayer<T>(n_, m_, zero_){
			this->id = CONVOLUTION_LAYER_ID;
			this->fft = fft_;
			this->is_first_layer = ifl_;
			
			this->init_config();
		}
		
		ConvolutionLayer(int32_t n_, FFT<T> *fft_, T zero_ = (T)0, bool ifl_ = 0) : 
				ReversibleLayer<T>(n_, zero_){
			this->id = CONVOLUTION_LAYER_ID;
			this->fft = fft_;
			this->is_first_layer = ifl_;
			
			this->init_config();
		}
		
		ConvolutionLayer(ifstream &get_in, FFT<T> *fft_){
			this->fft = fft_;
			this->variables_in(get_in);
		}

		~ConvolutionLayer(){}

		void init_config(){
			this->configClar = {"convolution_change_speed:"};
			this->config = {(T)0.01};
		}

		void connect_next(int32_t m_){
			this->m = m_;
			this->cn.resize(this->n+this->m-1, this->zero);
			this->cnC.resize(this->n+this->m-1, this->zero);
		}
		
		void project_next(Layer<T> *next){
			
			vector<T> conv = this->fft->convolution(this->v, this->cn, this->n+this->m-1);
			for(int32_t i=0; i<this->m; i++) next->v[i] = conv[i+this->n-1];
			
		}
		
		void variables_in(ifstream &get_in){
			
			get_in >> this->id >> this->n >> this->m >> this->zero;

			if(!get_in.good()) return;

			this->v.resize(this->n, this->zero);
			this->vC.resize(this->n, this->zero);
			this->connect_next(this->m);
			
			this->config_in(get_in);
			
			for(int32_t i=0; i<this->n+this->m-1; i++) get_in >> this->cn[i];
		}

		void variables_out(ofstream &get_out){

			get_out << this->id << ' ' << this->id << '\n';
			get_out << this->n << ' ' << this->m << ' ' << this->zero << '\n';
			
			this->config_out(get_out, 0);
			
			for(int32_t i=0; i<this->n+this->m-1; i++) get_out << this->cn[i] << ' ';
			get_out << '\n';
		}

		void set_variables(T val){
			for(int32_t i=0; i<this->n+this->m-1; i++) this->cn[i] = val;
		}
		
		void random_variables(T (*random_func)(void)){
			for(int32_t i=0; i<this->n+this->m-1; i++) this->cn[i] = random_func();
		}

		void downscale_changes(T down){
			for(int32_t i=0; i<this->n+this->m-1; i++) this->cnC[i] /= down;
		}
		
		void zero_changes(){
			for(T &i : this->vC) i = this->zero;
			for(int32_t i=0; i<this->n+this->m-1; i++) this->cnC[i] = this->zero;
		}

		void evaluate(vector<T> &feedback){
			
			/*
			   The logic here is the exact same as in the matrix layer.

			   Here's the brute force implementation, I think:

			   for(int i=0; i<n; i++){
					for(int j=0; j<m; j++){
						cnC[n-1-i+j] += v[i]*feedback[j];
					}
			   }
			   
			   for(int i=0; i<n; i++){
			   		for(int j=0; j<m; j++){
						vC[i] += nc[n-i-1+j]*feedback[j]
					}
			   }

			   See the FFT class for info on spesifics.
			*/

			// optimization for cases when convolution doesn't change
			if(this->cn_speed == (T)0.0) return;

			vector<T> ccnC = this->fft->convolution(this->v, feedback, this->n+this->m-1, 1, 0);
			for(int32_t i=0; i<this->n+this->m-1; i++) this->cnC[i] += ccnC[i];
			
			if(!this->is_first_layer){
				vector<T> cvC = this->fft->convolution(this->cn, feedback, this->n+this->m-1, 1, 0);
				for(int32_t i=0; i<this->n; i++) this->vC[i] += cvC[i+this->m-1];
			}
		}

		void adjust(){
			for(int32_t i=0; i<this->n+this->m-1; i++) this->cn[i] += this->cnC[i]*this->config[0];
		}	
};

#endif
