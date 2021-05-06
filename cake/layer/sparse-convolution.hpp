#ifndef SPARSE_CONVOLUTION_LAYER_HPP_
#define SPARSE_CONVOLUTION_LAYER_HPP_

#include <vector>
#include <algorithm>
#include <fstream>
#include <math.h>

#include "base.hpp"
#include "base-reversible.hpp"
#include "../func/fft.hpp"

using std::vector;
using std::ifstream;
using std::ifstream;

const int32_t SPARSE_CONVOLUTION_LAYER_ID = 0x0031;

template<class T> class SparseConvolutionLayer: public ReversibleLayer<T>{

	protected:

		bool is_first_layer = 0;
		FFT<T> *fft;

		vector<T> cn, cnC;

	public:

		/*
		   almost the same as ConvolutionLayer, but
		   the convolution layer size is 2*n-1 and m values
		   are gathered sparsely from the output convolution.
		*/

		SparseConvolutionLayer(){ this->id = SPARSE_CONVOLUTION_LAYER_ID; }
		
		SparseConvolutionLayer(int32_t n_, int32_t m_, FFT<T> *fft_, T zero_, bool ifl_ = 0) : 
				ReversibleLayer<T>(n_, m_, zero_){
			this->id = SPARSE_CONVOLUTION_LAYER_ID;
			this->fft = fft_;
			this->is_first_layer = ifl_;
			
			this->init_config();
		}
		
		SparseConvolutionLayer(int32_t n_, FFT<T> *fft_, T zero_ = (T)0, bool ifl_ = 0) : 
				ReversibleLayer<T>(n_, zero_){
			this->id = SPARSE_CONVOLUTION_LAYER_ID;
			this->fft = fft_;
			this->is_first_layer = ifl_;
			
			this->init_config();
		}
		
		SparseConvolutionLayer(ifstream &get_in, FFT<T> *fft_){
			this->fft = fft_;
			this->variables_in(get_in);
		}

		~SparseConvolutionLayer(){}

		void init_config(){
			this->configClar = {"convolution_change_speed:"};
			this->config = {(T)0.01};
		}

		void connect_next(int32_t m_){
			this->m = m_;
			this->cn.resize(2*this->n-1, this->zero);
			this->cnC.resize(2*this->n-1, this->zero);
		}
		
		void project_next(Layer<T> *next){
			
			vector<T> conv = this->fft->convolution(this->v, this->cn, 2*this->n-1);
			
			float jump = (float)this->n/this->m, pos = 0;

			for(int32_t i=0; i<this->m; i++){
				next->v[i] = conv[this->n-1+(int32_t)std::floor(pos)];
				pos += jump;
			}
			
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

			   Here's the brute force implementation:

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

			if(this->cn_speed == (T)0.0) return;

			vector<T> changes(this->n, (T)0);
			
			float jump = (float)this->n/this->m, pos = 0;

			for(int32_t i=0; i<this->m; i++){
				this->v[this->n-1+(int32_t)std::floor(pos)] = feedback[i];
				pos += jump;
			}
			
			vector<T> ccnC = this->fft->convolution(this->v, changes, this->n+this->m-1, 1, 0);
			for(int32_t i=0; i<this->n+this->m-1; i++) this->cnC[i] += ccnC[i];
			
			if(!this->is_first_layer){
				vector<T> cvC = this->fft->convolution(this->cn, changes, this->n+this->m-1, 1, 0);
				for(int32_t i=0; i<this->n; i++) this->vC[i] += cvC[i+this->m-1];
			}
			
		}

		void adjust(){
			for(int32_t i=0; i<this->n+this->m-1; i++) this->cn[i] += this->cnC[i]*this->config[0];
		}	
};

#endif
