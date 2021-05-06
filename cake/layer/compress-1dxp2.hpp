#ifndef COMPRESS_1DXP2_LAYER_HPP_
#define COMPRESS_1DXP2_LAYER_HPP_

#include <fstream>
#include <vector>
#include <algorithm>

#include "base.hpp"
#include "base-reversible.hpp"
#include "../func/compressions.hpp"

using std::vector;
using std::ifstream;
using std::ofstream;

const int32_t C_1DXP2_LAYER_ID = 0x0022;

template<class T> class C1dxp2Layer: public ReversibleLayer<T>{

	protected:

		vector<T> slope;

	public:

		/*
		   The values get compressed to the range [0, 1] with the
		   following function:
		   
			v[i] *= config[0];
			if(v[i] > 0){
				v[i] += 1;
				v[i] = 1 - 0.5/(v[i]*v[i]);
			} else {
				v[i] -= 1;
				v[i] = -0.5/(v[i]*v[i]);
			}

		*/
		
		C1dxp2Layer(){ this->id = C_1DXP2_LAYER_ID; }
		
		C1dxp2Layer(int32_t n_, int32_t m_, T zero_) : ReversibleLayer<T>(n_, m_, zero_){
			this->id = C_1DXP2_LAYER_ID;
			this->slope.resize(n_);
			this->init_config();
		}
		
		C1dxp2Layer(int32_t n_, T zero_ = (T)0) : ReversibleLayer<T>(n_, zero_){
			this->id = C_1DXP2_LAYER_ID;
			this->slope.resize(n_);
			this->init_config();
		}
		
		C1dxp2Layer(ifstream &get_in){
			this->variables_in(get_in);
		}
		
		~C1dxp2Layer(){}

		void init_config(){
			this->configClar = {"x-axis_compression:"};
			this->config = {(T)1};
		}

		void project_next(Layer<T> *next){

			compress::div_xp2<T>(this->v, this->slope, {this->config[0]});

			for(int32_t i=0; i<std::min(this->n, this->m); i++){
				next->v[i] = this->v[i];
			}
		}
		
		virtual void variables_in(ifstream &get_in){

			if(!get_in.good()) return;

			get_in >> this->id >> this->n >> this->m >> this->zero;
			
			this->config_in(get_in);
			
			this->v.resize(this->n, this->zero);
			this->connect_next(this->m);
			this->vC.resize(this->n, this->zero);
			this->slope.resize(this->n, this->zero);
		}

		virtual void variables_out(ofstream &get_out){
			get_out << this->id << ' ' << this->id << '\n';
			get_out << this->n << ' ' << this->m << ' ' << this->zero << '\n';
			
			this->config_out(get_out, 0);
		}

		void evaluate(vector<T> feedback){
			
			for(int32_t i=0; i<std::min(this->n, this->m); i++){
				this->vC[i] = this->slope[i]*feedback[i];
			}
		}
};

#endif
