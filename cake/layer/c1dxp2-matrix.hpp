#ifndef COMPRESS_1DXP2_MATRIX_LAYER_HPP_
#define COMPRESS_1DXP2_MATRIX_LAYER_HPP_

#include <fstream>
#include <vector>
#include <algorithm>

#include "base.hpp"
#include "base-reversible.hpp"
#include "matrix.hpp"
#include "../func/compressions.hpp"

using std::vector;
using std::ifstream;
using std::ofstream;

const int32_t C_1DXP2_MATRIX_LAYER_ID = 0x0012;

template<class T> class C1dxp2MatrixLayer: public MatrixLayer<T>{

	protected:

		vector<T> slope;

	public:

		/*
		   The values get compressed to the range [0, 1] with the
		   following function:
		   
			v[i] *= config[0];
			if(v[i] > 0){
				v[i] += 1;
				v[i] = 1 - 0.5/v[i];
			} else {
				v[i] -= 1;
				v[i] = -0.5/v[i];
			}

			then they are shuffled with a matrix on to the next layer.

		*/
		
		C1dxp2MatrixLayer(){ this->id = C_1DXP2_MATRIX_LAYER_ID; }
		
		C1dxp2MatrixLayer(int32_t n_, int32_t m_, T zero_) : MatrixLayer<T>(n_, m_, zero_){
			this->id = C_1DXP2_MATRIX_LAYER_ID;
			this->slope.resize(n_);
			this->init_config();
		}
		
		C1dxp2MatrixLayer(int32_t n_, T zero_ = (T)0) : MatrixLayer<T>(n_, zero_){
			this->id = C_1DXP2_MATRIX_LAYER_ID;
			this->slope.resize(n_);
			this->init_config();
		}
		
		C1dxp2MatrixLayer(ifstream &get_in){
			this->variables_in(get_in);
		}
		
		~C1dxp2MatrixLayer(){}

		void init_config(){
			this->configClar = {"matrix_change_speed:", "x-axis_compression:"};
			this->config = {(T)0.01, (T)1};
		}

		void project_next(Layer<T> *next){

			compress::div_xp2<T>(this->v, this->slope, {this->config[1]});
			
			next->set_vector_all(this->zero);

			for(int32_t i=0; i<this->n; i++){
				for(int32_t j=0; j<this->m; j++){
					next->v[j] += this->mx[i][j]*this->v[i];
				}
			}

		}
		
		void variables_in(ifstream &get_in){
			
			if(!get_in.good()) return;

			get_in >> this->id >> this->n >> this->m >> this->zero;
			
			this->config_in(get_in);
			
			this->v.resize(this->n, this->zero);
			this->connect_next(this->m);
			this->vC.resize(this->n, this->zero);
			this->slope.resize(this->n, this->zero);
			
			for(int32_t i=0; i<this->n; i++){
				for(int32_t j=0; j<this->m; j++) get_in >> this->mx[i][j];
			}

		}

		void variables_out(ofstream &get_out){
			
			get_out << this->id << ' ' << this->id << '\n';
			get_out << this->n << ' ' << this->m << ' ' << this->zero << '\n';
			
			this->config_out(get_out, 0);
			
			for(int32_t i=0; i<this->n; i++){
				for(int32_t j=0; j<this->m; j++){
					get_out << this->mx[i][j] << ' ';
				} get_out << '\n';
			}
		}
		
		void evaluate(vector<T> feedback){

			this->set_vector_changes(this->zero);

			for(int32_t i=0; i<this->n; i++){
				for(int32_t j=0; j<this->m; j++){
					this->mxC[i][j] += this->v[i]*feedback[j];
					this->vC[i] += this->mx[i][j]*feedback[j];
				}
			}

			for(int32_t i=0; i<this->n; i++) this->vC[i] *= this->slope[i];

		}

		void adjust(){
			for(int32_t i=0; i<this->n; i++){
				for(int32_t j=0; j<this->m; j++){
					this->mx[i][j] += this->config[0]*this->mxC[i][j];
				}
			}
		}	

};

#endif
