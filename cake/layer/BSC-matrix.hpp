#ifndef BSC_MATRIX_LAYER_HPP_
#define BSC_MATRIX_LAYER_HPP_

#include <vector>
#include <algorithm>
#include <fstream>

#include "base.hpp"
#include "base-reversible.hpp"

using std::vector;
using std::ifstream;
using std::ofstream;

const int32_t BSC_MATRIX_LAYER_ID = 0x0040;

template<class T> class BSCMatrixLayer: public ReversibleLayer<T>{
	
	protected:

		T one;
		vector<vector<T> > mx, mxC;
		vector<T> bias, sens, biasC, sensC, slope, ucv;

	public:

		/*

		   First, a bias is added to the vector:
		   v[i] += bias[i];

		   Then it's multiplied with sensetivity coefficients:
		   v[i] *= sens[i];

		   Then it's run through some compression funtion f:
		   v[i] = f(v[i])
		   
		   Last, the data is shuffled with a matrix on to the next layer.
		   Each node i from this layer effects each node j in the
		   next layer with some coefficient mx[i][j].

		*/

		BSCMatrixLayer(){ this->id = BSC_MATRIX_LAYER_ID; }
		
		BSCMatrixLayer(int32_t n_, int32_t m_, T zero_, T one_) : ReversibleLayer<T>(n_, m_, zero_){
			this->id = BSC_MATRIX_LAYER_ID;

			this->one = one_;
			
			this->bias.resize(n_, zero_);
			this->biasC.resize(n_, zero_);
			this->sens.resize(n_, one_);
			this->sensC.resize(n_, one_);

			this->slope.resize(n_, one_);
			this->ucv.resize(n_, zero_);

			this->connect_next(m_);
		
			this->init_config();
		}
		
		BSCMatrixLayer(int32_t n_, T zero_ = (T)0, T one_ = (T)1) : ReversibleLayer<T>(n_, zero_){
			this->id = BSC_MATRIX_LAYER_ID;
			
			this->one = one_;

			this->bias.resize(n_, zero_);
			this->biasC.resize(n_, zero_);
			this->sens.resize(n_, one_);
			this->sensC.resize(n_, one_);
			
			this->slope.resize(n_, one_);
			this->ucv.resize(n_, zero_);

			this->m = 0;
		
			this->init_config();
		}
		
		BSCMatrixLayer(ifstream &get_in){
			this->variables_in(get_in);
		}

		~BSCMatrixLayer(){}

		void init_config(){	
			this->configClar = {
				"matrix_change_speed:",
				"bias_change_speed:",
				"sensetivity_change_speed:"
			};
			this->config = {(T)0.01, (T)0.01, (T)0.001};
		}

		void connect_next(int32_t m_){
			this->m = m_;
			vector<T> tmp(m_, this->zero);
			this->mx.resize(this->n, tmp);
			this->mxC.resize(this->n, tmp);
		}
		
		void project_next(Layer<T> *next){
			
			next->set_vector_all(this->zero);
			
			for(int32_t i=0; i<this->n; i++){
				
				this->v[i] = (this->v[i]+this->bias[i])*this->sens[i];
				this->ucv[i] = this->v[i];

				for(int32_t j=0; j<this->m; j++){
					next->v[j] += this->mx[i][j]*this->v[i];
				}
			}
		}

		void variables_in(ifstream &get_in){
			
			if(!get_in.good()) return;

			get_in >> this->id >> this->n >> this->m >> this->zero >> this->one;
			
			this->config_in(get_in);
			
			this->connect_next(this->m);
			
			this->v.resize(this->n, this->zero);
			this->vC.resize(this->n, this->zero);
			this->bias.resize(this->n, this->zero);
			this->biasC.resize(this->n, this->zero);
			this->sens.resize(this->n, this->one);
			this->sensC.resize(this->n, this->one);

			this->slope.resize(this->n, this->one);
			this->ucv.resize(this->n, this->zero);
			
			for(int32_t i=0; i<this->n; i++) get_in >> this->bias[i];
			for(int32_t i=0; i<this->n; i++) get_in >> this->sens[i];	
			for(int32_t i=0; i<this->n; i++){
				for(int32_t j=0; j<this->m; j++) get_in >> this->mx[i][j];
			}

		}

		void variables_out(ofstream &get_out){
			
			get_out << this->id << ' ' << this->id << '\n';
			get_out << this->n << ' ' << this->m << ' ' << this->zero << ' ' << this->one << '\n';
			
			this->config_out(get_out, 0);
			
			for(int32_t i=0; i<this->n; i++) get_out << this->bias[i] << ' ';
			get_out << '\n';
			for(int32_t i=0; i<this->n; i++) get_out << this->sens[i] << ' ';
			get_out << '\n';
			
			for(int32_t i=0; i<this->n; i++){
				for(int32_t j=0; j<this->m; j++){
					get_out << this->mx[i][j] << ' ';
				} get_out << '\n';
			}
		}

		void set_variables(T val){
			
			for(int32_t i=0; i<this->n; i++){
				this->bias[i] = val;
				this->sens[i] = val;
				for(int32_t j=0; j<this->m; j++){
					this->mx[i][j] = val;
				}
			}
		}
		
		void random_variables(T (*random_func)(void)){
			
			for(int32_t i=0; i<this->n; i++){
				for(int32_t j=0; j<this->m; j++){
					this->mx[i][j] = random_func();
				}
			}
		}

		void downscale_changes(T down){
			for(int32_t i=0; i<this->n; i++){
				this->biasC[i] /= down;
				this->sensC[i] /= down;
				for(int32_t j=0; j<this->m; j++){
					this->mxC[i][j] /= down;
				}
			}
		}
		
		void zero_changes(){

			for(int32_t i=0; i<this->n; i++){

				this->vC[i] = this->zero;
				this->biasC[i] = this->zero;
				this->sensC[i] = this->zero;

				for(int32_t j=0; j<this->m; j++){
					this->mxC[i][j] = this->zero;
				}
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

			for(int32_t i=0; i<this->n; i++){
				this->vC[i] *= this->slope[i];
				this->biasC[i] += this->vC[i];
				this->sensC[i] += this->vC[i]*this->ucv[i];
				this->vC[i] *= this->sens[i];
			}
		}

		void adjust(){
			
			for(int32_t i=0; i<this->n; i++){
				this->bias[i] += this->config[1]*this->biasC[i];
				this->sens[i] += this->config[2]*this->sensC[i];
				for(int32_t j=0; j<this->m; j++){
					this->mx[i][j] += this->config[0]*this->mxC[i][j];
				}
			}
		}	
};

#endif
