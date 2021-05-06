#ifndef MATRIX_LAYER_HPP_
#define MATRIX_LAYER_HPP_

#include <vector>
#include <algorithm>
#include <fstream>

#include "base.hpp"
#include "base-reversible.hpp"

using std::vector;
using std::ifstream;
using std::ofstream;

const int32_t MATRIX_LAYER_ID = 0x0010;

template<class T> class MatrixLayer: public ReversibleLayer<T>{
	
	protected:

		vector<vector<T> > mx, mxC;

	public:

		/*
		   
		   The data is shuffled with a matrix on to the next layer.
		   Each node i from this layer effects each node j in the
		   next layer with some coefficient mx[i][j].

		*/

		MatrixLayer(){ this->id = MATRIX_LAYER_ID; }
		
		MatrixLayer(int32_t n_, int32_t m_, T zero_) : ReversibleLayer<T>(n_, m_, zero_){
			this->connect_next();
			this->id = MATRIX_LAYER_ID;
			this->init_config();
		}
		
		MatrixLayer(int32_t n_, T zero_ = (T)0) : ReversibleLayer<T>(n_, zero_){
			this->id = MATRIX_LAYER_ID;			
			this->init_config();
		}
		
		MatrixLayer(ifstream &get_in){
			this->variables_in(get_in);
		}

		virtual ~MatrixLayer(){}

		void init_config(){
			this->configClar = {"matrix_change_speed:"};
			this->config = {(T)0.01};
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

		void set_variables(T val){
			for(int32_t i=0; i<this->n; i++){
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
				for(int32_t j=0; j<this->m; j++){
					this->mxC[i][j] /= down;
				}
			}
		}
		
		void zero_changes(){

			for(T &i : this->vC) i = this->zero;

			for(int32_t i=0; i<this->n; i++){
				for(int32_t j=0; j<this->m; j++){
					this->mxC[i][j] = this->zero;
				}
			}
		}

		void evaluate(vector<T> feedback){

			/*
			   feedback shows the desired changes to variables in the
			   next layer.

			   how to read this:
			   The change to the edge between nodes i and j should
			   linearly correlate to the value of i and the desired change of j.

			   The desired change to node i should correlate to the desired
			   change of node j and the value of the edge between i and j.

			   This intuitively makes sense to me, so it's good enough.
			*/

			this->set_vector_changes(this->zero);

			for(int32_t i=0; i<this->n; i++){
				for(int32_t j=0; j<this->m; j++){
					this->mxC[i][j] += this->v[i]*feedback[j];
					this->vC[i] += this->mx[i][j]*feedback[j];
				}
			}
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
