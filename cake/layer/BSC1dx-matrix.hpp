#ifndef BSC1DX_MATRIX_LAYER_HPP_
#define BSC1DX_MATRIX_LAYER_HPP_

#include <vector>
#include <algorithm>
#include <fstream>

#include "base.hpp"
#include "base-reversible.hpp"
#include "BSC-matrix.hpp"
#include "../func/compressions.hpp"

using std::vector;
using std::ifstream;

const int32_t BSC1DX_MATRIX_LAYER_ID = 0x0041;

template<class T> class BSC1dxMatrixLayer: public BSCMatrixLayer<T>{
	
	public:

		/*

		   First, a bias is added to the vector:
		   v[i] += bias[i];

		   Then it's multiplied with sensetivity coefficients:
		   v[i] *= sens[i];

		   Then it's run through the div_x compression function.
		   see func/compressions
		   
		   Last, the data is shuffled with a matrix on to the next layer.
		   Each node i from this layer effects each node j in the
		   next layer with some coefficient mx[i][j].

		*/

		BSC1dxMatrixLayer(){ this->id = BSC1DX_MATRIX_LAYER_ID; }
		
		BSC1dxMatrixLayer(int32_t n_, int32_t m_, T zero_, T one_)
			: BSCMatrixLayer<T>(n_, m_, zero_, one_){	
			this->id = BSC1DX_MATRIX_LAYER_ID;
			this->init_config();
		}
		
		BSC1dxMatrixLayer(int32_t n_, T zero_ = (T)0, T one_ = (T)1)
			: BSCMatrixLayer<T>(n_, zero_, one_){	
			
			this->id = BSC1DX_MATRIX_LAYER_ID;
			this->init_config();
		}
		
		BSC1dxMatrixLayer(ifstream &get_in){
			this->variables_in(get_in);
		}

		~BSC1dxMatrixLayer(){}

		void init_config(){	
			this->configClar = {
				"matrix_change_speed:",
				"bias_change_speed:",
				"sensetivity_change_speed:",
				"x-axis_compression:"
			};
			this->config = {(T)0.01, (T)0.01, (T)0.001, (T)1};
		}

		void project_next(Layer<T> *next){
			
			next->set_vector_all(this->zero);
		
			compress::div_x<T>(this->v, this->slope, {this->config[3]});

			for(int32_t i=0; i<this->n; i++){
				
				this->v[i] = (this->v[i]+this->bias[i])*this->sens[i];
				this->ucv[i] = this->v[i];

				for(int32_t j=0; j<this->m; j++){
					next->v[j] += this->mx[i][j]*this->v[i];
				}
			}
		}
};

#endif
