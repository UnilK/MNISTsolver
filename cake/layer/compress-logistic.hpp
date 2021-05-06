#ifndef COMPRESS_LOGISTIC_LAYER_HPP_
#define COMPRESS_LOGISTIC_LAYER_HPP_

#include <fstream>
#include <vector>
#include <algorithm>

#include "base.hpp"
#include "base-reversible.hpp"
#include "../func/compressions.hpp"

using std::vector;
using std::ifstream;
using std::ofstream;

const int32_t C_LOGISTIC_LAYER_ID = 0x0023;

template<class T> class CLogisticLayer: public ReversibleLayer<T>{

	protected:

		vector<T> slope;

	public:

		/*
		   The values get compressed to the range [0, 1] with:
		   
			v[i] = 1/(1+exp(-v[i]*config[0]))

		*/
		
		CLogisticLayer(){ this->id = C_LOGISTIC_LAYER_ID; }
		
		CLogisticLayer(int32_t n_, int32_t m_, T zero_) : ReversibleLayer<T>(n_, m_, zero_){
			this->id = C_LOGISTIC_LAYER_ID;
			this->slope.resize(n_);
			this->init_config();
		}
		
		CLogisticLayer(int32_t n_, T zero_ = (T)0) : ReversibleLayer<T>(n_, zero_){
			this->id = C_LOGISTIC_LAYER_ID;
			this->slope.resize(n_);
			this->init_config();
		}
		
		CLogisticLayer(ifstream &get_in){
			this->variables_in(get_in);
		}
		
		~CLogisticLayer(){}

		void init_config(){
			this->configClar = {"x-axis_compression:"};
			this->config = {(T)1};
		}

		void project_next(Layer<T> *next){

			compress::logistic<T>(this->v, this->slope, {this->config[0]});

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
			
			/*

			   Just taking an inverse from the desired change
			   doesn't really make sense. Inspecting how the
			   function behaves in the area feels a bit more natural.

			   The slope of the function at some point x is:

			   d/dx(1-0.5/x) = 0.5/x^2

			   I'll use that to do something.

			   Yeah I have no Idea how this works.
			   Feedback/slope makes more sense to me,
			   but it quickly results in unstable values.
			   Feedback*slope should at least make the values stable,
			   though I fear a "beached whale" effect of sorts when
			   the values get large - once it gets up there,
			   there's no bringing it back.

			   Actually, now that I've thought more about this,
			   nothing really states that the values going
			   into a particluar node would always be large.
			   So the slope can be interpreted as:

			   Steep slope: I'm volatile, you get big changes from me!
			   Flat slope: don't bother changing me, it won't do anything anyways.

			   An illustration of the intended shape of the function:

			                  .......
			               ...
			              .
			           ...
			   ........


			*/


			for(int32_t i=0; i<std::min(this->n, this->m); i++){
				this->vC[i] = this->slope[i]*feedback[i];
			}
		}
};

#endif
