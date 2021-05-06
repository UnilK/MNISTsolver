#ifndef BASE_REVERSIBLE_LAYER_HPP_
#define BASE_REVERSIBLE_LAYER_HPP_

#include <vector>
#include <algorithm>
#include <fstream>

using std::vector;
using std::ifstream;
using std::ofstream;

const int32_t REVERSIBLE_LAYER_ID = 0x0001;

template<class T> class ReversibleLayer : public Layer<T>{

	/*
	   Reverse engineerable layers have the property that
	   given some feedback, the layer can figure
	   out itself how it should change.

	   The desired changes are stored in variables
	   marked with a capital C suffix.
	*/

	protected:

		vector<T> vC;

	public:

		ReversibleLayer(){ this->id = REVERSIBLE_LAYER_ID; }
		
		ReversibleLayer(int32_t n_, int32_t m_, T zero_) : Layer<T>(n_, m_, zero_){	
			this->vC.resize(n_, zero_);
			this->id = REVERSIBLE_LAYER_ID;
		}
		
		ReversibleLayer(int32_t n_, T zero_ = (T)0) : Layer<T>(n_, zero_){
			this->vC.resize(n_, zero_);
			this->id = REVERSIBLE_LAYER_ID;
		}

		ReversibleLayer(ifstream &get_in){
			this->variables_in(get_in);
		}
		
		virtual ~ReversibleLayer(){}

		virtual void variables_in(ifstream &get_in){

			if(!get_in.good()) return;

			get_in >> this->id >> this->n >> this->m >> this->zero;
			
			this->config_in(get_in);
			
			this->v.resize(this->n, this->zero);
			this->connect_next(this->m);
			this->vC.resize(this->n, this->zero);
		}

		virtual void variables_out(ofstream &get_out){
			get_out << this->id << ' ' << this->id << '\n';
			get_out << this->n << ' ' << this->m << ' ' << this->zero << '\n';
			
			this->config_out(get_out, 0);
		}

		// set all variables to a specific value
		virtual void set_variables(){}

		// for taking the average change over multiple training cases
		virtual void random_variables(T (*random_func)(void)){}

		void set_vector_changes(T val){
			for(T &i : this->vC) i = val;
		}

		vector<T> get_vector_changes(){
			return this->vC;
		}

		virtual void downscale_changes(T down){}
		
		// for resetting the changes between training batches
		virtual void zero_changes(){
			for(T &i : this->vC) i = this->zero;
		}

		// accumulate desired changes from feedback.
		virtual void evaluate(vector<T> feedback){
			for(int32_t i=0; i<std::min(this->n, this->m); i++){
				this->vC[i] = feedback[i];
			}
		}

		// desired changes are implemented.
		virtual void adjust(){}

};

#endif
