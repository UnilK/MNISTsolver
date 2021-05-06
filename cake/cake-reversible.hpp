#ifndef CAKE_REVERSIBLE_HPP_
#define CAKE_REVERSIBLE_HPP_

#include <vector>
#include <fstream>

#include "layer/base.hpp"
#include "layer/base-reversible.hpp"

using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;

const int32_t REVERSIBLE_CAKE_ID = 0x00010000;

template<class T> class ReversibleCake{

	/*

	   A cake solely consisting of reversible layers.

	   This class provides a simple API for running
	   data through the network of layers and changing
	   the network based on the results.

	*/

	protected:

		int32_t n = 0;
		T zero;
		vector<ReversibleLayer<T>*>  layer;

	public:

		int32_t id = REVERSIBLE_CAKE_ID;

		ReversibleCake(){}

		ReversibleCake(T zero_){
			this->zero = zero_;
			this->n = 0;
		}

		~ReversibleCake(){
			for(auto i : this->layer) delete i;
		}

		// The construction of a cake runs in 2 parts.

		// 1. add the desired layers
		void add_layer(ReversibleLayer<T> *new_layer){
			this->layer.push_back(new_layer);
			this->n++;
		}

		// 2. connect the layers to a cake.
		void connect_layers(){
			for(int32_t i=0; i<n-1; i++) this->layer[i]->connect_next(this->layer[i+1]->n);
			this->layer[n-1]->connect_next(this->layer[n-1]->n);
		}

		// randomize the variables used in the layers, varible = random_func().
		// Note that the return value of random_func doesn't have to be random.
		void random_variables(T (*random_func)(void)){
			for(auto i : this->layer) i->random_variables(random_func);
		}

		/*
		   Doing changes based on just one case doesn't really make
		   sense. It's better to take in a large batch of inputs,
		   average the desired changes and then adjust.

		   general structure for training:

		   cake.zero_changes();
		   
		   for(int i=0; i<10; i++){
				vector<T> processResult = cake.process(test_data[i]);
				vector<T> feedback = compare(processResult, desired_result[i]);
				cake.evaluate(feedback);
		   }

		   cake.downscale_changes(10);
		   cake.adjust();

		*/

		// changes should be reset using this function
		// between training batches.
		void zero_changes(){
			for(auto i : this->layer) i->zero_changes();
		}

		// For averaging the accumulated changes.
		void downscale_changes(T down){
			for(auto i : this->layer) i->downscale_changes(down);
		}

		// Runs the input data through the cake.
		vector<T> process(const vector<T> data_in){
			this->layer[0]->set_vector_values(data_in);	
			for(int32_t i=0; i<n-1; i++) this->layer[i]->project_next(this->layer[i+1]);
			this->layer[n-1]->project_next(this->layer[n-1]);
			return this->layer[n-1]->get_vector();
		}

		// Evaluates how successfull the last process run was
		// and accumulates the desired changes
		void evaluate(vector<T> &feedback){
			this->layer[n-1]->evaluate(feedback);
			for(int32_t i=n-2; i>=0; i--){
				this->layer[i]->evaluate(this->layer[i+1]->get_vector_changes());
			}
		}

		// apply the desired changes
		void adjust(){
			for(auto i : this->layer) i->adjust();
		}

		/*
		   Modifying the configuration of the cake:

		   cake.write_config();

		   // change the file config.ckc as you wish

		   cake.read_config();

		*/

		void write_config(){
			ofstream get_out("config.ckc");

			get_out << "CONFIG " << this->id << ' ' << this->n << '\n'; 
			
			int32_t counter = 0;
			for(auto i : this->layer){
				get_out << "Layer" << counter << ' ';
				counter++;
				i->config_out(get_out);
			}
			get_out.close();
		}

		void read_config(){
			
			ifstream get_in("config.ckc");
			
			int32_t idt, nt;
			string clarification;

			get_in >> clarification >> idt >> nt;

			if(idt == this->id && nt == this->n){
				for(auto i : this->layer){
					get_in >> clarification >> idt;
					if(i->id == idt) i->config_in(get_in);
					else break;
				}
			}
			get_in.close();
		}

		/*
		   Only a writing function is provided in this class,
		   as the layers are conceptually only loosely related
		   to the cake. A read function is implemented
		   based on the layers that are used.
		*/

		void write_file(std::string filename){
			
			ofstream get_out(filename);

			get_out << this->id << ' ' << this->n << '\n';
			for(auto i : this->layer) i->variables_out(get_out);

			get_out.close();
		}

};

#endif
