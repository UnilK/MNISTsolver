#ifndef BASE_LAYER_HPP_
#define BASE_LAYER_HPP_

#include <vector>
#include <algorithm>
#include <fstream>

using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;

const int32_t BASE_LAYER_ID = 0x0000;

template<class T> class Layer{

	/*
	   vector v holds the data nodes (numbers) of the layer.
	   config holds some configuration details. What the values in
	   this array mean depends on the desires of the derived class.
	   config_clar is a list of clarifications on the variables in config.
	   n is the size of v.
	   m is the size of the next layer.
	   id is an (id)entifying integer (preferably unique for every layer class).
	   "zero" is the default value in the array.

	   Most of the functions in this class are
	   meant to be overridden by it's child classes,
	   you can tell them apart from the "virtual" keyword.

	   The Layer class itself is meant to be used
	   for just calculating the result for some
	   input, not training.

	   The following should work with most classes
	   derived from Layer:

	   int n=3, x=4;
	   float z = 0.0;
	   vector<Layer<float>*> cake;
	   cake.push_back(new DerivedLayer0<type>(x, x, z));
	   cake.push_back(new DerivedLayer1<type>(x, x, z));
	   cake.push_back(new DerivedLayer2<type>(x, x, z));

	   for(int i=0; i<n-1; i++) cake[i]->project_next(cake[i+1]);
	   cake[n-1]->project_next(cake[n-1]);

	*/

	protected:

		int32_t m=0;
		T zero;
		vector<T> config;
		vector<string> configClar;

	public:
		
		int32_t n=0, id=BASE_LAYER_ID;
		vector<T> v;

		Layer(){ this->id = BASE_LAYER_ID; }
		
		// layers initialized with n and m will (and must) be
		// connected to the next layer from the get-go.
		Layer(int32_t n_, int32_t m_, T zero_){
			this->n = n_;
			this->zero = zero_;
			this->v.resize(n_, zero_);
			this->connect_next(m_);
			this->id = BASE_LAYER_ID;
		}
		
		// layers initialized with just n must be connected later on.
		Layer(int32_t n_, T zero_=(T)0){
			this->n = n_;
			this->zero = zero_;
			this->v.resize(n_, zero_);
			this->m = 0;
			this->id = BASE_LAYER_ID;
		}

		// initializes layer from a file with n and m.
		Layer(ifstream &get_in){
			this->variables_in(get_in);
		}
		
		virtual ~Layer(){}

		void set_vector_all(T val){
			for(T &i : v) i = val;
		}

		void set_vector_values(vector<T> v_){
			this->v = v_;
		}

		vector<T> get_vector(){
			return this->v;
		}

		virtual void connect_next(int32_t m_){
			this->m = m_;
		}
		
		// the base layer just copies the data to the next layer
		virtual void project_next(Layer<T> *next){
			for(int32_t i=0; i<std::min(this->n, this->m); i++){
				next->v[i] = this->v[i];
			}
		}
		
		/*
		   These functions are for changing things about the class
		   when it's being trained. Like how fast variables change & such.
		   The format is:
		   config_size
		   clarification variable
		   clarification variable
		   ...
		*/
		void config_in(ifstream &get_in){

			if(!get_in.good()) return;

			int32_t config_size;
			get_in >> config_size;
			this->config.resize(config_size);
			this->configClar.resize(config_size);
			for(int32_t i=0; i<config_size; i++){
				get_in >> this->configClar[i] >> this->config[i];
			}
		}
		
		void config_out(ofstream &get_out, bool write_id=1){
			/*
			   The id is written, but read elsewhere.
			   - it's too late to ask what sort of class
			   one would like to configure if a specific class
			   has already been chosen!
			*/
			if(write_id) get_out << this->id << ' ';
			get_out << config.size() << '\n';
			for(int32_t i=0; i<(int32_t)config.size(); i++){
				get_out << configClar[i] << ' ' << config[i] << '\n';
			}
		}

		// These are for saving & loading layers
		virtual void variables_in(ifstream &get_in){

			if(!get_in.good()) return;

			get_in >> this->id >> this->n >> this->m >> this->zero;
			this->config_in(get_in);
			this->v.resize(this->n, this->zero);
			this->connect_next(this->m);
		}

		virtual void variables_out(ofstream &get_out){
			/*
			   Same thing as with the config, the id is read elsewhere
			   first. In this case it's convenient to leave a copy to be read
			   here too.
			*/
			get_out << this->id << ' ' << this->id << '\n';
			get_out << this->n << ' ' << this->m << ' ' << this->zero << '\n';
			this->config_out(get_out, 0);

		}
};

#endif
