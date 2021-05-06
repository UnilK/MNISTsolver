#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <time.h>
#include <filesystem>

#include "cake/func/randoms.hpp"
#include "cake/func/fft.hpp"

#include "cake/cake-reversible.hpp"

#include "cake/layer/base.hpp"
#include "cake/layer/base-reversible.hpp"
#include "cake/layer/matrix.hpp"
#include "cake/layer/compress-1dx.hpp"
#include "cake/layer/c1dx-matrix.hpp"
#include "cake/layer/compress-1dxp2.hpp"
#include "cake/layer/compress-logistic.hpp"
#include "cake/layer/c1dxp2-matrix.hpp"
#include "cake/layer/BSC-matrix.hpp"
#include "cake/layer/BSC1dx-matrix.hpp"
#include "cake/layer/convolution.hpp"
#include "cake/layer/sparse-convolution.hpp"

using std::cin;
using std::cout;
using std::string;
using std::vector;
using std::pair;
using std::ifstream;

/*

   this is and ad hoc test platform for the MNIST
   database. For some sensible code, see the cake library.

*/

ReversibleCake<float> *solution = new ReversibleCake<float>(0.0);
FFT<float> *fft = new FFT<float>();

void read_mnist_cake(string, ReversibleCake<float>*&);

class TrainProtocol{
	
	public:

		int32_t train_size = 0, test_size = 0;
		int32_t image_height = 0, image_width = 0, image_size = 0;

		ReversibleCake<float> *trainee = NULL;
		vector<pair<vector<uint8_t>, int32_t> > train_data, test_data;
		
		TrainProtocol(){
			train_size = 0;
			test_size = 0;
		}

		TrainProtocol(ReversibleCake<float> *trainee_){
			trainee = trainee_;
			train_size = 0;
			test_size = 0;
		}

		TrainProtocol(string filepath, ReversibleCake<float> *trainee_){
			trainee = trainee_;
			data_from_file(filepath);
		}

		void data_from_file(string filepath){

			ifstream header_in(filepath);

			string directory;
			for(int32_t i=0; i<(int32_t)filepath.size(); i++){
				if(filepath[i] == '/') directory = directory = filepath.substr(0, i+1);
			}
			
			string training_images, training_labels, test_images, test_labels;
			
			header_in >> training_images >> training_labels;
			header_in >> test_images >> test_labels;

			header_in.close();

			read_images(directory+training_images, train_data);
			read_labels(directory+training_labels, train_data);
			train_size = train_data.size();
			
			read_images(directory+test_images, test_data);
			read_labels(directory+test_labels, test_data);
			test_size = test_data.size();
		}

		void flip_bytes(char *reg, int32_t len){
			for(int32_t i=0; i<len/2; i++) std::swap(reg[i], reg[len-i-1]);
		}

		int32_t get_little_int(char *reg){
			int32_t num;
			flip_bytes(reg, 4);
			memcpy(&num, reg, 4);	
			return num;
		}

		int32_t read_little_int(ifstream &data_in){
			char reg[4];
			data_in.read(reg, 4);
			return get_little_int(reg);
		}

		void read_images(string filename, vector<pair<vector<uint8_t>, int32_t> > &data_out){

			ifstream data_in(filename);

			read_little_int(data_in); // this contains unnecessary information
			int32_t data_size = read_little_int(data_in);
			
			if(!data_in.good()) return;

			data_out.resize(data_size);

			image_height = read_little_int(data_in);
			image_width = read_little_int(data_in);
			image_size = image_height*image_width;
			
			if(!data_in.good()) return;

			for(int32_t i=0; i<data_size; i++){
				char image[image_size];
				data_in.read(image, image_size);
				data_out[i].first.resize(image_size);
				for(int32_t j=0; j<image_size; j++){
					data_out[i].first[j] = (uint8_t)image[j];
				}
			}

			data_in.close();
		}
		
		void read_labels(string filename, vector<pair<vector<uint8_t>, int32_t> > &data_out){

			ifstream data_in(filename);

			read_little_int(data_in); // same unnecessary information
			int32_t data_size = read_little_int(data_in);
			
			if(!data_in.good()) return;

			data_out.resize(data_size);

			for(int32_t i=0; i<data_size; i++){
				char label;
				data_in.read(&label, 1);
				data_out[i].second = (int32_t)label;
			}

			data_in.close();
		}

		vector<float> tofloat(vector<uint8_t> &data_in){
			vector<float> ret(data_in.size());
			for(int32_t i=0; i<(int32_t)ret.size(); i++) ret[i] = (float)((int32_t)data_in[i]);
			return ret;
		}
		
		void train_batch(int32_t size){
		
			trainee->zero_changes();

			for(int32_t i=0; i<size; i++){
				int32_t sample = rand()%train_size;

				vector<float> feedback = trainee->process(tofloat(train_data[sample].first));
				
				for(float &i : feedback) i = -i;

				feedback[train_data[sample].second] += 1.0;

				trainee->evaluate(feedback);
			}

			trainee->downscale_changes((float)size);
			trainee->adjust();

		}

		void train_batches(int32_t amount, int32_t size){
			for(int32_t i=0; i<amount; i++){
				train_batch(size);
			}
		}


		int32_t test(){

			int32_t score = 0;

			for(int32_t i=0; i<test_size; i++){
				vector<float> result = trainee->process(tofloat(test_data[i].first));
				int32_t ans = 0;
				float max = -1e9;
				for(int32_t j=0; j<10; j++){
					if(max < result[j]){
						ans = j;
						max = result[j];
					}
				}
				if(ans == test_data[i].second) score++;
			}

			return score;
		}

		void train_program(string name, int32_t amount, int32_t size){
			

			int32_t count = 0, best = -1;

			string dir = "saves/"+name, prevBest;

			std::filesystem::create_directory(dir);

			dir += "/";

			ofstream confo(dir+"train_state");
			confo << "continue\n";
			confo.close();

			while(1){
				
				train_batches(amount, size);
				int32_t score = test();
				
				if(score > best){
					best = score;
					prevBest = dir+"#"+std::to_string(count)+"-"+std::to_string(best);
					trainee->write_file(prevBest);
				} else {
					read_mnist_cake(prevBest, trainee);
					solution = trainee;
				}

				ifstream confi(dir+"train_state");
				string state;
				confi >> state;
				confi.close();

				count++;

				if(state != "continue") break;
			}
			

		}
};



void read_mnist_cake(string filename, ReversibleCake<float>* &cake){
	
	ifstream get_in(filename);

	if(!get_in.good()) return;

	delete cake;
	cake = new ReversibleCake<float>(0.0);

	int32_t n;
	get_in >> cake->id >> n;

	for(int i=0; i<n; i++){
		int32_t id;
		get_in >> id;
		
		ReversibleLayer<float> *layer = NULL;

		switch(id){
			case REVERSIBLE_LAYER_ID:
				layer = new ReversibleLayer<float>(get_in);
				break;
			case MATRIX_LAYER_ID:
				layer = new MatrixLayer<float>(get_in);
				break;
			case BSC_MATRIX_LAYER_ID:
				layer = new BSCMatrixLayer<float>(get_in);
				break;
			case BSC1DX_MATRIX_LAYER_ID:
				layer = new BSC1dxMatrixLayer<float>(get_in);
				break;
			case C_1DX_LAYER_ID:
				layer = new C1dxLayer<float>(get_in);
				break;
			case C_1DX_MATRIX_LAYER_ID:
				layer = new C1dxMatrixLayer<float>(get_in);
				break;
			case C_1DXP2_LAYER_ID:
				layer = new C1dxp2Layer<float>(get_in);
				break;
			case C_1DXP2_MATRIX_LAYER_ID:
				layer = new C1dxp2MatrixLayer<float>(get_in);
				break;
			case C_LOGISTIC_LAYER_ID:
				layer = new CLogisticLayer<float>(get_in);
				break;
			case CONVOLUTION_LAYER_ID:
				layer = new ConvolutionLayer<float>(get_in, fft);
				break;
			case SPARSE_CONVOLUTION_LAYER_ID:
				layer = new SparseConvolutionLayer<float>(get_in, fft);
				break;
		}

		if(layer != NULL) cake->add_layer(layer);
	}
	
	cake->connect_layers();
}

int main(){

	srand(time(0));

	solution->add_layer(new SparseConvolutionLayer<float>(784, fft));
	solution->add_layer(new C1dxMatrixLayer<float>(240));
	solution->add_layer(new C1dxMatrixLayer<float>(16));
	solution->add_layer(new C1dxMatrixLayer<float>(16));
	solution->add_layer(new C1dxp2Layer<float>(10));

	solution->connect_layers();
	solution->random_variables(randoms::float24<float>);

	TrainProtocol protocol("data/MNIST/input_files", solution);

	while(1){
		
		string inst;
		cin >> inst;

		if(inst == "data"){

			string filename;
			cin >> filename;

			protocol.data_from_file(filename);
			cout << "done\n";

		} else if(inst == "train"){

			int32_t amount, size;
			cin >> amount >> size;

			protocol.train_batches(amount, size);

			cout << "done\n";
		} else if(inst == "program"){

			string name;
			int32_t amount, size;
			cin >> name >> amount >> size;

			protocol.train_program(name, amount, size);

			cout << "done\n";

		} else if(inst == "test"){

			cout << protocol.test() << '\n';

		} else if(inst == "config"){

			cin >> inst;
			if(inst == "in") solution->read_config();
			else if(inst == "out") solution->write_config();
			else cout << "unknown instruction, try \"help\" \n";

		} else if(inst == "save"){

			string filename;
			cin >> filename;

			solution->write_file("saves/"+filename);

			cout << "done\n";

		} else if(inst == "load"){

			string filename;
			cin >> filename;

			read_mnist_cake("saves/"+filename, solution);
			protocol.trainee = solution;

			cout << "done\n";

		} else if(inst == "help"){

			cout
				<< "\nsupported commands are:\n"
				<< "data filename(string)\n"
				<< "train amount(int) batch_size(int)\n"
				<< "test\n"
				<< "config in/out\n"
				<< "save filename(string)\n"
				<< "load filename(string)\n"
				<< "help (duh)\n"
				<< "exit\n\n";

		} else if(inst == "exit"){
			break;
		} else {
			cout << "unknown instruction, try \"help\" \n";
		}
	}

	return 0;

}
