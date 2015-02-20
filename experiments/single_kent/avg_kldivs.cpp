#include <iostream>
#include <fstream>
#include <sstream>
#include <vector> 
#include <cstdlib>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::program_options;
using namespace boost::filesystem;

typedef std::vector<double> Vector;

#define NUM_METHODS 4 
#define MOMENT 0 
#define MLE 1
#define MAP 2
#define MML 3 

std::vector<Vector> load_table(string &file_name, int D)
{
  std::vector<Vector> sample;
  ifstream file(file_name.c_str());
  string line;
  Vector numbers(D,0);
  int i;
  while (getline(file,line)) {
    boost::char_separator<char> sep(" \t");
    boost::tokenizer<boost::char_separator<char> > tokens(line,sep);
    i = 0;
    BOOST_FOREACH(const string &t, tokens) {
      istringstream iss(t);
      double x;
      iss >> x;
      numbers[i++] = x;
    }
    sample.push_back(numbers);
  }
  file.close();
  return sample;
}

int partition(Vector &list, std::vector<int> &index, int left, int right)
{
	double temp,pivotPoint = list[right];
	int storeIndex = left,temp_i;
	for(int i=left; i<right; i++) {
		if(list[i] < pivotPoint) {
			temp = list[i];
			list[i] = list[storeIndex];
			list[storeIndex] = temp;
			temp_i = index[i];
			index[i] = index[storeIndex];
			index[storeIndex] = temp_i;
			storeIndex += 1;	
		}
	}
	temp = list[storeIndex];
	list[storeIndex] = list[right];
	list[right] = temp;
	temp_i = index[storeIndex];
	index[storeIndex] = index[right];
	index[right] = temp_i;
	return storeIndex;
}

void quicksort(Vector &list, std::vector<int> &index, int left, int right)
{
	if(left < right)
	{
		int pivotNewIndex = partition(list,index,left,right);
		quicksort(list,index,left,pivotNewIndex-1);
		quicksort(list,index,pivotNewIndex+1,right);
	}
}

Vector sort(Vector &list)
{
  int num_samples = list.size();
	Vector sortedList(list);
  std::vector<int> index(num_samples,0);
	for(int i=0; i<num_samples; i++) {
			index[i] = i;
  }
	quicksort(sortedList,index,0,num_samples-1);
  return sortedList;
}

std::vector<Vector> flip(std::vector<Vector> &table)
{
  int num_rows = table.size();
  Vector empty_vector(num_rows,0);
  int num_cols = table[0].size();
  std::vector<Vector> inverted_table(num_cols,empty_vector);
  for (int i=0; i<num_cols; i++) {
    for (int j=0; j<num_rows; j++) {
      inverted_table[i][j] = table[j][i];
    }
  }
  return inverted_table;
}

double computeMedian(Vector &list)
{
  Vector sorted_list = sort(list);
  int n = sorted_list.size();
  if (n % 2 == 1) {
    return sorted_list[n/2];
  } else {
    return (sorted_list[n/2-1]+sorted_list[n/2])/2;
  }
}

Vector computeMedians(std::vector<Vector> &table)
{
  std::vector<Vector> inverted_table = flip(table);
  int num_cols = table[0].size();
  Vector medians(num_cols,0);
  for (int i=0; i<num_cols; i++) {
    medians[i] = computeMedian(inverted_table[i]);
  }
  return medians;
}

double computeMean(Vector &list)
{
  double sum = 0;
  for (int i=0; i<list.size(); i++) {
    sum += list[i];
  }
  return sum / (double)list.size();
}

Vector computeMeans(std::vector<Vector> &table)
{
  std::vector<Vector> inverted_table = flip(table);
  int num_cols = table[0].size();
  Vector means(num_cols,0);
  for (int i=0; i<num_cols; i++) {
    means[i] = computeMean(inverted_table[i]);
  }
  return means;
}

double computeVariance(Vector &list)
{
  double mean = computeMean(list);
  double sum = 0;
  for (int i=0; i<list.size(); i++) {
    sum += (list[i]-mean) * (list[i]-mean);
  }
  return sum / (double) (list.size()-1);
}

int minimumIndex(Vector &values)
{
  int min_index = 0;
  double min_val = values[0];
  for (int i=1; i<values.size(); i++) { 
    if (values[i] <= min_val) {
      min_index = i;
      min_val = values[i];
    } // if()
  } // for()
  return min_index;
}

int maximumIndex(Vector &values)
{
  int max_index = 0;
  double max_val = values[0];
  for (int i=1; i<values.size(); i++) { 
    if (values[i] > max_val) {
      max_index = i;
      max_val = values[i];
    }
  }
  return max_index;
}

Vector computeMeans(ostream &out, std::vector<Vector> &p_est_all)
{
  Vector means = computeMeans(p_est_all);
  for (int i=0; i<p_est_all[0].size(); i++) {
    out << fixed << scientific << setprecision(6) << means[i] << "\t";
  }
  out << endl;
  return means;
}

void computeWins(
  ostream &out, std::vector<Vector> &values
) {
  int num_elements = values.size();
  std::vector<int> wins(values[0].size(),0);

  int min_index;
  for (int i=0; i<num_elements; i++) {
    min_index = minimumIndex(values[i]);
    wins[min_index]++;
  }
  double percent_wins;
  for (int j=0; j<wins.size(); j++) {
    percent_wins = wins[j] * 100.0 / num_elements;
    out << fixed << setw(10) << setprecision(2) << percent_wins;
  }
  out << endl;
}

struct Parameters
{
  int N;
  double kappa,eccentricity;
  int quantity;
};

struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;
  string quantity;

  options_description desc("Allowed options");
  desc.add_options()
       ("n",value<int>(&parameters.N),"value to be scaled")
       ("kappa",value<double>(&parameters.kappa),"kappa")
       ("eccentricity",value<double>(&parameters.eccentricity),"eccentricity")
       ("all",value<string>(&quantity),"all (kappa/eccentricity) ?")
  ;
  variables_map vm;
  store(command_line_parser(argc,argv).options(desc).run(),vm);
  notify(vm);

  if (vm.count("all")) {
    if (quantity.compare("k") == 0) { // all k fixed e
      parameters.quantity = 1;
    } else if (quantity.compare("e") == 0) { // all e fixed k
      parameters.quantity = 2;
    }
  } 
  return parameters;
}

// Usage: ./avg_kldivs --n 100 --kappa 10 --all e
// Usage: ./avg_kldivs --n 100 --eccentricity 0.1 --all k
int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  double kappa,beta,eccentricity;

  //string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_kappa_until_50/";
  //string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_uniform_prior/";
  //string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_vmf_prior/";
  //string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_beta_prior/";
  //string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_new2_prior/";
  string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_new32_prior/";
  //string n_str = "./N_" + boost::lexical_cast<string>(parameters.N) + "/";

  string current_dir,kappa_str,eccentricity_str;
  string kldivs_file;
  std::vector<Vector> kldivs_table;

  if (parameters.quantity == 1) { // all k fixed e
    eccentricity = parameters.eccentricity;
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << eccentricity;
    string eccentricity_str = sse.str();
    string avg_kldivs_file = n_str + "e_" + eccentricity_str + "_avg_kldivs.dat";
    string wins_kldivs_file = n_str + "e_" + eccentricity_str + "_wins_kldivs.dat";
    ofstream wins_kldivs(wins_kldivs_file.c_str());
    ofstream avg_kldivs(avg_kldivs_file.c_str());
    kappa = 10;
    while (kappa <= 100) {
      ostringstream ssk;
      ssk << fixed << setprecision(0);
      ssk << kappa;
      kappa_str = ssk.str();
      current_dir = n_str + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
      kldivs_file = current_dir + "kldivs";
      kldivs_table = load_table(kldivs_file,NUM_METHODS);
      avg_kldivs << fixed << setw(10) << setprecision(1) << eccentricity;
      avg_kldivs << fixed << setw(10) << kappa << "\t";
      computeMeans(avg_kldivs,kldivs_table);
      wins_kldivs << fixed << setw(10) << setprecision(1) << eccentricity;
      wins_kldivs << fixed << setw(10) << kappa << "\t";
      computeWins(wins_kldivs,kldivs_table);
      kappa += 10;
    } // while()
    wins_kldivs.close();
    avg_kldivs.close();
  } else if (parameters.quantity == 2) { // all e fixed k
    kappa = parameters.kappa;
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string avg_kldivs_file = n_str + "k_" + kappa_str + "_avg_kldivs.dat";
    string wins_kldivs_file = n_str + "k_" + kappa_str + "_wins_kldivs.dat";
    ofstream wins_kldivs(wins_kldivs_file.c_str());
    ofstream avg_kldivs(avg_kldivs_file.c_str());
    eccentricity = 0.1;
    while (eccentricity < 0.95) {
      ostringstream sse;
      sse << fixed << setprecision(1);
      sse << eccentricity;
      eccentricity_str = sse.str();
      current_dir = n_str + "k_" + kappa_str + "_e_" + eccentricity_str + "/";
      kldivs_file = current_dir + "kldivs";
      kldivs_table = load_table(kldivs_file,NUM_METHODS);
      avg_kldivs << fixed << setw(10) << kappa;
      avg_kldivs << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
      computeMeans(avg_kldivs,kldivs_table);
      wins_kldivs << fixed << setw(10) << kappa;
      wins_kldivs << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
      computeWins(wins_kldivs,kldivs_table);
      eccentricity += 0.1;
    } // while() 
    wins_kldivs.close();
    avg_kldivs.close();
  } // if()
}

