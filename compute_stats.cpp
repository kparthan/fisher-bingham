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

std::vector<Vector> load_matrix(string &file_name, int D)
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

int minimumIndex(int ignore_map, Vector &values)
{
  int min_index = 0;
  double min_val = values[0];
  if (ignore_map) {
    for (int i=1; i<values.size(); i++) { 
      if (i != MAP) {
        if (values[i] <= min_val) {
          min_index = i;
          min_val = values[i];
        } // if()
      } 
    } // for()
  } else { // don't ignore_map
    for (int i=1; i<values.size(); i++) { 
      if (values[i] <= min_val) {
        min_index = i;
        min_val = values[i];
      } // if()
    } // for()
  } // if(ignore_map)
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

Vector computeEstimateMedians(ostream &out, std::vector<Vector> &p_est_all)
{
  Vector medians = computeMedians(p_est_all);
  for (int i=0; i<p_est_all[0].size(); i++) {
    out << scientific << medians[i] << "\t";
  }
  out << endl;
  return medians;
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

void computeMeanSquaredError(
  ostream &out, double p, std::vector<Vector> &p_est_all
) {
  int num_elements = p_est_all.size();

  Vector error(p_est_all[0].size(),0);
  for (int i=0; i<num_elements; i++) {
    for (int j=0; j<p_est_all[0].size(); j++) {
      error[j] += (p - p_est_all[i][j]) * (p - p_est_all[i][j]);
    }
  }
  for (int j=0; j<p_est_all[0].size(); j++) {
    error[j] /= num_elements;
    out << fixed << scientific << setprecision(6) << error[j] << "\t";
  }
  out << endl;
}

void computeWins(
  int ignore_map, ostream &out, std::vector<Vector> &values
) {
  int num_elements = values.size();
  std::vector<int> wins(values[0].size(),0);

  int min_index;
  for (int i=0; i<num_elements; i++) {
    min_index = minimumIndex(ignore_map,values[i]);
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
  int ignore_map;
};

struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;

  options_description desc("Allowed options");
  desc.add_options()
       ("n",value<int>(&parameters.N),"value to be scaled")
       ("ignore_map","to ignore map estimate or not ?")
  ;
  variables_map vm;
  store(command_line_parser(argc,argv).options(desc).run(),vm);
  notify(vm);

  if (vm.count("ignore_map")) {
    parameters.ignore_map = 1;
  } else {
    parameters.ignore_map = 0;
  }

  return parameters;
}

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  double kappa,beta,eccentricity;

  //string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_uniform_prior/";
  //string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_vmf_prior/";
  string n_str = "N_" + boost::lexical_cast<string>(parameters.N) + "_beta_prior/";
  string parent_dir = "experiments/single_kent/" + n_str + "/";
  string current_dir,kappa_str,eccentricity_str;
  string kappas_file,betas_file,negloglike_file,kldivs_file,msglens_file;
  std::vector<Vector> kappas_table,betas_table,negloglike_table,kldivs_table,msglens_table;

  string wins_kldivs_file = parent_dir + "wins_kldivs";
  string wins_negloglike_file = parent_dir + "wins_negloglike";
  string wins_msglens_file = parent_dir + "wins_msglens";
  string mse_kappa_file = parent_dir + "mse_kappas";
  string mse_beta_file = parent_dir + "mse_betas";
  string avg_kldivs_file = parent_dir + "avg_kldivs";
  string avg_negloglike_file = parent_dir + "avg_negloglike";
  string avg_msglens_file = parent_dir + "avg_msglens";

  ofstream wins_kldivs(wins_kldivs_file.c_str());
  ofstream wins_negloglike(wins_negloglike_file.c_str());
  ofstream wins_msglens(wins_msglens_file.c_str());
  ofstream mse_kappa(mse_kappa_file.c_str());
  ofstream mse_beta(mse_beta_file.c_str());
  ofstream avg_kldivs(avg_kldivs_file.c_str());
  ofstream avg_negloglike(avg_negloglike_file.c_str());
  ofstream avg_msglens(avg_msglens_file.c_str());

  kappa = 10;
  while (kappa <= 100) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    kappa_str = ssk.str();
    eccentricity = 0.1;
    while (eccentricity < 0.95) {
      beta = 0.5 * kappa * eccentricity;
      ostringstream sse;
      sse << fixed << setprecision(1);
      sse << eccentricity;
      eccentricity_str = sse.str();
      current_dir = parent_dir + "k_" + kappa_str + "_e_" + eccentricity_str + "/";

      kappas_file = current_dir + "kappas";
      kappas_table = load_matrix(kappas_file,NUM_METHODS);
      mse_kappa << fixed << setw(10) << kappa;
      mse_kappa << fixed << setw(10) << beta;
      mse_kappa << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
      computeMeanSquaredError(mse_kappa,kappa,kappas_table);

      betas_file = current_dir + "betas";
      betas_table = load_matrix(betas_file,NUM_METHODS);
      mse_beta << fixed << setw(10) << kappa;
      mse_beta << fixed << setw(10) << beta;
      mse_beta << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
      computeMeanSquaredError(mse_beta,beta,betas_table);

      negloglike_file = current_dir + "negloglike";
      negloglike_table = load_matrix(negloglike_file,NUM_METHODS);
      avg_negloglike << fixed << setw(10) << kappa;
      avg_negloglike << fixed << setw(10) << beta;
      avg_negloglike << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
      computeMeans(avg_negloglike,negloglike_table);
      wins_negloglike << fixed << setw(10) << kappa;
      wins_negloglike << fixed << setw(10) << beta;
      wins_negloglike << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
      computeWins(parameters.ignore_map,wins_negloglike,negloglike_table);

      kldivs_file = current_dir + "kldivs";
      kldivs_table = load_matrix(kldivs_file,NUM_METHODS);
      avg_kldivs << fixed << setw(10) << kappa;
      avg_kldivs << fixed << setw(10) << beta;
      avg_kldivs << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
      computeMeans(avg_kldivs,kldivs_table);
      wins_kldivs << fixed << setw(10) << kappa;
      wins_kldivs << fixed << setw(10) << beta;
      wins_kldivs << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
      computeWins(parameters.ignore_map,wins_kldivs,kldivs_table);

      msglens_file = current_dir + "msglens";
      msglens_table = load_matrix(msglens_file,NUM_METHODS);
      avg_msglens << fixed << setw(10) << kappa;
      avg_msglens << fixed << setw(10) << beta;
      avg_msglens << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
      computeMeans(avg_msglens,msglens_table);
      wins_msglens << fixed << setw(10) << kappa;
      wins_msglens << fixed << setw(10) << beta;
      wins_msglens << fixed << setw(10) << setprecision(1) << eccentricity << "\t";
      computeWins(parameters.ignore_map,wins_msglens,msglens_table);
      
      eccentricity += 0.1;
    } // eccentricity
    kappa += 10;
  } // kappa
  wins_kldivs.close();
  wins_negloglike.close();
  wins_msglens.close();
  mse_kappa.close();
  mse_beta.close();
  avg_kldivs.close();
  avg_negloglike.close();
  avg_msglens.close();
}

