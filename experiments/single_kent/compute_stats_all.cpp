#include <iostream>
#include <fstream>
#include <sstream>
#include <vector> 
#include <cstdlib>
#include <iomanip>
#include <sys/stat.h>

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

#define PSI M_PI/4.0
#define ALPHA M_PI/2.0
#define ETA M_PI/2.0

struct stat st = {0};

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

/* S **********************************************************************************/ 
void plot_script1_kldivs(string &kldivs_folder)
{
  string wins_kldivs_file = kldivs_folder + "wins_kldivs.dat";
  string script_file = kldivs_folder + "plot1_wins_kldivs.p";
  string plot_file = kldivs_folder + "plot1_wins_kldivs.eps";
  ofstream out(script_file.c_str());
  out << "set terminal postscript eps enhanced color\n\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "set grid y\n";
  out << "set style data histograms\n";
  out << "set style fill solid 1.0 noborder\n";
  out << "set ytics 10 nomirror\n";
  out << "set yrange [:100]\n";
  out << "set ylabel \"\% of wins\"\n";
  out << "set ytics 10\n\n"; 
  out << "plot \"" << wins_kldivs_file << "\" using 2 t \"MOM\", \\\n"
      << "\"\" using 3 t \"MLE\", \\\n"
      << "\"\" using 4 t \"MAP\", \\\n"
      << "\"\" using 5:xtic(1) t \"MML\"";
  out.close();
  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

void plot_script2_kldivs(string &kldivs_folder)
{
  string wins_kldivs_file = kldivs_folder + "wins_kldivs.dat";
  string script_file = kldivs_folder + "plot2_wins_kldivs.p";
  string plot_file = kldivs_folder + "plot2_wins_kldivs.eps";
  ofstream out(script_file.c_str());
  out << "set terminal postscript eps enhanced color\n\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "set style data linespoints\n";
  out << "set style fill solid 1.0 noborder\n";
  out << "set ytics 10 nomirror\n";
  out << "set yrange [:100]\n";
  out << "set ylabel \"\% of wins\"\n";
  out << "set ytics 10\n\n"; 
  out << "plot \"" << wins_kldivs_file << "\" using 1:2 t \"MOM\", \\\n"
      << "\"\" using 1:3 t \"MLE\", \\\n"
      << "\"\" using 1:4 t \"MAP\", \\\n"
      << "\"\" using 1:5 t \"MML\"";
  out.close();
  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

void plot_stack_wins_kldivs(string &kldivs_folder)
{
  string wins_kldivs_file = kldivs_folder + "wins_kldivs.dat";
  string script_file = kldivs_folder + "stack_wins_kldivs.p";
  string plot_file = kldivs_folder + "stack_wins_kldivs.eps";
  ofstream out(script_file.c_str());
  out << "set terminal postscript eps enhanced color\n\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "set key invert reverse left top\n\n";
  out << "set grid y\n";
  out << "set style data histograms\n";
  out << "set style histogram rowstacked\n";
  out << "set boxwidth 0.5\n";
  out << "set style fill solid 1.0 noborder\n";
  out << "set ytics 10 nomirror\n";
  out << "set yrange [:100]\n";
  out << "set ylabel \"\% of wins\"\n";
  out << "set ytics 10\n\n"; 
  out << "plot \"" << wins_kldivs_file << "\" using 2 t \"MOM\", \\\n"
      << "\"\" using 3 t \"MLE\", \\\n"
      << "\"\" using 4 t \"MAP\", \\\n"
      << "\"\" using 5:xtic(1) t \"MML\"";
  out.close();
  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

void boxplot_kldivs_fixed_ecc(string &n_str, string &ecc_str, string &kldivs_folder)
{
  string kldivs_file;
  string script_file = kldivs_folder + "boxplot_kldivs.p";
  string plot_file = kldivs_folder + "boxplot_kldivs.eps";

  ofstream out(script_file.c_str());
  out << "set terminal post eps color enhanced\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "box_width=0.20\n";
  out << "set style fill solid 0.25 noborder\n";
  out << "set style boxplot outliers pointtype 7\n";
  out << "set style data boxplot\n";
  out << "set boxwidth box_width #relative\n";
  out << "set pointsize 0.5\n";
  out << "unset key\n";
  out << "set border 2\n";
  out << "set xtics nomirror\n";
  out << "set ytics nomirror\n";
  out << "set xlabel \"{/Symbol k}\\n\"\n";
  out << "set ylabel \"KL-divergence\"\n";
  out << "set xlabel font \"Times-Roman, 25\"\n";
  out << "set ylabel font \"Times-Roman, 25\"\n";
  out << "set xtics font \"Times-Roman, 20\"\n";
  out << "set ytics font \"Times-Roman, 20\"\n";
  out << "set xtics (\"10\" 1, \"50\" 2, \"100\" 3) scale 0.0\n";
  out << "d_width=0.5*box_width\n\n";

  kldivs_file = n_str + "k_10_e_" + ecc_str + "/kldivs";
  out << "plot \"" << kldivs_file << "\" using ((1)-3*d_width):1 lt 1 lc rgb \"red\", \\\n"
      << "\"" << kldivs_file << "\" using ((1)-d_width):2 lt 1 lc rgb \"blue\", \\\n"
      << "\"" << kldivs_file << "\" using ((1)+d_width):2 lt 1 lc rgb \"dark-green\", \\\n"
      << "\"" << kldivs_file << "\" using ((1)+3*d_width):2 lt 1 lc rgb \"black\", \\\n";
  kldivs_file = n_str + "k_50_e_" + ecc_str + "/kldivs";
  out << "\"" << kldivs_file << "\" using ((2)-3*d_width):1 lt 1 lc rgb \"red\", \\\n"
      << "\"" << kldivs_file << "\" using ((2)-d_width):2 lt 1 lc rgb \"blue\", \\\n"
      << "\"" << kldivs_file << "\" using ((2)+d_width):2 lt 1 lc rgb \"dark-green\", \\\n"
      << "\"" << kldivs_file << "\" using ((2)+3*d_width):2 lt 1 lc rgb \"black\", \\\n";
  kldivs_file = n_str + "k_100_e_" + ecc_str + "/kldivs";
  out << "\"" << kldivs_file << "\" using ((3)-3*d_width):1 lt 1 lc rgb \"red\", \\\n"
      << "\"" << kldivs_file << "\" using ((3)-d_width):2 lt 1 lc rgb \"blue\", \\\n"
      << "\"" << kldivs_file << "\" using ((3)+d_width):2 lt 1 lc rgb \"dark-green\", \\\n"
      << "\"" << kldivs_file << "\" using ((3)+3*d_width):2 lt 1 lc rgb \"black\"\n";
  out.close();

  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

void boxplot_kldivs_fixed_kappa(string &n_str, string &kappa_str, string &kldivs_folder)
{
  string kldivs_file;
  string script_file = kldivs_folder + "boxplot_kldivs.p";
  string plot_file = kldivs_folder + "boxplot_kldivs.eps";

  ofstream out(script_file.c_str());
  out << "set terminal post eps color enhanced\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "box_width=0.20\n";
  out << "set style fill solid 0.25 noborder\n";
  out << "set style boxplot outliers pointtype 7\n";
  out << "set style data boxplot\n";
  out << "set boxwidth box_width #relative\n";
  out << "set pointsize 0.5\n";
  out << "unset key\n";
  out << "set border 2\n";
  out << "set xtics nomirror\n";
  out << "set ytics nomirror\n";
  out << "set xlabel \"{/Symbol k}\\n\"\n";
  out << "set ylabel \"KL-divergence\"\n";
  out << "set xlabel font \"Times-Roman, 25\"\n";
  out << "set ylabel font \"Times-Roman, 25\"\n";
  out << "set xtics font \"Times-Roman, 20\"\n";
  out << "set ytics font \"Times-Roman, 20\"\n";
  out << "set xtics (\"0.1\" 1, \"0.5\" 2, \"0.9\" 3) scale 0.0\n";
  out << "d_width=0.5*box_width\n\n";

  kldivs_file = n_str + "k_" + kappa_str + "_e_0.1/kldivs";
  out << "plot \"" << kldivs_file << "\" using ((1)-3*d_width):1 lt 1 lc rgb \"red\", \\\n"
      << "\"" << kldivs_file << "\" using ((1)-d_width):2 lt 1 lc rgb \"blue\", \\\n"
      << "\"" << kldivs_file << "\" using ((1)+d_width):2 lt 1 lc rgb \"dark-green\", \\\n"
      << "\"" << kldivs_file << "\" using ((1)+3*d_width):2 lt 1 lc rgb \"black\", \\\n";
  kldivs_file = n_str + "k_" + kappa_str + "_e_0.5/kldivs";
  out << "\"" << kldivs_file << "\" using ((2)-3*d_width):1 lt 1 lc rgb \"red\", \\\n"
      << "\"" << kldivs_file << "\" using ((2)-d_width):2 lt 1 lc rgb \"blue\", \\\n"
      << "\"" << kldivs_file << "\" using ((2)+d_width):2 lt 1 lc rgb \"dark-green\", \\\n"
      << "\"" << kldivs_file << "\" using ((2)+3*d_width):2 lt 1 lc rgb \"black\", \\\n";
  kldivs_file = n_str + "k_" + kappa_str + "_e_0.9/kldivs";
  out << "\"" << kldivs_file << "\" using ((3)-3*d_width):1 lt 1 lc rgb \"red\", \\\n"
      << "\"" << kldivs_file << "\" using ((3)-d_width):2 lt 1 lc rgb \"blue\", \\\n"
      << "\"" << kldivs_file << "\" using ((3)+d_width):2 lt 1 lc rgb \"dark-green\", \\\n"
      << "\"" << kldivs_file << "\" using ((3)+3*d_width):2 lt 1 lc rgb \"black\"\n";
  out.close();

  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
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

// all k fixed e
void tabulate_eccentricity_kldivs(
  string &output_file, string &n_str, string &ecc_str
) {
  std::vector<Vector> kldivs;
  ofstream out(output_file.c_str());
  double kappa = 10;
  while (kappa < 101) {
    out << fixed << setw(10) << setprecision(0) << kappa;
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";
    string kldivs_file = current_dir + "kldivs";
    kldivs = load_table(kldivs_file,NUM_METHODS);
    computeWins(out,kldivs);
    kappa += 10;
  } // while()
  out.close();
}

// all e fixed k
void tabulate_kappa_kldivs(
  string &output_file, string &n_str, string &kappa_str
) {
  std::vector<Vector> kldivs;
  string kldivs_file;
  ofstream out(output_file.c_str());
  double ecc = 0.1;
  while (ecc < 0.95) {
    out << fixed << setw(10) << setprecision(1) << ecc;
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    string ecc_str = sse.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";
    kldivs_file = current_dir + "kldivs";
    kldivs = load_table(kldivs_file,NUM_METHODS);
    computeWins(out,kldivs);
    ecc += 0.1;
  } // while()
  out.close();
}

// all k fixed e
void compute_average_kldivs_fixed_ecc(
  string &n_str, string &ecc_str, string &kldivs_folder
) {
  std::vector<Vector> kldivs;
  string output_file = kldivs_folder + "kldivs_avg.dat";
  ofstream out(output_file.c_str());
  double kappa = 10;
  while (kappa < 101) {
    out << fixed << setw(10) << setprecision(0) << kappa << "\t\t";
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";
    string kldivs_file = current_dir + "kldivs";
    kldivs = load_table(kldivs_file,NUM_METHODS);
    Vector avg_kldivs = computeMeans(kldivs);
    for (int i=0; i<NUM_METHODS; i++) {
      out << fixed << scientific << setprecision(6) << avg_kldivs[i] << "\t\t";
    }
    out << endl;

    kappa += 10;
  } // while()
  out.close();
}

void plot_avg_kldivs(string &kldivs_folder)
{
  string avg_kldivs_file = kldivs_folder + "kldivs_avg.dat";
  string script_file = kldivs_folder + "kldivs_avg.p";
  string plot_file = kldivs_folder + "kldivs_avg.eps";
  ofstream out(script_file.c_str());
  out << "set terminal postscript eps enhanced color\n\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "set style data linespoints\n";
  out << "set style fill solid 1.0 noborder\n";
  out << "set ylabel \"Average KL-divergence\"\n";
  out << "plot \"" << avg_kldivs_file << "\" using 1:2 t \"MOM\", \\\n"
      << "\"\" using 1:3 t \"MLE\", \\\n"
      << "\"\" using 1:4 t \"MAP\", \\\n"
      << "\"\" using 1:5 t \"MML\"";
  out.close();
  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

// all e fixed k
void compute_average_kldivs_fixed_kappa(
  string &n_str, string &kappa_str, string &kldivs_folder
) {
  std::vector<Vector> kldivs;
  string kldivs_file;
  string output_file = kldivs_folder + "kldivs_avg.dat";
  ofstream out(output_file.c_str());
  double ecc = 0.1;
  while (ecc < 0.95) {
    out << fixed << setw(10) << setprecision(1) << ecc << "\t\t";
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    string ecc_str = sse.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";
    kldivs_file = current_dir + "kldivs";
    kldivs = load_table(kldivs_file,NUM_METHODS);
    Vector avg_kldivs = computeMeans(kldivs);
    for (int i=0; i<NUM_METHODS; i++) {
      out << fixed << scientific << setprecision(6) << avg_kldivs[i] << "\t\t";
    }
    out << endl;
    ecc += 0.1;
  } // while()
  out.close();
}

void compute_kldivs_diff(
  std::vector<Vector> &kldivs, 
  Vector &kldivs_wins, 
  Vector &kldivs_losses
) {
  double diff;
  kldivs_wins = Vector(3,0);
  kldivs_losses = Vector(3,0);
  for (int i=0; i<kldivs.size(); i++) {
    for (int j=0; j<NUM_METHODS-1; j++) {
      diff = kldivs[i][j] - kldivs[i][MML];
      if (diff < 0) {  // MML loses ...
        kldivs_wins[j] -= diff;
      } else {  // MML wins ...
        kldivs_losses[j] -= diff;
      } // if ()
    } // j ()
  } // i ()
}

double compute_kldivs_diff_fixed_ecc(string &n_str, string &ecc_str, string &kldivs_folder)
{
  double max = 0;
  Vector kldivs_wins,kldivs_losses;
  string output = kldivs_folder + "kldivs_diff.dat";
  ofstream out(output.c_str());
  double kappa = 10;
  while (kappa < 101) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";
    string kldivs_file = current_dir + "kldivs";
    std::vector<Vector> kldivs_table = load_table(kldivs_file,NUM_METHODS);
    compute_kldivs_diff(kldivs_table,kldivs_wins,kldivs_losses);
    out << fixed << setw(10) << setprecision(0) << kappa << "\t\t";
    for (int i=0; i<3; i++) {
      out << fixed << scientific << setprecision(6) << kldivs_wins[i] << "\t\t";
      if (kldivs_wins[i] > max) max = kldivs_wins[i];
    }
    for (int i=0; i<3; i++) {
      out << fixed << scientific << setprecision(6) << kldivs_losses[i] << "\t\t";
      if (fabs(kldivs_losses[i]) > max) max = fabs(kldivs_losses[i]);
    }
    out << endl;
    kappa += 10;
  } // while()
  out.close();
  return max;
}

double compute_kldivs_diff_fixed_kappa(string &n_str, string &kappa_str, string &kldivs_folder)
{
  double max = 0;
  Vector kldivs_wins,kldivs_losses;
  string output = kldivs_folder + "kldivs_diff.dat";
  ofstream out(output.c_str());
  double ecc = 0.1;
  while (ecc < 0.95) {
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    string ecc_str = sse.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";
    string kldivs_file = current_dir + "kldivs";
    std::vector<Vector> kldivs_table = load_table(kldivs_file,NUM_METHODS);
    compute_kldivs_diff(kldivs_table,kldivs_wins,kldivs_losses);
    out << fixed << setw(10) << setprecision(1) << ecc << "\t\t";
    for (int i=0; i<3; i++) {
      out << fixed << scientific << setprecision(6) << kldivs_wins[i] << "\t\t";
      if (kldivs_wins[i] > max) max = kldivs_wins[i];
    }
    for (int i=0; i<3; i++) {
      out << fixed << scientific << setprecision(6) << kldivs_losses[i] << "\t\t";
      if (fabs(kldivs_losses[i]) > max) max = fabs(kldivs_losses[i]);
    }
    out << endl;
    ecc += 0.1;
  } // while()
  out.close();
  return max;
}

void plot_kldivs_diff(string &kldivs_folder, double max)
{
  double ytics;
  if (max >= 1) ytics = 1;
  else if (max > 0.5)  ytics = 0.5;
  else if (max > 0.1) ytics = 0.1;
  else ytics = 0.05;
  string kldivs_diff_file = kldivs_folder + "kldivs_diff.dat";
  string script_file = kldivs_folder + "plot_kldivs_diff.p";
  string plot_file = kldivs_folder + "plot_kldivs_diff.eps";
  ofstream out(script_file.c_str());
  out << "set terminal postscript eps enhanced color\n\n";
  out << "set output \"" << plot_file << "\"\n\n";
  out << "set ylabel \"Difference in KL-divergence\" offset 0,-7,0\n";
  out << "set multiplot layout 2,1\n";
  out << "set lmargin 7\n";
  out << "set rmargin 5\n";
  out << "set ytics " << ytics << "\n";
  out << "unset xtics\n";
  out << "set bmargin 0\n";
  out << "set grid y\n";
  out << "set style data histograms\n";
  out << "set style fill solid 1.0 noborder\n";
  out << "set yr [0:" << max+0.1 << "]\n";
  out << "plot \"" << kldivs_diff_file << "\" using 2 t \"MOM\" lc rgb \"red\", \\\n"
      << "\"\" using 3 t \"MLE\" lc rgb \"green\", \\\n"
      << "\"\" using 4 t \"MAP\" lc rgb \"blue\"\n";
  out << "unset ylabel\n";
  out << "set tmargin 0\n";
  out << "set bmargin at screen 0.1\n";
  out << "set xtics nomirror\n";
  out << "set yr [" << -(max+0.1) << ":0]\n";
  out << "plot \"" << kldivs_diff_file << "\" using 5 notitle lc rgb \"red\", \\\n"
      << "\"\" using 6 notitle lc rgb \"green\", \\\n"
      << "\"\" using 7:xtic(1) notitle lc rgb \"blue\", \\\n"
      << "0 lt 1 lw 10 lc rgb \"black\" notitle\n";
  out << "unset multiplot";
  out.close();
  string cmd = "gnuplot -persist " + script_file;
  if(system(cmd.c_str()));
}

void process_kldivs(string &n_str)
{
  string kldivs_folder;
  string wins_kldivs_file;

  string all_kappas = n_str + "fixed_ecc/";
  string ecc_str,ecc_folder;
  double ecc = 0.1;
  while (ecc < 0.95) {
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    ecc_str = sse.str();
    ecc_folder = all_kappas + "ecc_" + ecc_str + "/";
    kldivs_folder = ecc_folder + "kldivs/";

    // kldivs wins
    wins_kldivs_file = kldivs_folder + "wins_kldivs.dat";
    tabulate_eccentricity_kldivs(wins_kldivs_file,n_str,ecc_str);
    plot_stack_wins_kldivs(kldivs_folder);
    plot_script1_kldivs(kldivs_folder);
    plot_script2_kldivs(kldivs_folder);

    // kldivs
    boxplot_kldivs_fixed_ecc(n_str,ecc_str,kldivs_folder);
    compute_average_kldivs_fixed_ecc(n_str,ecc_str,kldivs_folder);
    plot_avg_kldivs(kldivs_folder);

    // differences
    double max = compute_kldivs_diff_fixed_ecc(n_str,ecc_str,kldivs_folder);
    plot_kldivs_diff(kldivs_folder,max);

    ecc += 0.1;
  } // while(ecc)

  string all_ecc = n_str + "fixed_kappa/";
  string kappa_str,kappa_folder;
  double kappa = 10;
  while (kappa < 101) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    kappa_str = ssk.str();
    kappa_folder = all_ecc + "kappa_" + kappa_str + "/";
    kldivs_folder = kappa_folder + "kldivs/";

    // kldivs wins
    wins_kldivs_file = kldivs_folder + "wins_kldivs.dat";
    tabulate_kappa_kldivs(wins_kldivs_file,n_str,kappa_str);
    plot_stack_wins_kldivs(kldivs_folder);
    plot_script1_kldivs(kldivs_folder);
    plot_script2_kldivs(kldivs_folder);

    // kldivs
    boxplot_kldivs_fixed_kappa(n_str,kappa_str,kldivs_folder);
    compute_average_kldivs_fixed_kappa(n_str,kappa_str,kldivs_folder);
    plot_avg_kldivs(kldivs_folder);

    // differences
    double max = compute_kldivs_diff_fixed_kappa(n_str,kappa_str,kldivs_folder);
    plot_kldivs_diff(kldivs_folder,max);

    kappa += 10;
  } // while(kappa)
}
/* E **********************************************************************************/ 

/* S **********************************************************************************/ 

Vector computeBiasSquared(std::vector<Vector> &p_est, double p)
{
  Vector p_est_means = computeMeans(p_est);
  Vector biassq(NUM_METHODS,0);
  for (int i=0; i<NUM_METHODS; i++) {
    double diff = p_est_means[i] - p;
    biassq[i] = diff * diff;
  }
  return biassq;
}

Vector computeVariance(std::vector<Vector> &p_est, double p)
{
  Vector p_est_means = computeMeans(p_est);
  Vector variance(NUM_METHODS,0);
  for (int i=0; i<p_est.size(); i++) {
    for (int j=0; j<NUM_METHODS; j++) {
      double diff = p_est_means[j] - p_est[i][j];
      variance[j] += (diff * diff);
    }
  }
  for (int j=0; j<NUM_METHODS; j++) {
    variance[j] /= p_est.size();
  }
  return variance;
}

Vector computeMeanSquaredError(std::vector<Vector> &p_est, double p)
{
  Vector p_est_means = computeMeans(p_est);
  Vector mse(NUM_METHODS,0);
  for (int i=0; i<p_est.size(); i++) {
    for (int j=0; j<NUM_METHODS; j++) {
      double diff = p - p_est[i][j];
      mse[j] += (diff * diff);
    }
  }
  for (int j=0; j<NUM_METHODS; j++) {
    mse[j] /= p_est.size();
  }
  return mse;
}

void compute_psi_errors(
  string &errors_folder, string &n_str, string &ecc_str
) {
  std::vector<Vector> psi_est;
  Vector biassq_est,variance_est,mse_est;

  string biassq_file = errors_folder + "biassq_psi";
  string variance_file = errors_folder + "variance_psi";
  string mse_file = errors_folder + "mse_psi";
  ofstream biassq(biassq_file.c_str());
  ofstream variance(variance_file.c_str());
  ofstream mse(mse_file.c_str());
  double kappa = 10;
  while (kappa < 101) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";

    string psi_file = current_dir + "psi_est";
    psi_est = load_table(psi_file,NUM_METHODS);
    biassq_est = computeBiasSquared(psi_est,PSI);
    variance_est = computeVariance(psi_est,PSI);
    mse_est = computeMeanSquaredError(psi_est,PSI);

    biassq << fixed << setw(10) << setprecision(0) << kappa << "\t";
    variance << fixed << setw(10) << setprecision(0) << kappa << "\t";
    mse << fixed << setw(10) << setprecision(0) << kappa << "\t";
    for(int i=0; i<NUM_METHODS; i++) {
      biassq << scientific << setprecision(6) << biassq_est[i] << "\t";
      variance << scientific << setprecision(6) << variance_est[i] << "\t";
      mse << scientific << setprecision(6) << mse_est[i] << "\t";
    }
    biassq << endl;
    variance << endl;
    mse << endl;

    kappa += 10;
  } // while()
  biassq.close();
  variance.close();
  mse.close();
}

void compute_alpha_errors(
  string &errors_folder, string &n_str, string &ecc_str
) {
  std::vector<Vector> alpha_est;
  Vector biassq_est,variance_est,mse_est;

  string biassq_file = errors_folder + "biassq_alpha";
  string variance_file = errors_folder + "variance_alpha";
  string mse_file = errors_folder + "mse_alpha";
  ofstream biassq(biassq_file.c_str());
  ofstream variance(variance_file.c_str());
  ofstream mse(mse_file.c_str());
  double kappa = 10;
  while (kappa < 101) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";

    string alpha_file = current_dir + "alpha_est";
    alpha_est = load_table(alpha_file,NUM_METHODS);
    biassq_est = computeBiasSquared(alpha_est,ALPHA);
    variance_est = computeVariance(alpha_est,ALPHA);
    mse_est = computeMeanSquaredError(alpha_est,ALPHA);

    biassq << fixed << setw(10) << setprecision(0) << kappa << "\t";
    variance << fixed << setw(10) << setprecision(0) << kappa << "\t";
    mse << fixed << setw(10) << setprecision(0) << kappa << "\t";
    for(int i=0; i<NUM_METHODS; i++) {
      biassq << scientific << setprecision(6) << biassq_est[i] << "\t";
      variance << scientific << setprecision(6) << variance_est[i] << "\t";
      mse << scientific << setprecision(6) << mse_est[i] << "\t";
    }
    biassq << endl;
    variance << endl;
    mse << endl;

    kappa += 10;
  } // while()
  biassq.close();
  variance.close();
  mse.close();
}

void compute_eta_errors(
  string &errors_folder, string &n_str, string &ecc_str
) {
  std::vector<Vector> eta_est;
  Vector biassq_est,variance_est,mse_est;

  string biassq_file = errors_folder + "biassq_eta";
  string variance_file = errors_folder + "variance_eta";
  string mse_file = errors_folder + "mse_eta";
  ofstream biassq(biassq_file.c_str());
  ofstream variance(variance_file.c_str());
  ofstream mse(mse_file.c_str());
  double kappa = 10;
  while (kappa < 101) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";

    string eta_file = current_dir + "eta_est";
    eta_est = load_table(eta_file,NUM_METHODS);
    biassq_est = computeBiasSquared(eta_est,ETA);
    variance_est = computeVariance(eta_est,ETA);
    mse_est = computeMeanSquaredError(eta_est,ETA);

    biassq << fixed << setw(10) << setprecision(0) << kappa << "\t";
    variance << fixed << setw(10) << setprecision(0) << kappa << "\t";
    mse << fixed << setw(10) << setprecision(0) << kappa << "\t";
    for(int i=0; i<NUM_METHODS; i++) {
      biassq << scientific << setprecision(6) << biassq_est[i] << "\t";
      variance << scientific << setprecision(6) << variance_est[i] << "\t";
      mse << scientific << setprecision(6) << mse_est[i] << "\t";
    }
    biassq << endl;
    variance << endl;
    mse << endl;

    kappa += 10;
  } // while()
  biassq.close();
  variance.close();
  mse.close();
}

void compute_kappa_errors(
  string &errors_folder, string &n_str, string &ecc_str
) {
  std::vector<Vector> kappas_est;
  Vector biassq_est,variance_est,mse_est;

  string biassq_file = errors_folder + "biassq_kappa";
  string variance_file = errors_folder + "variance_kappa";
  string mse_file = errors_folder + "mse_kappa";
  ofstream biassq(biassq_file.c_str());
  ofstream variance(variance_file.c_str());
  ofstream mse(mse_file.c_str());
  double kappa = 10;
  while (kappa < 101) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";

    string kappas_file = current_dir + "kappa_est";
    kappas_est = load_table(kappas_file,NUM_METHODS);
    biassq_est = computeBiasSquared(kappas_est,kappa);
    variance_est = computeVariance(kappas_est,kappa);
    mse_est = computeMeanSquaredError(kappas_est,kappa);

    biassq << fixed << setw(10) << setprecision(0) << kappa << "\t";
    variance << fixed << setw(10) << setprecision(0) << kappa << "\t";
    mse << fixed << setw(10) << setprecision(0) << kappa << "\t";
    for(int i=0; i<NUM_METHODS; i++) {
      biassq << scientific << setprecision(6) << biassq_est[i] << "\t";
      variance << scientific << setprecision(6) << variance_est[i] << "\t";
      mse << scientific << setprecision(6) << mse_est[i] << "\t";
    }
    biassq << endl;
    variance << endl;
    mse << endl;

    kappa += 10;
  } // while()
  biassq.close();
  variance.close();
  mse.close();
}

void compute_beta_errors(
  string &errors_folder, string &n_str, double ecc, string &ecc_str
) {
  std::vector<Vector> betas_est;
  Vector biassq_est,variance_est,mse_est;

  string biassq_file = errors_folder + "biassq_beta";
  string variance_file = errors_folder + "variance_beta";
  string mse_file = errors_folder + "mse_beta";
  ofstream biassq(biassq_file.c_str());
  ofstream variance(variance_file.c_str());
  ofstream mse(mse_file.c_str());
  double kappa = 10;
  while (kappa < 101) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    string kappa_str = ssk.str();
    string current_dir = n_str + "k_" + kappa_str + "_e_" + ecc_str + "/";

    double beta = 0.5 * ecc * kappa;

    string betas_file = current_dir + "beta_est";
    betas_est = load_table(betas_file,NUM_METHODS);
    biassq_est = computeBiasSquared(betas_est,beta);
    variance_est = computeVariance(betas_est,beta);
    mse_est = computeMeanSquaredError(betas_est,beta);

    biassq << fixed << setw(10) << setprecision(0) << kappa << "\t";
    variance << fixed << setw(10) << setprecision(0) << kappa << "\t";
    mse << fixed << setw(10) << setprecision(0) << kappa << "\t";
    for(int i=0; i<NUM_METHODS; i++) {
      biassq << scientific << setprecision(6) << biassq_est[i] << "\t";
      variance << scientific << setprecision(6) << variance_est[i] << "\t";
      mse << scientific << setprecision(6) << mse_est[i] << "\t";
    }
    biassq << endl;
    variance << endl;
    mse << endl;

    kappa += 10;
  } // while()
  biassq.close();
  variance.close();
  mse.close();
}

void combine(std::vector<string> &files, string &output_file)
{
  std::vector<std::vector<Vector> > all_tables;

  for (int i=0; i<files.size(); i++) {
    std::vector<Vector> table = load_table(files[i],NUM_METHODS+1);
    all_tables.push_back(table);
  }

  int num_kappas = all_tables[0].size();
  ofstream out(output_file.c_str());

  for (int i=0; i<num_kappas; i++) {
    Vector sum(NUM_METHODS,0);
    for (int j=0; j<all_tables.size(); j++) {
      for (int k=0; k<NUM_METHODS; k++) {
        sum[k] += all_tables[j][i][k+1];
      } // k
    } // j
    out << fixed << setw(10) << setprecision(0) << all_tables[0][i][0] << "\t\t";
    for (int k=0; k<NUM_METHODS; k++) {
      out << scientific << setprecision(6) << sum[k] << "\t\t";
    } // k
    out << endl;
  } // i

  out.close();
}

void compute_all_errors(
  string &errors_folder, string &n_str, string &ecc_str
) {
  std::vector<string> biassq_files,variance_files,mse_files;
  string biassq_file,variance_file,mse_file,output_file;

  biassq_file = errors_folder + "biassq_psi";
  biassq_files.push_back(biassq_file);
  biassq_file = errors_folder + "biassq_alpha";
  biassq_files.push_back(biassq_file);
  biassq_file = errors_folder + "biassq_eta";
  biassq_files.push_back(biassq_file);
  biassq_file = errors_folder + "biassq_kappa";
  biassq_files.push_back(biassq_file);
  biassq_file = errors_folder + "biassq_beta";
  biassq_files.push_back(biassq_file);
  output_file = errors_folder + "biassq_all";
  combine(biassq_files,output_file);

  variance_file = errors_folder + "variance_psi";
  variance_files.push_back(variance_file);
  variance_file = errors_folder + "variance_alpha";
  variance_files.push_back(variance_file);
  variance_file = errors_folder + "variance_eta";
  variance_files.push_back(variance_file);
  variance_file = errors_folder + "variance_kappa";
  variance_files.push_back(variance_file);
  variance_file = errors_folder + "variance_beta";
  variance_files.push_back(variance_file);
  output_file = errors_folder + "variance_all";
  combine(variance_files,output_file);

  mse_file = errors_folder + "mse_psi";
  mse_files.push_back(mse_file);
  mse_file = errors_folder + "mse_alpha";
  mse_files.push_back(mse_file);
  mse_file = errors_folder + "mse_eta";
  mse_files.push_back(mse_file);
  mse_file = errors_folder + "mse_kappa";
  mse_files.push_back(mse_file);
  mse_file = errors_folder + "mse_beta";
  mse_files.push_back(mse_file);
  output_file = errors_folder + "mse_all";
  combine(mse_files,output_file);
}

void process_estimates(string &n_str)
{
  string errors_folder;
  string errors_file;

  string all_kappas = n_str + "fixed_ecc/";
  string ecc_str,ecc_folder;
  double ecc = 0.1;
  while (ecc < 0.95) {
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    ecc_str = sse.str();
    ecc_folder = all_kappas + "ecc_" + ecc_str + "/";
    errors_folder = ecc_folder + "errors/";

    // psi errors
    compute_psi_errors(errors_folder,n_str,ecc_str);
    // alpha errors
    compute_alpha_errors(errors_folder,n_str,ecc_str);
    // eta errors
    compute_eta_errors(errors_folder,n_str,ecc_str);

    // kappa errors
    compute_kappa_errors(errors_folder,n_str,ecc_str);
    // beta errors
    compute_beta_errors(errors_folder,n_str,ecc,ecc_str);

    // combined errors
    compute_all_errors(errors_folder,n_str,ecc_str);

    ecc += 0.1;
  } // while(ecc)

}

/* E **********************************************************************************/ 

/* S **********************************************************************************/ 
void check_and_create_directory(string &directory)
{
  if (stat(directory.c_str(), &st) == -1) {
    mkdir(directory.c_str(), 0700);
  }
}

void create_required_folders(string &n_str)
{
  string kldivs_folder,errors_folder;

  string all_kappas = n_str + "fixed_ecc/";
  check_and_create_directory(all_kappas);
  string ecc_str,ecc_folder;
  double ecc = 0.1;
  while (ecc < 0.95) {
    ostringstream sse;
    sse << fixed << setprecision(1);
    sse << ecc;
    ecc_str = sse.str();
    ecc_folder = all_kappas + "ecc_" + ecc_str + "/";
    check_and_create_directory(ecc_folder);

    kldivs_folder = ecc_folder + "kldivs/";
    check_and_create_directory(kldivs_folder);

    errors_folder = ecc_folder + "errors/";
    check_and_create_directory(errors_folder);

    ecc += 0.1;
  } // while(ecc)

  string all_ecc = n_str + "fixed_kappa/";
  check_and_create_directory(all_ecc);
  string kappa_str,kappa_folder;
  double kappa = 10;
  while (kappa < 101) {
    ostringstream ssk;
    ssk << fixed << setprecision(0);
    ssk << kappa;
    kappa_str = ssk.str();
    kappa_folder = all_ecc + "kappa_" + kappa_str + "/";
    check_and_create_directory(kappa_folder);

    kldivs_folder = kappa_folder + "kldivs/";
    check_and_create_directory(kldivs_folder);

    errors_folder = kappa_folder + "errors/";
    check_and_create_directory(errors_folder);

    kappa += 10;
  } // while(kappa)
}
/* E **********************************************************************************/

struct Parameters
{
  int N;
};

struct Parameters parseCommandLineInput(int argc, char **argv)
{
  struct Parameters parameters;

  options_description desc("Allowed options");
  desc.add_options()
       ("n",value<int>(&parameters.N),"sample size")
  ;
  variables_map vm;
  store(command_line_parser(argc,argv).options(desc).run(),vm);
  notify(vm);

  return parameters;
}

int main(int argc, char **argv)
{
  struct Parameters parameters = parseCommandLineInput(argc,argv);

  //string n_str = "./N_" + boost::lexical_cast<string>(parameters.N) + "_prior1/";
  string n_str = "./N_" + boost::lexical_cast<string>(parameters.N) + "_prior2/";
  create_required_folders(n_str);

  process_kldivs(n_str);

  process_estimates(n_str);

}

