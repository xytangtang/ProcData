#include <Rcpp.h>
using namespace Rcpp;

double sum_abs_diff(NumericVector x, NumericVector y) {
  int l = std::min(x.length(), y.length());
  double res = 0;
  for(int i=0; i<l; i++){
    int diff = x[i] - y[i];
    res += diff > 0? diff : -diff;
  }
  return res;
}

int sum_abs_diff_int(IntegerVector x, IntegerVector y) {
  int l = std::min(x.length(), y.length());
  int res = 0;
  for(int i=0; i<l; i++){
    int diff = x[i] - y[i];
    res += diff > 0? diff : -diff;
  }
  return res;
}

std::map<String, IntegerVector> count_events(StringVector seq) {
  int l = seq.length();
  std::map<String, IntegerVector> dict;
  for(int i=0; i<l; i++){
    String event = seq[i];
    dict[event].push_back(i);
  }
  return dict;
}

double calculate_f1_cpp(StringVector seq1, StringVector seq2) {
  StringVector common_events = intersect(seq1, seq2);
  int n_common_events = common_events.length();
  double f_score = 0;
  std::map<String, IntegerVector> dict1 = count_events(seq1);
  std::map<String, IntegerVector> dict2 = count_events(seq2);
  for(int i=0; i<n_common_events; i++){
    String event = common_events[i];
    int l = std::min(dict1[event].length(), dict2[event].length());
    f_score += sum_abs_diff_int(dict1[event], dict2[event]);
  }
  return f_score/std::max(seq1.length(), seq2.length());
}

int calculate_g_cpp(StringVector seq1, StringVector seq2) {
  return setdiff(seq1, seq2).length() + setdiff(seq2, seq1).length();
}


double calculate_time_f1_cpp(StringVector seq1, StringVector seq2, NumericVector ts1, NumericVector ts2) {
  StringVector common_events = intersect(seq1, seq2);
  int n_common_events = common_events.length();
  double f_score = 0;
  std::map<String, IntegerVector> dict1 = count_events(seq1);
  std::map<String, IntegerVector> dict2 = count_events(seq2);
  for(int i=0; i<n_common_events; i++){
    String event = common_events[i];
    int l = std::min(dict1[event].length(), dict2[event].length());
    f_score += sum_abs_diff(ts1[dict1[event]], ts2[dict2[event]]);
  }
  return f_score/std::max(seq1.length(), seq2.length());
}


// [[Rcpp::export]]
double calculate_dissimilarity_cpp(StringVector seq1, StringVector seq2) {
  int l1 = seq1.length();
  int l2 = seq2.length();
  double f_score = calculate_f1_cpp(seq1, seq2);
  double g_score = calculate_g_cpp(seq1, seq2);
  return (f_score+g_score)/(l1+l2);
}

// [[Rcpp::export]]
double group_score_cpp(List list1, List list2) {
  double ave = 0;
  int l1 = list1.length();
  int l2 = list2.length();
  for(int i=0; i<l1; i++)
    for(int j=0; j<l2; j++)
      ave += calculate_dissimilarity_cpp(list1[i], list2[j]);
  return ave/(l1*l2);
}

// [[Rcpp::export]]
NumericMatrix calculate_dist_cpp(List seqs) {
  int n = seqs.length();
  NumericMatrix Dmat(n,n);
  for(int i=1; i<n; i++){
    for(int j=0; j<i; j++){
      Dmat(i,j) = calculate_dissimilarity_cpp(seqs[i],seqs[j]);
      Dmat(j,i) = Dmat(i,j);
    }
  }
  return Dmat;
}

// [[Rcpp::export]]
NumericMatrix calculate_group_dist_cpp(List seqs) {
  int n = seqs.length();
  NumericMatrix Dmat(n,n);
  for(int i=1; i<n; i++){
    for(int j=0; j<i; j++){
      Dmat(i,j) = group_score_cpp(seqs[i],seqs[j]);
      Dmat(j,i) = Dmat(i,j);
    }
  }
  return Dmat;
}