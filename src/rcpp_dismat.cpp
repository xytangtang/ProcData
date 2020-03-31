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

double sum_matched(StringVector x, StringVector y, NumericVector xx){
  double res = 0;
  int l = x.length();
  LogicalVector a = is_na(match(x,y));
  for(int i=0; i<l; i++) res+= a[i]? xx[i] : 0;
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
    f_score += sum_abs_diff(ts1[dict1[event]], ts2[dict2[event]]);
  }
  return f_score/std::max(ts1[seq1.length()-1], ts2[seq2.length()-1]);
}

double calculate_time_g_cpp(StringVector seq1, StringVector seq2, NumericVector ts1, NumericVector ts2) {
  double res = 0;
  NumericVector ti1 = diff(ts1), ti2 = diff(ts2);
  ti1.push_front(0);
  ti2.push_front(0);
  res += sum_matched(seq1, seq2, ti1) + sum_matched(seq2, seq1, ti2);
  return res;
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
double calculate_tdissimilarity_cpp(StringVector seq1, StringVector seq2, NumericVector ts1, NumericVector ts2) {
  int l1 = seq1.length();
  int l2 = seq2.length();
  double f_score = calculate_time_f1_cpp(seq1, seq2, ts1, ts2);
  double g_score = calculate_time_g_cpp(seq1, seq2, ts1, ts2);
  return (f_score+g_score)/(ts1[l1-1]+ts2[l2-1]);
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

//' calculate dissimilarity matrix
//' @param seqs a list of action sequences
//' @export
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
NumericMatrix calculate_tdist_cpp(List action_seqs, List time_seqs) {
  int n = action_seqs.length();
  NumericMatrix Dmat(n,n);
  for(int i=1; i<n; i++){
    for(int j=0; j<i; j++){
      Dmat(i,j) = calculate_tdissimilarity_cpp(action_seqs[i],action_seqs[j], time_seqs[i], time_seqs[j]);
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

// [[Rcpp::export]]
double calculate_common_event_proportion(StringVector seq1, StringVector seq2) {
  StringVector common_events = intersect(seq1, seq2);
  int n_common_events = common_events.size();
  StringVector all_events = union_(seq1, seq2);
  int n_all_events = all_events.size();
  
  return (0.0 + n_common_events) / (n_all_events + 0.0);
}

// [[Rcpp::export]]
StringVector paste_seq(StringVector seq, int k) {
  int L = seq.size();
  StringVector pasted_seq(L - k + 1);
  for (int i = 0; i < L - k + 1; i++) {
    String new_event;
    for (int j = 0; j < k; j++) {
      new_event += seq[i+j];
      if (j != (k-1)) {
        new_event += "\t";
      }
    }
    pasted_seq[i] = new_event;
  }
  
  return pasted_seq;
}

// [[Rcpp::export]]
double calculate_common_ngram_proportion(StringVector seq1, StringVector seq2, int k) {
  StringVector seq1_pasted = paste_seq(seq1, k);
  StringVector seq2_pasted = paste_seq(seq2, k);
  return calculate_common_event_proportion(seq1_pasted, seq2_pasted);
}

// [[Rcpp::export]]
double calculate_ngram_dissimilarity(StringVector seq1, StringVector seq2) {
  
  int L = 3;
  NumericVector p_vec(L);
  for (int l=0; l < L; l++) {
    p_vec[l] = calculate_common_ngram_proportion(seq1, seq2, l+1);
  }
  
  return exp(mean(log(1.0 - p_vec)));
  
}

// [[Rcpp::export]]
NumericMatrix calculate_ngram_dist_cpp(List seqs) {
  int n = seqs.length();
  NumericMatrix Dmat(n,n);
  for(int i=1; i<n; i++){
    for(int j=0; j<i; j++){
      Dmat(i,j) = calculate_ngram_dissimilarity(seqs[i],seqs[j]);
      Dmat(j,i) = Dmat(i,j);
    }
  }
  return Dmat;
}