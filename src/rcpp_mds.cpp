#include <Rcpp.h>
using namespace Rcpp;

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

Rcpp::IntegerVector randomShuffle(Rcpp::IntegerVector a) {
    
    // clone a into b to leave a alone
    Rcpp::IntegerVector b = Rcpp::clone(a);
    
    std::random_shuffle(b.begin(), b.end(), randWrapper);
    
    return b;
}

double myabs(double x) {
    return x > 0 ? x : -x;
}

void center_mds (NumericMatrix Theta)
{
    int N = Theta.nrow();
    int K = Theta.ncol();
    
    double theta_mean = 0.0;
    
    for (int k = 0; k < K; k++)
    {
        theta_mean = 0.0;
        for (int i = 0; i < N; i++)
        {
            theta_mean += Theta(i, k);
        }
        
        theta_mean /= N;
        
        for (int i = 0; i < N; i++)
        {
            Theta(i, k) -= theta_mean;
        }
    }
}

void normalize_mmds(NumericMatrix Theta, NumericMatrix W)
{
    int N = Theta.nrow();
    int K = Theta.ncol();
    int L = W.nrow();
    
    double norm = 0.0;
    double sum_theta = 0.0;
    
    for (int k = 0; k < K; k++)
    {
        sum_theta = 0.0;
        for (int i = 0; i < N; i++)
        {
            sum_theta += Theta(i, k);
        }
        sum_theta /= N;
        
        for (int i = 0; i < N; i++)
        {
            Theta(i, k) -= sum_theta;
        }
    }
    
    for (int k = 0; k < K; k++)
    {
        norm = 0.0;
        
        for (int i = 0; i < N; i++)
        {
            norm += Theta(i, k) * Theta(i, k);
        }
        
        for (int i = 0; i < N; i++)
        {
            Theta(i, k) = Theta(i, k) / sqrt(norm / N);
        }
        
        for (int l = 0; l < L; l++)
        {
            W(l, k) = W(l, k) * norm / N;
        }
    }
}

double calculate_dist_l2(NumericMatrix Theta, int i, int j) {
    
    double dist = 0.0;
    int K = Theta.ncol();
    
    for (int k = 0; k < K; k++) {
        dist += (Theta(i, k) - Theta(j, k)) * (Theta(i, k) - Theta(j, k));
    }
    
    dist = sqrt(dist);
    
    return dist;
    
}

double calculate_loss(NumericMatrix D, NumericMatrix Theta) {
    
    double delta = 0.0, r = 0.0;
    int i, j;
    int N = D.nrow();
    
    for (i = 0; i < N - 1; i++) {
        for (j = i+1; j < N; j++) {
            r = calculate_dist_l2(Theta, i, j);
            delta += (r - D(i, j)) * (r - D(i, j));
        }
    }
    
    return delta/N/(N-1)/2.0;
}

double calculate_loss_subset(NumericMatrix D, NumericMatrix Theta, IntegerMatrix sample_set) {
    
    double delta = 0.0, r = 0.0;
    int i, j, index_sample;
    int n_sample = sample_set.nrow();
    
    for (index_sample = 0; index_sample < n_sample; index_sample++)
    {
        i = sample_set(index_sample, 0);
        j = sample_set(index_sample, 1);
        r = calculate_dist_l2(Theta, i, j);
        delta += (r - D(i, j)) * (r - D(i, j));
    }
    
    return delta/n_sample;
}

// [[Rcpp::export]]
List MDS(NumericMatrix D, NumericMatrix Theta, int n_epoch, double step_size, double tot, unsigned int seed) {
    int N = Theta.nrow();
    int K = Theta.ncol();
    int i = 0, j = 0, k = 0;
    int convergence_tag = 0;
    double d_data = 0.0, d_est = 0.0, part1 = 0.0, part2 = 0.0, part3 = 0.0;
    
    IntegerVector sample_order(N);
    double loss_old = 0.0, loss_new = 0.0;
    
    for (i = 0; i < N; i++) sample_order[i] = i;
    
    for (int index_epoch = 0; index_epoch < n_epoch; index_epoch++)
    {
        loss_old = loss_new;
        
        sample_order = randomShuffle(sample_order);
        
        for (int index_sample = 0; index_sample < N; index_sample++)
        {
            i = sample_order[index_sample];
            
            j = floor(runif(1)[0] * (N - 1));
            if (j >= i) j += 1;
            
            d_data = D(i, j);
            d_est = calculate_dist_l2(Theta, i, j);
            
            part1 = 1.0 - d_data / d_est;
            
            for (k = 0; k < K; k++)
            {
                part2 = Theta(i, k) - Theta(j, k);
                part3 = 2.0 * part1 * part2 * step_size;
                Theta(i, k) -= part3;
                Theta(j, k) += part3;
            }
            
        }
        
        center_mds(Theta);
        
        loss_new = calculate_loss(D, Theta);
        
        if (myabs(loss_new - loss_old) < tot) {
            convergence_tag = 1;
            break;
            
        }
    }
    
    return List::create(Named("convergence") = convergence_tag,
                        Named("loss") = loss_new);

}

// [[Rcpp::export]]
List MDS_subset(NumericMatrix D, NumericMatrix Theta, int n_epoch, double step_size, double tot, IntegerMatrix train_set, IntegerMatrix valid_set) {
    int K = Theta.ncol();
    int i = 0, j = 0, k = 0, index_epoch = 0;
    int convergence_tag = 0;
    double d_data = 0.0, d_est = 0.0, part1 = 0.0, part2 = 0.0, part3 = 0.0;
    
    int n_train = train_set.nrow();
    IntegerVector train_idx = Range(0, n_train-1);
    
    double loss_old = 0.0, loss_new = 0.0, valid_loss = 0.0;
    NumericVector loss_trace(n_epoch);
    for (index_epoch = 0; index_epoch < n_epoch; index_epoch++)
    {
        loss_old = loss_new;
        
        train_idx = randomShuffle(train_idx);
        
        for (int index_sample = 0; index_sample < n_train; index_sample++)
        {
            i = train_set(train_idx[index_sample], 0);
            j = train_set(train_idx[index_sample], 1);
            
            d_data = D(i, j);
            d_est = calculate_dist_l2(Theta, i, j);
            
            part1 = 1.0 - d_data / d_est;
            
            for (k = 0; k < K; k++)
            {
                part2 = Theta(i, k) - Theta(j, k);
                part3 = 2.0 * part1 * part2 * step_size;
                Theta(i, k) -= part3;
                Theta(j, k) += part3;
            }
            
        }
        
        center_mds(Theta);
        
        loss_new = calculate_loss_subset(D, Theta, train_set);
        loss_trace[index_epoch] = loss_new;
        if (myabs(loss_new - loss_old) < tot) {
            convergence_tag = 1;
            break;
        }
    }
    
    valid_loss = calculate_loss_subset(D, Theta, valid_set);
    
    return List::create(Named("convergence") = convergence_tag,
                        Named("train_loss") = loss_new,
                        Named("valid_loss") = valid_loss);
}

// place center algorithm for mds
//double calculate_dist1(NumericMatrix Theta, int i, int j) {
//    
//    double dist = 0.0;
//    int K = Theta.ncol();
//    
//    for (int k = 0; k < K; k++) {
//        dist += (Theta(i, k) - Theta(j, k)) * (Theta(i, k) - Theta(j, k));
//    }
//    
//    dist = sqrt(dist);
//    
//    return dist;
//    
//}
//
//double calculate_dist2(NumericMatrix Theta1, int i, NumericMatrix Theta2, int j) {
//    
//    double dist = 0.0;
//    int K = Theta1.ncol();
//    
//    for (int k = 0; k < K; k++) {
//        dist += (Theta1(i, k) - Theta2(j, k)) * (Theta1(i, k) - Theta2(j, k));
//    }
//    
//    dist = sqrt(dist);
//    
//    return dist;
//}
//
//
//double calculate_loss1(NumericMatrix D, NumericMatrix Theta) {
//    
//    double delta = 0.0, r = 0.0;
//    int i, j;
//    int N = D.nrow();
//    
//    for (i = 0; i < N - 1; i++) {
//        for (j = i+1; j < N; j++) {
//            r = calculate_dist1(Theta, i, j);
//            delta += (r - D(i, j)) * (r - D(i, j));
//        }
//    }
//    
//    return delta;
//}
//
//double calculate_loss2(NumericMatrix D_new2old, NumericMatrix D_new2new, NumericMatrix Theta_old, NumericMatrix Theta_new) {
//    
//    double delta = 0.0, r = 0.0;
//    int i, j;
//    int N_new = D_new2old.nrow(), N_old = D_new2old.ncol();
//    
//    for (i = 0; i < N_new - 1; i++) {
//        for (j = i+1; j < N_new; j++) {
//            r = calculate_dist1(Theta_new, i, j);
//            delta += (r - D_new2new(i, j)) * (r - D_new2new(i, j));
//        }
//    }
//    
//    for (i = 0; i < N_new; i++) {
//        for (j = 0; j < N_old; j++) {
//            r = calculate_dist2(Theta_new, i, Theta_old, j);
//            delta += (r - D_new2old(i, j)) * (r - D_new2old(i, j));
//        }
//    }
//    
//    return delta;
//}
//
//
//
//void fMDS(NumericMatrix D, NumericMatrix Theta, double tot) {
//    int N = D.nrow();
//    int K = Theta.ncol();
//    int i, j, k;
//    
//    NumericMatrix xhat(N, K);
//    
//    double delta = 1.0, loss_old = 0.0, loss_new = 0.0, alpha;
//    loss_new = calculate_loss1(D, Theta);
//    while (delta > tot) {
//        loss_old = loss_new;
//        for (i = 0; i < N; i++)
//        {
//            // calculate xhat
//            for (j = 0; j < N; j++)
//            {
//                if (j != i) {
//                    alpha = 1.0 - D(i,j) / calculate_dist1(Theta, i, j);
//                    for (k = 0; k < K; k++) xhat(j, k) = Theta(j, k) * alpha + (1.0 - alpha) * Theta(i, k);
//                }
//            }
//            
//            // update x
//            for (k = 0; k < K; k++) {
//                Theta(i, k) = 0.0;
//                for (j = 0; j < N; j++) {
//                    if (j != i) Theta(i, k) += xhat(j, k);
//                }
//                Theta(i, k) = Theta(i, k) / (N - 1);
//            }
//        }
//        loss_new = calculate_loss1(D, Theta);
//        delta = myabs(loss_old - loss_new);
//    }
//}
