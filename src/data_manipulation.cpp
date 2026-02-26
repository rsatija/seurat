#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <Rinternals.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]



// [[Rcpp::export]]
Eigen::SparseMatrix<double> RunUMISampling(Eigen::SparseMatrix<double> data, int sample_val, bool upsample = false, bool display_progress=true){
    Progress p(data.outerSize(), display_progress);
    Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
    for (int k=0; k < data.outerSize(); ++k){
      p.increment();
      for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
        double entry = it.value();
        if( (upsample) || (colSums[k] > sample_val)){
          entry = entry * double(sample_val) / colSums[k];
          if (fmod(entry, 1) != 0){
            double rn = R::runif(0,1);
            if(fmod(entry, 1) <= rn){
              it.valueRef() = floor(entry);
            }
            else{
              it.valueRef() = ceil(entry);
            }
          }
          else{
            it.valueRef() = entry;
          }
        }
      }
    }
  return(data);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> RunUMISamplingPerCell(Eigen::SparseMatrix<double> data, NumericVector sample_val, bool upsample = false, bool display_progress=true){
  Progress p(data.outerSize(), display_progress);
  Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      double entry = it.value();
      if( (upsample) || (colSums[k] > sample_val[k])){
        entry = entry * double(sample_val[k]) / colSums[k];
        if (fmod(entry, 1) != 0){
          double rn = R::runif(0,1);
          if(fmod(entry, 1) <= rn){
            it.valueRef() = floor(entry);
          }
          else{
            it.valueRef() = ceil(entry);
          }
        }
        else{
          it.valueRef() = entry;
        }
      }
    }
  }
  return(data);
}


typedef Eigen::Triplet<double> T;
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> RowMergeMatrices(Eigen::SparseMatrix<double, Eigen::RowMajor> mat1, Eigen::SparseMatrix<double, Eigen::RowMajor> mat2, std::vector< std::string > mat1_rownames,
                                             std::vector< std::string > mat2_rownames, std::vector< std::string > all_rownames){


  // Set up hash maps for rowname based lookup
  std::unordered_map<std::string, int> mat1_map;
  for(unsigned int i = 0; i < mat1_rownames.size(); i++){
    mat1_map[mat1_rownames[i]] = i;
  }
  std::unordered_map<std::string, int> mat2_map;
  for(unsigned int i = 0; i < mat2_rownames.size(); i++){
    mat2_map[mat2_rownames[i]] = i;
  }

  // set up tripletList for new matrix creation
  std::vector<T> tripletList;
  int num_rows = all_rownames.size();
  int num_col1 = mat1.cols();
  int num_col2 = mat2.cols();


  tripletList.reserve(mat1.nonZeros() + mat2.nonZeros());
  for(int i = 0; i < num_rows; i++){
    std::string key = all_rownames[i];
    if (mat1_map.count(key)){
      for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it1(mat1, mat1_map[key]); it1; ++it1){
        tripletList.emplace_back(i, it1.col(), it1.value());
      }
    }
    if (mat2_map.count(key)){
      for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it2(mat2, mat2_map[key]); it2; ++it2){
        tripletList.emplace_back(i, num_col1 + it2.col(), it2.value());
      }
    }
  }
  Eigen::SparseMatrix<double> combined_mat(num_rows, num_col1 + num_col2);
  combined_mat.setFromTriplets(tripletList.begin(), tripletList.end());
  return combined_mat;
}

// Log-normalize sparse matrix columns for Seurat's NormalizeData.fast path.
//
// For each column k:
//   data[k, i] <- log1p(data[k, i] / colSums[k] * scale_factor)
//
// Notes:
// - Input is passed by value, so caller receives a new object.
// - Uses RcppProgress progress bar when display_progress is TRUE.
// - Assumes column sums are strictly positive for sparse normalized input.
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> LogNorm(Eigen::SparseMatrix<double> data, int scale_factor, bool display_progress = true){
  Progress p(data.outerSize(), display_progress);
  const int ncols = data.cols();
  const double sf = static_cast<double>(scale_factor);
  const int *outer = data.outerIndexPtr();
  double *x = data.valuePtr();
  bool has_nonfinite = false;
  std::vector<double> col_sums(ncols, 0.0);
  for (int k = 0; k < ncols; ++k) {
    const int start = outer[k];
    const int end = outer[k + 1];
    double col_sum = 0.0;
    for (int j = start; j < end; ++j) {
      if (!std::isfinite(x[j])) {
        has_nonfinite = true;
      }
      col_sum += x[j];
    }
    col_sums[k] = col_sum;
  }
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    const int start = outer[k];
    const int end = outer[k + 1];
    const double scale = sf / col_sums[k];
    for (int j = start; j < end; ++j) {
      if (has_nonfinite) {
        if (!std::isfinite(x[j])) {
          continue;
        }
      }
      x[j] = log1p(x[j] * scale);
    }
  }
  return data;
}

//' Calculate sparse row means and unbiased variances in a single pass.
//'
//' This provides a fused C++ implementation for the VST branch:
//'   - mean[x] = mean value across cells for each feature
//'   - variance[x] = unbiased variance across cells for each feature
//'
//' Equivalent to the legacy pair of rowMeans + SparseRowVar2 for nonzero-aware
//' sparse inputs when input has no nonfinite values.
// [[Rcpp::export(rng = false)]]
List SparseRowStats(Eigen::SparseMatrix<double> mat, bool display_progress = true){
  mat = mat.transpose();
  const int ncols = mat.cols();
  const int nrows = mat.rows();
  NumericVector means(ncols);
  NumericVector variances(ncols);
  if (display_progress == true) {
    Rcpp::Rcerr << "Calculating gene means and variances" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k = 0; k < mat.outerSize(); ++k) {
    p.increment();
    double col_sum = 0;
    double col_sq_sum = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
      col_sum += it.value();
      col_sq_sum += it.value() * it.value();
    }
    const double feature_mean = col_sum / nrows;
    means[k] = feature_mean;
    variances[k] = (col_sq_sum - (nrows * feature_mean * feature_mean)) / (nrows - 1);
  }
  return List::create(
    _["mean"] = means,
    _["variance"] = variances
  );
}

// CLR-normalize sparse matrix by margin.
// Margin 1 = rows, Margin 2 = columns.
// For each vector v:
//   v <- log1p(v / exp(sum(log1p(v[v > 0], na.rm = TRUE) / length(v)))
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> CLRNorm(Eigen::SparseMatrix<double> data, int margin, bool display_progress = true){
  if ((margin != 1) && (margin != 2)) {
    stop("`margin` must be 1 or 2");
  }
  const int nrows = data.rows();
  const int ncols = data.cols();
  const int vector_len = (margin == 1) ? ncols : nrows;
  const int progress_n = (margin == 1) ? ncols * 2 : ncols;
  Progress p(progress_n, display_progress);
  const int *outer = data.outerIndexPtr();
  const int *inner = data.innerIndexPtr();
  double *x = data.valuePtr();
  bool has_nonfinite = false;

  if (margin == 1) {
    // Row-wise CLR without explicit transpose:
    // 1) accumulate row-wise log1p-sums from non-zero entries.
    // 2) convert each row sum to its geometric mean denominator.
    // 3) apply `log1p(value / denom[row])` on the original column-major matrix.
    std::vector<double> row_log_sums(nrows, 0.0);
    for (int k = 0; k < ncols; ++k) {
      p.increment();
      const int start = outer[k];
      const int end = outer[k + 1];
      for (int j = start; j < end; ++j) {
        const double value = x[j];
        if ((value > 0) && std::isfinite(value)) {
          row_log_sums[inner[j]] += log1p(value);
        } else if (!std::isfinite(value)) {
          has_nonfinite = true;
        }
      }
    }
    std::vector<double> row_denoms(nrows);
    const double inv_nrows = 1.0 / static_cast<double>(vector_len);
    for (int i = 0; i < nrows; ++i) {
      row_denoms[i] = exp(row_log_sums[i] * inv_nrows);
    }
    for (int k = 0; k < ncols; ++k) {
      p.increment();
      const int start = outer[k];
      const int end = outer[k + 1];
      for (int j = start; j < end; ++j) {
        const double value = x[j];
        if (has_nonfinite && !std::isfinite(value)) {
          continue;
        }
        x[j] = log1p(value / row_denoms[inner[j]]);
      }
    }
  } else {
    for (int k = 0; k < data.outerSize(); ++k) {
      p.increment();
      double log_sum = 0;
      const int start = outer[k];
      const int end = outer[k + 1];
      for (int j = start; j < end; ++j) {
        const double value = x[j];
        if ((value > 0) && std::isfinite(value)) {
          log_sum += log1p(value);
        } else if (!std::isfinite(value)) {
          has_nonfinite = true;
        }
      }
      const double denom = exp(log_sum / static_cast<double>(vector_len));
      for (int j = start; j < end; ++j) {
        const double value = x[j];
        if (has_nonfinite && !std::isfinite(value)) {
          continue;
        }
        x[j] = log1p(value / denom);
      }
    }
  }
  return data;
}

// Relative-counts normalize sparse columns in-place:
//   x <- x / colSums(x) * scale_factor.
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> RelativeCountsNorm(Eigen::SparseMatrix<double> data,
                                              double scale_factor,
                                              bool display_progress = true) {
  const int ncols = data.cols();
  const int *outer = data.outerIndexPtr();
  double *x = data.valuePtr();
  Progress p(ncols, display_progress);
  bool has_nonfinite = false;
  std::vector<double> col_sums(ncols, 0.0);
  for (int k = 0; k < ncols; ++k) {
    const int start = outer[k];
    const int end = outer[k + 1];
    double col_sum = 0.0;
    for (int j = start; j < end; ++j) {
      if (!std::isfinite(x[j])) {
        has_nonfinite = true;
      }
      col_sum += x[j];
    }
    col_sums[k] = col_sum;
  }
  for (int k = 0; k < ncols; ++k) {
    p.increment();
    const double col_sum = col_sums[k];
    const int start = outer[k];
    const int end = outer[k + 1];
    const double scale = scale_factor / col_sum;
    for (int j = start; j < end; ++j) {
      const double value = x[j];
      if (has_nonfinite && !std::isfinite(value)) {
        continue;
      }
      x[j] = value * scale;
    }
  }
  return data;
}

/* Performs column scaling and/or centering. Equivalent to using scale(mat, TRUE, apply(x,2,sd)) in R.
 Note: Doesn't handle NA/NaNs in the same way the R implementation does, */

// [[Rcpp::export(rng = false)]]
NumericMatrix Standardize(Eigen::Map<Eigen::MatrixXd> mat, bool display_progress = true){
  Progress p(mat.cols(), display_progress);
  NumericMatrix std_mat(mat.rows(), mat.cols());
  for(int i=0; i < mat.cols(); ++i){
    p.increment();
    Eigen::ArrayXd r = mat.col(i).array();
    double colMean = r.mean();
    double colSdev = sqrt((r - colMean).square().sum() / (mat.rows() - 1));
    NumericMatrix::Column new_col = std_mat(_, i);
    for(int j=0; j < new_col.size(); j++) {
      new_col[j] = (r[j] - colMean) / colSdev;
    }
  }
  return std_mat;
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastSparseRowScale(Eigen::SparseMatrix<double> mat, bool scale = true, bool center = true,
                                   double scale_max = 10, bool display_progress = true){
  mat = mat.transpose();
  Progress p(mat.outerSize(), display_progress);
  Eigen::MatrixXd scaled_mat(mat.rows(), mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double colMean = 0;
    double colSdev = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      colMean += it.value();
    }
    colMean = colMean / mat.rows();
    if (scale == true){
      int nnZero = 0;
      if(center == true){
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          nnZero += 1;
          colSdev += pow((it.value() - colMean), 2);
        }
        colSdev += pow(colMean, 2) * (mat.rows() - nnZero);
      }
      else{
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          colSdev += pow(it.value(), 2);
        }
      }
      colSdev = sqrt(colSdev / (mat.rows() - 1));
    }
    else{
      colSdev = 1;
    }
    if(center == false){
      colMean = 0;
    }
    Eigen::VectorXd col = Eigen::VectorXd(mat.col(k));
    scaled_mat.col(k) = (col.array() - colMean) / colSdev;
    for(int s=0; s<scaled_mat.col(k).size(); ++s){
      if(scaled_mat(s,k) > scale_max){
        scaled_mat(s,k) = scale_max;
      }
    }
  }
  return scaled_mat.transpose();
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastSparseRowScaleWithKnownStats(Eigen::SparseMatrix<double> mat, NumericVector mu, NumericVector sigma, bool scale = true, bool center = true,
                                   double scale_max = 10, bool display_progress = true){
    mat = mat.transpose();
    Progress p(mat.outerSize(), display_progress);
    Eigen::MatrixXd scaled_mat(mat.rows(), mat.cols());
    for (int k=0; k<mat.outerSize(); ++k){
        p.increment();
        double colMean = 0;
        double colSdev = 1;
        if (scale == true){
            colSdev = sigma[k];
        }
        if(center == true){
            colMean = mu[k];
        }
        Eigen::VectorXd col = Eigen::VectorXd(mat.col(k));
        scaled_mat.col(k) = (col.array() - colMean) / colSdev;
        for(int s=0; s<scaled_mat.col(k).size(); ++s){
            if(scaled_mat(s,k) > scale_max){
                scaled_mat(s,k) = scale_max;
            }
        }
    }
    return scaled_mat.transpose();
}

/* Note: May not handle NA/NaNs in the same way the R implementation does, */

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastCov(Eigen::MatrixXd mat, bool center = true){
  if (center) {
    mat = mat.rowwise() - mat.colwise().mean();
  }
  Eigen::MatrixXd cov = (mat.adjoint() * mat) / double(mat.rows() - 1);
  return(cov);
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastCovMats(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool center = true){
  if(center){
    mat1 = mat1.rowwise() - mat1.colwise().mean();
    mat2 = mat2.rowwise() - mat2.colwise().mean();
  }
  Eigen::MatrixXd cov = (mat1.adjoint() * mat2) / double(mat1.rows() - 1);
  return(cov);
}

/* Note: Faster than the R implementation but is not in-place */
// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastRBind(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2){
  Eigen::MatrixXd mat3(mat1.rows() + mat2.rows(), mat1.cols());
  mat3 << mat1, mat2;
  return(mat3);
}

/* Calculates the row means of the logged values in non-log space */
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd FastExpMean(Eigen::SparseMatrix<double> mat, bool display_progress){
  int ncols = mat.cols();
  Eigen::VectorXd rowmeans(mat.rows());
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene means" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += expm1(it.value());
    }
    rm = rm / ncols;
    rowmeans[k] = log1p(rm);
  }
  return(rowmeans);
}


/* use this if you know the row means */
// [[Rcpp::export(rng = false)]]
NumericVector SparseRowVar2(Eigen::SparseMatrix<double> mat,
                            NumericVector mu,
                            bool display_progress){
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene variances" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  NumericVector allVars = no_init(mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
      nZero -= 1;
      colSum += pow(it.value() - mu[k], 2);
    }
    colSum += pow(mu[k], 2) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}

/* standardize matrix rows using given mean and standard deviation,
   clip values larger than vmax to vmax,
   then return variance for each row */
// [[Rcpp::export(rng = false)]]
NumericVector SparseRowVarStd(Eigen::SparseMatrix<double> mat,
                              NumericVector mu,
                              NumericVector sd,
                              double vmax,
                              bool display_progress){
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating feature variances of standardized and clipped values" << std::endl;
  }
  mat = mat.transpose();
  NumericVector allVars(mat.cols());
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    if (sd[k] == 0) continue;
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      nZero -= 1;
      colSum += pow(std::min(vmax, (it.value() - mu[k]) / sd[k]), 2);
    }
    colSum += pow((0 - mu[k]) / sd[k], 2) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}

/* Calculate the variance to mean ratio (VMR) in non-logspace (return answer in
log-space) */
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd FastLogVMR(Eigen::SparseMatrix<double> mat,  bool display_progress){
  int ncols = mat.cols();
  Eigen::VectorXd rowdisp(mat.rows());
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene variance to mean ratios" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    double v = 0;
    int nnZero = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += expm1(it.value());
    }
    rm = rm / ncols;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      v += pow(expm1(it.value()) - rm, 2);
      nnZero += 1;
    }
    v = (v + (ncols - nnZero) * pow(rm, 2)) / (ncols - 1);
    rowdisp[k] = log(v/rm);

  }
  return(rowdisp);
}

/* Calculates the variance of rows of a matrix */
// [[Rcpp::export(rng = false)]]
NumericVector RowVar(Eigen::Map<Eigen::MatrixXd> x){
  NumericVector out(x.rows());
  for(int i=0; i < x.rows(); ++i){
    Eigen::ArrayXd r = x.row(i).array();
    double rowMean = r.mean();
    out[i] = (r - rowMean).square().sum() / (x.cols() - 1);
  }
  return out;
}

/* Calculate the variance in non-logspace (return answer in non-logspace) */
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd SparseRowVar(Eigen::SparseMatrix<double> mat, bool display_progress){
  int ncols = mat.cols();
  Eigen::VectorXd rowdisp(mat.rows());
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene variances" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    double v = 0;
    int nnZero = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += (it.value());
    }
    rm = rm / ncols;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      v += pow((it.value()) - rm, 2);
      nnZero += 1;
    }
    v = (v + (ncols - nnZero) * pow(rm, 2)) / (ncols - 1);
    rowdisp[k] = v;
  }
  return(rowdisp);
}

//cols_idx should be 0-indexed
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> ReplaceColsC(Eigen::SparseMatrix<double> mat, NumericVector col_idx, Eigen::SparseMatrix<double> replacement){
  int rep_idx = 0;
  for(auto const &ci : col_idx){
    mat.col(ci) = replacement.col(rep_idx);
    rep_idx += 1;
  }
  return(mat);
}

template <typename S>
std::vector<size_t> sort_indexes(const std::vector<S> &v) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}

// [[Rcpp::export(rng = false)]]
List GraphToNeighborHelper(Eigen::SparseMatrix<double> mat) {
  mat = mat.transpose();
  //determine the number of neighbors
  int n = 0;
  for(Eigen::SparseMatrix<double>::InnerIterator it(mat, 0); it; ++it) {
    n += 1;
  }
  Eigen::MatrixXd nn_idx(mat.rows(), n);
  Eigen::MatrixXd nn_dist(mat.rows(), n);

  for (int k=0; k<mat.outerSize(); ++k){
    int n_k = 0;
    std::vector<double> row_idx;
    std::vector<double> row_dist;
    row_idx.reserve(n);
    row_dist.reserve(n);
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
      if (n_k > (n-1)) {
        Rcpp::stop("Not all cells have an equal number of neighbors.");
      }
      row_idx.push_back(it.row() + 1);
      row_dist.push_back(it.value());
      n_k += 1;
    }
    if (n_k != n) {
      Rcpp::Rcout << n << ":::" << n_k << std::endl;
      Rcpp::stop("Not all cells have an equal number of neighbors.");
    }
    //order the idx based on dist
    std::vector<size_t> idx_order = sort_indexes(row_dist);
    for(int i = 0; i < n; ++i) {
      nn_idx(k, i) = row_idx[idx_order[i]];
      nn_dist(k, i) = row_dist[idx_order[i]];
    }
  }
  List neighbors = List::create(nn_idx, nn_dist);
  return(neighbors);
}
