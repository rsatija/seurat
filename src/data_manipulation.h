#ifndef DATA_MANIPULATION
#define DATA_MANIPULATION

#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>

using namespace Rcpp;

//----------------------------------------------------
Eigen::SparseMatrix<double> RunUMISampling(Eigen::SparseMatrix<double> data, int sample_val,
                                           bool upsample, bool display_progress);
Eigen::SparseMatrix<double> RunUMISamplingPerCell(Eigen::SparseMatrix<double> data,
                                                  NumericVector sample_val, bool upsample,
                                                  bool display_progress);
Eigen::SparseMatrix<double> RowMergeMatrices(Eigen::SparseMatrix<double, Eigen::RowMajor> mat1,
                                             Eigen::SparseMatrix<double, Eigen::RowMajor> mat2,
                                             std::vector< std::string > mat1_rownames,
                                             std::vector< std::string > mat2_rownames,
                                             std::vector< std::string > all_rownames);
// Log-normalize sparse columns (used by NormalizeData sparse rewrite).
// Formula: value <- log1p(value / colSum(column) * scale_factor)
Eigen::SparseMatrix<double> LogNorm(Eigen::SparseMatrix<double> data, int scale_factor,
                                    bool display_progress );
// CLR-normalize sparse columns/rows:
//   x <- log1p(x / exp(mean(log1p(x[x > 0], na.rm=TRUE))))
// computed over each vector of the selected margin.
Eigen::SparseMatrix<double> CLRNorm(Eigen::SparseMatrix<double> data, int margin,
                                   bool display_progress );
// Relative-counts normalization on sparse columns:
//   x <- x / colSum(column) * scale_factor.
Eigen::SparseMatrix<double> RelativeCountsNorm(Eigen::SparseMatrix<double> data,
                                              double scale_factor,
                                              bool display_progress);
NumericMatrix Standardize(const Eigen::Map<Eigen::MatrixXd> mat, bool display_progress);
Eigen::MatrixXd FastSparseRowScale(Eigen::SparseMatrix<double> mat, bool scale, bool center,
                                   double scale_max, bool display_progress);
Eigen::MatrixXd FastCov(Eigen::MatrixXd mat, bool center);
Eigen::MatrixXd FastCovMats(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool center);
Eigen::MatrixXd FastRBind(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2);
Eigen::VectorXd FastExpMean(Eigen::MatrixXd mat, bool display_progress);
List SparseRowStats(Eigen::SparseMatrix<double> mat, bool display_progress);
Eigen::VectorXd FastRowMean(Eigen::MatrixXd mat, bool display_progress);
Eigen::VectorXd FastLogVMR(Eigen::SparseMatrix<double> mat, bool display_progress);
Eigen::VectorXd FastExpVar(Eigen::SparseMatrix<double> mat, bool display_progress);
Eigen::VectorXd SparseRowVar(Eigen::SparseMatrix<double> mat, bool display_progress);
NumericVector SparseRowVar2(Eigen::SparseMatrix<double> mat,
                            NumericVector mu,
                            bool display_progress);
NumericVector SparseRowVarStd(Eigen::SparseMatrix<double> mat,
                              NumericVector mu,
                              NumericVector sd,
                              double vmax,
                              bool display_progress);
NumericVector RowVar(Eigen::Map<Eigen::MatrixXd> x);
template <typename S>
std::vector<size_t> sort_indexes(const std::vector<S> &v);
List GraphToNeighborHelper(Eigen::SparseMatrix<double> mat);
//----------------------------------------------------

#endif//DATA_MANIPULATION
