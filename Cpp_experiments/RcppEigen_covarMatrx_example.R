
##########################################################################


   fcprd <- cxxfunction(signature(AA = "matrix"),
                        'using Eigen::Map;
                         using Eigen::MatrixXd;
                         using Eigen::Lower;

                         const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
                         const int m(A.rows()), n(A.cols());
                         MatrixXd AAt(MatrixXd(m, m).setZero().selfadjointView<Lower>().rankUpdate(A));

                         return List::create(Named("crossprod(A)") = AAt);'
                        , "RcppEigen")
  str(crp <- fcprd(A))
  fcprd(x)

