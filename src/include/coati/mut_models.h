#ifndef MUT_MODELS_H
#define MUT_MODELS_H

#include <fst/fstlib.h>
#include <Eigen/Dense>
#include <iostream>
#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>
#include <unsupported/Eigen/MatrixFunctions>
#include "utils.h"

typedef Eigen::Matrix<double, 64, 64>Matrix64f;
typedef Eigen::Matrix<double, 64, 1>Vector64f;

using namespace fst;
using namespace std;

void nuc2pos(VectorFst<StdArc>& n2p);
void marg_mut(VectorFst<StdArc>& mut_fst, VectorFst<StdArc> marg_pos);
void toycoati(VectorFst<StdArc>& mut_fst);
void toy_marg(VectorFst<StdArc>& mut_fst);
void dna_mut(VectorFst<StdArc>& mut_fst);
void ecm_p(Matrix64f& P);
void ecm(VectorFst<StdArc>& mut_fst);
void ecm_marginal(VectorFst<StdArc>& mut_fst);
bool syn(cod c1, cod c2);
void nts_ntv(cod c1, cod c2, int& nts, int& ntv);
double k(cod c1, cod c2, int model=0);

#endif
