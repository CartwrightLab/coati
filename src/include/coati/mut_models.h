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

void nuc2pos(fst::VectorFst<fst::StdArc>& n2p);
void marg_mut(fst::VectorFst<fst::StdArc>& mut_fst, fst::VectorFst<fst::StdArc> marg_pos);
void toycoati(fst::VectorFst<fst::StdArc>& mut_fst);
void dna_mut(fst::VectorFst<fst::StdArc>& mut_fst);
void ecm_p(Matrix64f& P);
void ecm(fst::VectorFst<fst::StdArc>& mut_fst);
void ecm_marginal(fst::VectorFst<fst::StdArc>& ecm_m);
std::string cod2aa();
bool syn();
double k();
bool isStop();

#endif
