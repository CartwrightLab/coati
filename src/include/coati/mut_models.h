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
void marg_mut(fst::VectorFst<fst::StdArc>& mut_fst);
void toycoati(fst::VectorFst<fst::StdArc>& mut_fst);
void dna_mut(fst::VectorFst<fst::StdArc>& mut_fst);
void ecm(fst::VectorFst<fst::StdArc>& mut_fst);
std::string cod2aa();
bool syn();
float k();
bool isStop();

#endif
