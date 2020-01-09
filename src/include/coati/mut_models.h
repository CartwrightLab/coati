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
#include "optimize.h"
#include "utils.h"

typedef Eigen::Matrix<double, 64, 64>Matrix64f;
typedef Eigen::Vector<double, 64>Vector64f;

fst::VectorFst<fst::StdArc> marg_mut();
fst::VectorFst<fst::StdArc> nuc2pos();
fst::VectorFst<fst::StdArc> toycoati();
fst::VectorFst<fst::StdArc> dna_mut();
fst::VectorFst<fst::StdArc> ecm();
std::string cod2aa();
bool syn();
float k();
bool isStop();

#endif
