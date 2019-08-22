#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <fst/fstlib.h>

using namespace fst;

VectorFst<StdArc> optimize(VectorFst<StdArc>);

#endif
