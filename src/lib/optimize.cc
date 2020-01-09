
#include <fst/fstlib.h>

using namespace fst;

// TODO: write as template for more fst and arc type variability
VectorFst<StdArc> optimize(VectorFst<StdArc> fst_raw) {

	// encode FST
	SymbolTable *syms = SymbolTable::ReadText("fst/dna_syms.txt");
	EncodeMapper<StdArc> encoder(kEncodeLabels, ENCODE);
	encoder.SetInputSymbols(syms);
	encoder.SetOutputSymbols(syms);
	EncodeFst<StdArc> fst_enc = EncodeFst<StdArc>(fst_raw, encoder);

	// reduce to more efficient form
	// 1. epsilon removal
	VectorFst<StdArc> fst_rmep;
	fst_rmep = RmEpsilonFst<StdArc>(fst_enc);

	// 2. determinize
	VectorFst<StdArc> fst_det;
	fst_det = DeterminizeFst<StdArc>(fst_rmep);

	// 3. minimize
	VectorFst<StdArc> fst_min;
	Minimize(&fst_det, &fst_min);

	// decode
	DecodeFst<StdArc> fst_dec = DecodeFst<StdArc>(fst_det, encoder);

	// epsilon removal after decoding
	VectorFst<StdArc> fst_opt;
	fst_opt = RmEpsilonFst<StdArc>(fst_dec);

	return fst_opt;
}
