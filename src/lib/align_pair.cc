#include <doctest/doctest.h>

#include <coati/align_pair.hpp>

namespace coati {

// m -> m      (1-g)*(1-g)*P(b[j] | a[i])/P(b[j])
// m -> d      (1-g)*g
// m -> i      g
// m -> END    1-g
// d -> m      1-e
// d -> d      e
// d -> i      0
// d -> END    1
// i -> m      (1-e)*(1-g)
// i -> d      (1-e)*g
// i -> i      e
// i -> END    (1-e)

// a is "ancestor"
// b is "descendant"
// match matrix holds log(P(b[j]|a[i]))-log(P(b[j]))
// gap_open = log(g)
// gap_extend = log(e)
void align_pair(align_pair_work_t &work, const seq_view_t &a,
                const seq_view_t &b, const Matrixf &match,
                utils::args_t &args) {
    // calculate log(1-g) log(1-e) log(g) log(e)
    float_t no_gap = std::log1pf(-args.gap.open);
    float_t gap_stop = std::log1pf(-args.gap.extend);
    float_t gap_open = ::logf(args.gap.open);
    float_t gap_extend = ::logf(args.gap.extend);
    size_t look_back = args.gap.len;
    size_t start = look_back - 1;

    const float_t lowest = std::numeric_limits<float_t>::lowest();

    // create matrices
    size_t len_a = a.length() + look_back;  // length of ancestor
    size_t len_b = b.length() + look_back;  // length of descendant
    work.resize(len_a, len_b, lowest);      // resize work matrices

    // initialize the margins of the matrices
    work.mch(start, start) = 0.0;

    for(size_t i = start + look_back; i < len_a; i += look_back) {
        work.del(i, start) = work.del_del(i, start) =
            no_gap + gap_open + gap_extend * static_cast<float_t>(i - 1);
    }
    for(size_t j = start + look_back; j < len_b; j += look_back) {
        work.ins(start, j) = work.ins_ins(start, j) =
            gap_open + gap_extend * static_cast<float_t>(j - 1);
    }

    // fill the body of the matrices
    for(size_t i = look_back; i < len_a; ++i) {
        for(size_t j = look_back; j < len_b; ++j) {
            //  from match, ins, or del to match
            auto mch = match(a[i - look_back], b[j - look_back]);
            work.mch_mch(i, j) = work.mch(i - 1, j - 1) + 2 * no_gap + mch;
            work.del_mch(i, j) = work.del(i - 1, j - 1) + gap_stop + mch;
            work.ins_mch(i, j) =
                work.ins(i - 1, j - 1) + gap_stop + no_gap + mch;

            // from match or del to del
            work.mch_del(i, j) =
                work.mch(i - look_back, j) + no_gap + gap_open +
                gap_extend * static_cast<float_t>(look_back - 1);
            work.ins_del(i, j) =
                work.ins(i - look_back, j) + gap_stop + gap_open +
                gap_extend * static_cast<float_t>(look_back - 1);
            work.del_del(i, j) =
                work.del(i - look_back, j) + gap_extend * look_back;

            // from match, del, or ins to ins
            work.mch_ins(i, j) =
                work.mch(i, j - look_back) + gap_open +
                gap_extend * static_cast<float_t>(look_back - 1);
            work.ins_ins(i, j) =
                work.ins(i, j - look_back) +
                gap_extend * static_cast<float_t>(look_back - 1);

            // save score
            work.mch(i, j) = maximum(work.mch_mch(i, j), work.del_mch(i, j),
                                     work.ins_mch(i, j));
            work.del(i, j) = maximum(work.mch_del(i, j), work.del_del(i, j),
                                     work.ins_del(i, j));
            work.ins(i, j) = maximum(work.mch_ins(i, j), work.ins_ins(i, j));
        }
    }
    {
        // adjust the terminal state
        work.mch(len_a - 1, len_b - 1) += no_gap;
        work.ins(len_a - 1, len_b - 1) += gap_stop;
    }
}

void traceback(const align_pair_work_t &work, const std::string &a,
               const std::string &b, utils::alignment_t &aln,
               size_t look_back) {
    size_t i = static_cast<int>(work.mch_mch.rows() - 1);
    size_t j = static_cast<int>(work.mch_mch.cols() - 1);

    aln.fasta.seqs.resize(2);
    aln.fasta.seqs[0].reserve(i + j);
    aln.fasta.seqs[1].reserve(i + j);

    aln.weight = maximum(work.mch(i, j), work.del(i, j), work.ins(i, j));
    int m = max_index(work.mch(i, j), work.del(i, j), work.ins(i, j));

    while((j > (look_back - 1)) || (i > (look_back - 1))) {
        switch(m) {
        case 0:  // match
            aln.fasta.seqs[0].push_back(a[i - look_back]);
            aln.fasta.seqs[1].push_back(b[j - look_back]);
            m = max_index(work.mch_mch(i, j), work.del_mch(i, j),
                          work.ins_mch(i, j));
            i--;
            j--;
            break;
        case 1:  // deletion
            for(size_t k = i; k > (i - look_back); k--) {
                aln.fasta.seqs[0].push_back(a[k - look_back]);
                aln.fasta.seqs[1].push_back('-');
            }
            m = max_index(work.mch_del(i, j), work.del_del(i, j),
                          work.ins_del(i, j));
            i -= look_back;
            break;
        case 2:  // insertion
            for(size_t k = j; k > (j - look_back); k--) {
                aln.fasta.seqs[0].push_back('-');
                aln.fasta.seqs[1].push_back(b[k - look_back]);
            }
            m = work.mch_ins(i, j) > work.ins_ins(i, j) ? 0 : 2;
            j -= look_back;
            break;
        }
    }

    std::reverse(aln.fasta.seqs[0].begin(), aln.fasta.seqs[0].end());
    std::reverse(aln.fasta.seqs[1].begin(), aln.fasta.seqs[1].end());
}

}  // namespace coati
