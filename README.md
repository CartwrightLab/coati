# coati
Codon-Aware Multiple Sequence Alignments

## Pairwise Alignment: coati-alignpair

```
Usage:	coati-alignpair file.fasta [options]

Allowed options:
  -h [ --help ]                 Display this message
  -f [ --fasta ] arg            fasta file path
  -m [ --model ] arg (=m-coati) substitution model: coati, m-coati (default),
                                dna, ecm, m-ecm
  -w [ --weight ] arg           Write alignment score to file
  -o [ --output ] arg           Alignment output file
  -s [ --score ]                Calculate alignment score using marginal COATi
                                model
  -r [ --rate ] arg             Substitution rate matrix (CSV)
 ```
