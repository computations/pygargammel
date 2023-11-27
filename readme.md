# Introduction

PyGargammel is a reimplementation the ancient dna damage simulator GARGAMMEL focused on better exposing the model
parameters for simulation. To this end, it is quite barebones, and only implements the Brigg's model of ADNA damage.
Unlike GARGAMMEL, PyGargammel does not do contamination of samples, instead choosing to adopt a unix philosophy of doing
just one thing. In this case, it is implementing and exposing the parameters of the Brigg's Model.

# Usage

An example command for PyGargammel is 

```
pygargammel --fasta align.fasta --nick-freq 0.24 --overhang-parameter 0.36 --double-strand-deamination 0.0097
--single-strand-deamination 0.22 --output damage.fasta --log pygargammel.log
```

Each of the options specified in the example command are required. Here is a quick and informal explanation of the
options:

- `--fasta`: The source sequences to damage. If more than one sequence is present PyGargammel will attempt to damage
  each sequence individually.
- `--nick-freq`: Probability of nicks per base. Has the effect of controlling the length and number of fragments.
  Higher numbers mean more and shorter fragments. Accepts values between 0.0 and 1.0, where 0.0 is no chance of nicks,
  and 1.0 is a nick on every base.
- `--overhang-parameter`: Parameter to the geometric distribution determining the length of the overhangs. Somewhat
  counter intuitively, higher numbers mean shorter overhangs. Accepts values in the range `(0.0, 1.0]`. Value close to
  zero give very long overhangs, and 1.0 will give zero length overhangs.
- `--double-strand-deamination`: Rate at deamination will occur in double stranded regions. Accepts values between 0.0
  and 1.0, with 0.0 being no deamination, and 1.0 being total deamination.
- `--single-strand-deamination`: Rate at deamination will occur in single strand regions. Accepts values between 0.0
  and 1.0, with 0.0 being no deamination, and 1.0 being total deamination.
- `--reference`: Ensure that generated fragments have at least _one_ informative
  character in common with at least _one_ sequence in the alignment.
