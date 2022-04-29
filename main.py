#!/usr/bin/env python3

import numpy
import enum
import itertools
import json
import gzip
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq


class Strand(enum.Enum):
    POSITIVE = enum.auto()
    NEGATIVE = enum.auto()


class DnaNick:

    def __init__(self, index, strand=Strand.POSITIVE):
        self._index = index
        self._strand = strand

    def index(self):
        return self._index

    def __str__(self):
        return f"(index: {self._index})"


class SiteDamage:

    def __init__(self, index, old_char, new_char):
        self._relative_index = index
        self._old_character = old_char
        self._new_character = new_char

    def apply_to_seq(self, mutseq):
        if mutseq[self._relative_index] == self._old_character:
            mutseq[self._relative_index] = self._new_character
        else:
            raise RuntimeError(
                f"The old character ({self._old_character}) did not match the actual character ({mutseq[self._relative_index]})"
            )

    def __str__(self):
        return f"(index: {self._relative_index}, old: {self._old_character}, new: {self._new_character})"


class AdnaFragment:

    def __init__(self, sequence, start, end, sequence_name):
        self._start_index = start
        self._end_index = end
        self._sequence = sequence[self._start_index:self._end_index]
        self._sequence_name = sequence_name
        self._damage = []

    @property
    def sequence_name(self):
        return self._sequence_name

    @property
    def annotated_sequence_name(self):
        return self.sequence_name + "_" + str(self._start_index) + "-" + str(
            self._end_index)

    @property
    def left_overhang(self):
        return self._overhang_left_len

    @left_overhang.setter
    def left_overhang(self, length):
        if length > len(self._sequence):
            raise RuntimeError("Overhang is too long")
        self._overhang_left_len = length

    @property
    def right_overhang(self):
        return self._overhang_right_len

    @right_overhang.setter
    def right_overhang(self, length):
        if length > len(self._sequence):
            raise RuntimeError("Overhang is too long")
        self._overhang_right_len = length

    @property
    def right_overhang_start_index(self):
        return len(self) - self.right_overhang

    def __iter__(self, index):
        return self._sequence[index]

    def __len__(self):
        return len(self._sequence)

    def add_damage(self, damage):
        self._damage.append(damage)

    def _to_sequence_mut(self):
        mut_seq = MutableSeq(self._sequence)
        for d in self._damage:
            d.apply_to_seq(mut_seq)
        return mut_seq

    def to_sequence(self):
        return Seq(self._to_sequence_mut())

    def to_aligned_sequence(self, sequence_length):
        mut_seq = self._to_sequence_mut()
        return ('N' * self._start_index) + mut_seq + (
            'N' * (sequence_length - (len(mut_seq) + self._start_index)))

    def __str__(self):
        return f"(start: {self._start_index}, end: {self._end_index}, left overhang: {self._overhang_left_len}, right overhang: {self._overhang_right_len})"

    def add_overhang(self, overhang_length_param):
        sequence_len = len(self._sequence)
        left_overhang_len = sequence_len
        right_overhang_len = sequence_len

        while (left_overhang_len > sequence_len
               and right_overhang_len > sequence_len
               and right_overhang_len + left_overhang_len > sequence_len):
            left_overhang_len = numpy.random.geometric(
                overhang_length_param) - 1
            right_overhang_len = numpy.random.geometric(
                overhang_length_param) - 1

        self.left_overhang = left_overhang_len
        self.right_overhang = right_overhang_len

    def _roll_damage(self, index, seq_char, test_char, damage_char, rate):
        if seq_char.lower() == test_char.lower():
            if numpy.random.uniform() < rate:
                tmp_damage = SiteDamage(index, seq_char, damage_char)
                self.add_damage(tmp_damage)

    def _add_c_to_t_damage(self, single_strand_rate, double_strand_rate):
        for i, c in enumerate(self._sequence[:self.left_overhang]):
            self._roll_damage(i, c, 'C', 'T', single_strand_rate)

        for i, c in enumerate(self._sequence[self.left_overhang:],
                              start=self.left_overhang):
            self._roll_damage(i, c, 'C', 'T', double_strand_rate)

    def _add_a_to_g_damage(self, single_strand_rate, double_strand_rate):
        for i, c in enumerate(
                self._sequence[:self.right_overhang_start_index]):
            self._roll_damage(i, c, 'A', 'G', double_strand_rate)

        for i, c in enumerate(self._sequence[self.right_overhang_start_index:],
                              start=self.right_overhang_start_index):
            self._roll_damage(i, c, 'A', 'G', single_strand_rate)

    def add_deaminations(self, single_strand_rate, double_strand_rate):
        # add left overhand damage
        self._add_c_to_t_damage(single_strand_rate, double_strand_rate)
        self._add_a_to_g_damage(single_strand_rate, double_strand_rate)

    def to_dict(self):
        obj_dict = vars(self)
        obj_dict.pop("_sequence", None)
        return obj_dict


def produce_nicks(sequence, nick_freq):
    nicks = []
    for i, _ in enumerate(sequence):
        if numpy.random.uniform() < nick_freq:
            nicks.append(DnaNick(i))
    return nicks


def produce_fragments(nicks, sequence, sequence_name):
    fragments = []
    fragments.append(AdnaFragment(sequence, 0, nicks[0].index(),
                                  sequence_name))

    for n1, n2 in itertools.pairwise(nicks):
        fragments.append(
            AdnaFragment(sequence, n1.index(), n2.index(), sequence_name))

    fragments.append(
        AdnaFragment(sequence, nicks[-1].index(), len(sequence),
                     sequence_name))
    return fragments


def filter_fragments(frags, min_len=0):
    ret_frags = []
    for frag in frags:
        if len(frag) < min_len:
            continue
        if len(frag.to_sequence().ungap()) == 0:
            continue
        ret_frags.append(frag)

    return ret_frags


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=str, required=True)
    parser.add_argument("--nick-freq", type=float, required=True)
    parser.add_argument("--overhang-parameter", type=float, required=True)
    parser.add_argument("--min-length", type=int, required=False, default=None)
    parser.add_argument("--min-fragments", type=int, required=False, default=0)

    parser.add_argument("--double-strand-deamination",
                        type=float,
                        required=True)
    parser.add_argument("--single-strand-deamination",
                        type=float,
                        required=True)

    parser.add_argument("--align",
                        action=argparse.BooleanOptionalAction,
                        default=False)
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--log", type=str, required=True)
    parser.add_argument("--seed", type=int)
    args = parser.parse_args()

    alignment = SeqIO.parse(open(args.fasta), 'fasta')

    if not args.seed is None:
        numpy.random.default_rng(args.seed)

    with gzip.GzipFile(args.log + '.gz', 'w') as logfile:
        for sequence in alignment:
            total_fragments = 0
            while total_fragments < args.min_fragments:
                nicks = produce_nicks(sequence.seq, args.nick_freq)
                frags = produce_fragments(nicks, sequence.seq, sequence.id)
                for frag in frags:
                    frag.add_overhang(args.overhang_parameter)
                    frag.add_deaminations(args.single_strand_deamination,
                                          args.double_strand_deamination)

                frags = filter_fragments(
                    frags,
                    args.min_length if not args.min_length is None else 0)
                total_fragments += len(frags)
                with gzip.GzipFile(args.output + ".gz", 'a') as fastafile:
                    for f in frags:
                        fastafile.write((">" + f.annotated_sequence_name +
                                         "\n").encode('utf-8'))
                        fastafile.write(
                            bytes(
                                f.to_aligned_sequence(len(sequence)) if args.
                                align else f.to_sequence()))
                        fastafile.write(b"\n")
                logfile.write(
                    json.dumps([f.to_dict() for f in frags],
                               default=vars).encode('utf-8'))
