"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread

This module contains some tests for Badread. To run them, execute `python3 -m unittest` from the
root Badread directory.

This file is part of Badread. Badread is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Badread is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Badread.
If not, see <http://www.gnu.org/licenses/>.
"""

import collections
import math
import os
import random
import statistics
import unittest

import badread.misc
import badread.qscore_model
import badread.settings


class TestAlignSequences(unittest.TestCase):

    def test_align_sequences_1(self):
        in_seq_1 = 'ACGTCGTA'
        in_seq_2 = 'ACGTCGTA'
        in_cigar = '8='
        out_seq_1, out_seq_2, out_cigar = \
            badread.qscore_model.align_sequences_from_edlib_cigar(in_seq_1, in_seq_2, in_cigar)
        self.assertEqual(out_seq_1, 'ACGTCGTA')
        self.assertEqual(out_cigar, '========')
        self.assertEqual(out_seq_2, 'ACGTCGTA')

    def test_align_sequences_2(self):
        in_seq_1 = 'ACGTCGTA'
        in_seq_2 = 'ACGTGGTA'
        in_cigar = '4=1X3='
        out_seq_1, out_seq_2, out_cigar = \
            badread.qscore_model.align_sequences_from_edlib_cigar(in_seq_1, in_seq_2, in_cigar)
        self.assertEqual(out_seq_1, 'ACGTCGTA')
        self.assertEqual(out_cigar, '====X===')
        self.assertEqual(out_seq_2, 'ACGTGGTA')

    def test_align_sequences_3(self):
        in_seq_1 = 'ACGTCGTA'
        in_seq_2 = 'ACGTGTA'
        in_cigar = '4=1I3='
        out_seq_1, out_seq_2, out_cigar = \
            badread.qscore_model.align_sequences_from_edlib_cigar(in_seq_1, in_seq_2, in_cigar)
        self.assertEqual(out_seq_1, 'ACGTCGTA')
        self.assertEqual(out_cigar, '====I===')
        self.assertEqual(out_seq_2, 'ACGT-GTA')

    def test_align_sequences_4(self):
        in_seq_1 = 'ACGTCGTA'
        in_seq_2 = 'ACGTGTA'
        in_cigar = '4=1I3='
        out_seq_1, out_seq_2, out_cigar = \
            badread.qscore_model.align_sequences_from_edlib_cigar(in_seq_1, in_seq_2, in_cigar,
                                                                  gap_char=' ')
        self.assertEqual(out_seq_1, 'ACGTCGTA')
        self.assertEqual(out_cigar, '====I===')
        self.assertEqual(out_seq_2, 'ACGT GTA')

    def test_align_sequences_5(self):
        in_seq_1 = 'ACGTGTA'
        in_seq_2 = 'ACGTCGTA'
        in_cigar = '4=1D3='
        out_seq_1, out_seq_2, out_cigar = \
            badread.qscore_model.align_sequences_from_edlib_cigar(in_seq_1, in_seq_2, in_cigar)
        self.assertEqual(out_seq_1, 'ACGT-GTA')
        self.assertEqual(out_cigar, '====D===')
        self.assertEqual(out_seq_2, 'ACGTCGTA')


class TestUniformDist(unittest.TestCase):

    def test_uniform_dist_1(self):
        scores, probs = badread.qscore_model.uniform_dist_scores_and_probs(1, 2)
        self.assertEqual(scores, [1, 2])
        self.assertEqual(probs, [0.5, 0.5])

    def test_uniform_dist_2(self):
        scores, probs = badread.qscore_model.uniform_dist_scores_and_probs(10, 14)
        self.assertEqual(scores, [10, 11, 12, 13, 14])
        self.assertEqual(probs, [0.2, 0.2, 0.2, 0.2, 0.2])


class TestQScoreConversions(unittest.TestCase):

    def test_qscore_char_to_val(self):
        self.assertEqual(badread.qscore_model.qscore_char_to_val('!'), 0)
        self.assertEqual(badread.qscore_model.qscore_char_to_val('+'), 10)
        self.assertEqual(badread.qscore_model.qscore_char_to_val('5'), 20)
        self.assertEqual(badread.qscore_model.qscore_char_to_val('?'), 30)
        self.assertEqual(badread.qscore_model.qscore_char_to_val('I'), 40)

    def test_qscore_val_to_char(self):
        self.assertEqual(badread.qscore_model.qscore_val_to_char(0), '!')
        self.assertEqual(badread.qscore_model.qscore_val_to_char(10), '+')
        self.assertEqual(badread.qscore_model.qscore_val_to_char(20), '5')
        self.assertEqual(badread.qscore_model.qscore_val_to_char(30), '?')
        self.assertEqual(badread.qscore_model.qscore_val_to_char(40), 'I')


class TestRandomModel(unittest.TestCase):

    def setUp(self):
        self.trials = 10000
        null = open(os.devnull, 'w')
        self.model = badread.qscore_model.QScoreModel('random', output=null)
        null.close()

    def test_type(self):
        self.assertEqual(self.model.type, 'random')

    def test_scores_and_probs_cigars(self):
        self.assertEqual(sorted(self.model.scores.keys()), sorted(['=', 'X', 'I']))
        self.assertEqual(sorted(self.model.probabilities.keys()), sorted(['=', 'X', 'I']))

    def test_all_cigars_equal(self):
        self.assertEqual(self.model.scores['='], self.model.scores['X'])
        self.assertEqual(self.model.scores['='], self.model.scores['I'])
        self.assertEqual(self.model.probabilities['='], self.model.probabilities['X'])
        self.assertEqual(self.model.probabilities['='], self.model.probabilities['I'])

    def test_flat_distribution(self):
        self.assertEqual(len(set(self.model.probabilities['='])), 1)

    def test_mean_and_stdev(self):
        qscores = []
        for _ in range(self.trials):
            cigar = random.choice(['=', 'X', 'I'])  # cigar doesn't matter for random qscore model
            q = self.model.get_qscore(cigar)
            q = badread.qscore_model.qscore_char_to_val(q)
            qscores.append(q)
        dist_min = badread.settings.RANDOM_QSCORE_MIN
        dist_max = badread.settings.RANDOM_QSCORE_MAX
        target_mean = (dist_min + dist_max) / 2
        target_stdev = math.sqrt(((dist_max - dist_min + 1) ** 2 - 1) / 12)
        self.assertAlmostEqual(statistics.mean(qscores), target_mean, delta=0.5)
        self.assertAlmostEqual(statistics.stdev(qscores), target_stdev, delta=0.5)


class TestIdealModel(unittest.TestCase):

    def setUp(self):
        self.trials = 10000
        null = open(os.devnull, 'w')
        self.model = badread.qscore_model.QScoreModel('ideal', output=null)
        null.close()

    def test_type(self):
        self.assertEqual(self.model.type, 'ideal')

    def test_scores_and_probs_cigars(self):
        self.assertEqual(sorted(self.model.scores.keys()),
                         sorted(['=========', '=======', '=====', '===', '=', 'X', 'I']))
        self.assertEqual(sorted(self.model.probabilities.keys()),
                         sorted(['=========', '=======', '=====', '===', '=', 'X', 'I']))

    def one_cigar_test(self, cigar, dist_min, dist_max):
        qscores = []
        for _ in range(self.trials):
            q = self.model.get_qscore(cigar)
            q = badread.qscore_model.qscore_char_to_val(q)
            qscores.append(q)
        target_mean = (dist_min + dist_max) / 2
        target_stdev = math.sqrt(((dist_max - dist_min + 1) ** 2 - 1) / 12)
        self.assertAlmostEqual(statistics.mean(qscores), target_mean, delta=0.5)
        self.assertAlmostEqual(statistics.stdev(qscores), target_stdev, delta=0.5)

    def test_nine_matches_cigar(self):
        self.one_cigar_test('=========', badread.settings.IDEAL_QSCORE_RANK_6_MIN,
                            badread.settings.IDEAL_QSCORE_RANK_6_MAX)

    def test_seven_matches_cigar(self):
        self.one_cigar_test('=======', badread.settings.IDEAL_QSCORE_RANK_5_MIN,
                            badread.settings.IDEAL_QSCORE_RANK_5_MAX)

    def test_five_matches_cigar(self):
        self.one_cigar_test('=====', badread.settings.IDEAL_QSCORE_RANK_4_MIN,
                            badread.settings.IDEAL_QSCORE_RANK_4_MAX)

    def test_three_matches_cigar(self):
        self.one_cigar_test('===', badread.settings.IDEAL_QSCORE_RANK_3_MIN,
                            badread.settings.IDEAL_QSCORE_RANK_3_MAX)

    def test_one_match_cigar(self):
        self.one_cigar_test('=', badread.settings.IDEAL_QSCORE_RANK_2_MIN,
                            badread.settings.IDEAL_QSCORE_RANK_2_MAX)

    def test_insertion_cigar(self):
        self.one_cigar_test('I', badread.settings.IDEAL_QSCORE_RANK_1_MIN,
                            badread.settings.IDEAL_QSCORE_RANK_1_MAX)

    def test_mismatch_cigar(self):
        self.one_cigar_test('X', badread.settings.IDEAL_QSCORE_RANK_1_MIN,
                            badread.settings.IDEAL_QSCORE_RANK_1_MAX)


class TestLoadQScoreModel(unittest.TestCase):
    """
    Loads a simple qscore error model from a file and makes sure it looks okay.
    """
    def setUp(self):
        self.null = open(os.devnull, 'w')

    def tearDown(self):
        self.null.close()

    def test_k_size(self):
        model_filename = os.path.join(os.path.dirname(__file__), 'simple_qscore_model')
        model = badread.qscore_model.QScoreModel(model_filename, output=self.null)
        self.assertEqual(model.kmer_size, 3)

    def test_scores_and_probs_cigars(self):
        model_filename = os.path.join(os.path.dirname(__file__), 'simple_qscore_model')
        model = badread.qscore_model.QScoreModel(model_filename, output=self.null)
        self.assertEqual(sorted(model.scores.keys()),
                         sorted(model.probabilities.keys()))
        self.assertEqual(sorted(model.scores.keys()),
                         sorted(['=', 'I', 'X', '===', '=D==', '==D=', 'I==', '==I', '=I=',
                                 '==X', 'X==', '=X=', 'III', '=II', 'II=', '=DD==', '==DD=',
                                 'XX=', '=XX', 'X=X', '=D=D=', 'I=I', '=DDD==', '==DDD=', 'I=X',
                                 'X=I', 'X=D=', '=D=X', '=XI', 'IX=', 'XI=', '=IX']))

    def test_type(self):
        model_filename = os.path.join(os.path.dirname(__file__), 'simple_qscore_model')
        model = badread.qscore_model.QScoreModel(model_filename, output=self.null)
        self.assertEqual(model.type, 'model')

    def test_bad_format(self):
        model_filename = os.path.join(os.path.dirname(__file__), 'simple_qscore_model_bad_format')
        with self.assertRaises(SystemExit) as cm:
            badread.qscore_model.QScoreModel(model_filename, output=self.null)
        self.assertTrue('does not seem to be a valid qscore model file' in str(cm.exception))


class TestLoadBuiltInModels(unittest.TestCase):
    """
    Loads the Nanopore and PacBio models that come with Badread.
    """
    def setUp(self):
        self.null = open(os.devnull, 'w')

    def tearDown(self):
        self.null.close()

    def test_nanopore2018(self):
        model = badread.qscore_model.QScoreModel('nanopore2018', output=self.null)
        self.assertGreater(len(model.scores), 1000)
        self.assertEqual(model.kmer_size, 9)
        self.assertEqual(sorted(model.scores.keys()), sorted(model.probabilities.keys()))

    def test_nanopore2020(self):
        model = badread.qscore_model.QScoreModel('nanopore2020', output=self.null)
        self.assertGreater(len(model.scores), 1000)
        self.assertEqual(model.kmer_size, 9)
        self.assertEqual(sorted(model.scores.keys()), sorted(model.probabilities.keys()))

    def test_nanopore2023(self):
        model = badread.qscore_model.QScoreModel('nanopore2023', output=self.null)
        self.assertGreater(len(model.scores), 1000)
        self.assertEqual(model.kmer_size, 9)
        self.assertEqual(sorted(model.scores.keys()), sorted(model.probabilities.keys()))

    def test_pacbio2016(self):
        model = badread.qscore_model.QScoreModel('pacbio2016', output=self.null)
        self.assertGreater(len(model.scores), 1000)
        self.assertEqual(model.kmer_size, 9)
        self.assertEqual(sorted(model.scores.keys()), sorted(model.probabilities.keys()))

    def test_pacbio2021(self):
        model = badread.qscore_model.QScoreModel('pacbio2021', output=self.null)
        self.assertGreater(len(model.scores), 400)
        self.assertEqual(model.kmer_size, 9)
        self.assertEqual(sorted(model.scores.keys()), sorted(model.probabilities.keys()))


class TestLoadedQScoreModelDistributions(unittest.TestCase):
    """
    Loads a simple qscore error model and tests to make sure the distributions which come from
    particular CIGARs look okay.
    """
    def setUp(self):
        self.trials = 10000
        null = open(os.devnull, 'w')
        model_filename = os.path.join(os.path.dirname(__file__), 'simple_qscore_model')
        self.model = badread.qscore_model.QScoreModel(model_filename, output=null)
        null.close()

    def get_mean_stdev(self, cigar):
        qscores = []
        for _ in range(self.trials):
            q = self.model.get_qscore(cigar)
            q = badread.qscore_model.qscore_char_to_val(q)
            qscores.append(q)
        return statistics.mean(qscores), statistics.stdev(qscores)

    def test_cigar_1(self):
        mean, stdev = self.get_mean_stdev('=')
        self.assertAlmostEqual(mean, 10.97, delta=0.5)
        self.assertAlmostEqual(stdev, 3.22, delta=0.5)

    def test_cigar_2(self):
        mean, stdev = self.get_mean_stdev('=X=')
        self.assertAlmostEqual(mean, 6.64, delta=0.5)
        self.assertAlmostEqual(stdev, 4.26, delta=0.5)

    def test_cigar_3(self):
        mean, stdev = self.get_mean_stdev('==DDD=')
        self.assertAlmostEqual(mean, 9.12, delta=0.5)
        self.assertAlmostEqual(stdev, 3.76, delta=0.5)


class TestGetQScoresRandom(unittest.TestCase):
    """
    For the random qscore model, it shouldn't matter what the sequence is, the best and worst
    qscores can appear in any place.
    """
    def setUp(self):
        self.trials = 1000
        null = open(os.devnull, 'w')
        self.model = badread.qscore_model.QScoreModel('random', output=null)
        null.close()

    def check_min_max_positions(self, sequence, fragment):
        min_indices = set()
        max_indices = set()
        for _ in range(self.trials):
            qscores, _, _ = badread.qscore_model.get_qscores(sequence, fragment, self.model)
            self.assertEqual(len(qscores), len(sequence))
            qscores = [badread.qscore_model.qscore_char_to_val(q) for q in qscores]
            min_indices.add(qscores.index(min(qscores)))
            max_indices.add(qscores.index(max(qscores)))
        self.assertEqual(len(min_indices), len(sequence))
        self.assertEqual(len(max_indices), len(sequence))

    def test_get_qscores_1(self):
        self.check_min_max_positions('ACGACTACGTCAGACT',
                                     'ACGACTACGTCAGACT')

    def test_get_qscores_2(self):
        self.check_min_max_positions('ACGACTACCTCAGACT',
                                     'ACGACTACGTCAGACT')

    def test_get_qscores_3(self):
        self.check_min_max_positions('ACGACTACACT',
                                     'ACGACTACGTCAGACT')


class TestGetQScoresIdeal(unittest.TestCase):
    """
    For the ideal qscore model, the worst qscore should reliably show up at mismatch and insertion
    positions.
    """
    def setUp(self):
        self.trials = 100
        null = open(os.devnull, 'w')
        self.model = badread.qscore_model.QScoreModel('ideal', output=null)
        null.close()

    def check_min_positions(self, sequence, fragment, expected_min):
        min_indices = set()
        for _ in range(self.trials):
            qscores, _, _ = badread.qscore_model.get_qscores(sequence, fragment, self.model)
            self.assertEqual(len(qscores), len(sequence))
            qscores = [badread.qscore_model.qscore_char_to_val(q) for q in qscores]
            min_indices.add(qscores.index(min(qscores)))
        self.assertEqual(len(min_indices), 1)
        self.assertEqual(list(min_indices)[0], expected_min)

    def test_get_qscores_1(self):
        self.check_min_positions('ACGACTACCTCAGACT',
                                 'ACGACTACGTCAGACT', 8)

    def test_get_qscores_2(self):
        self.check_min_positions('ACGACTCACGTCAGACT',
                                 'ACGACTACGTCAGACT', 6)


class TestGetQScoresLoadedModel(unittest.TestCase):
    """
    For the loaded qscore model, the lowest qscore should usually show up at mismatch and insertion
    positions, but not always.
    """
    def setUp(self):
        self.trials = 1000
        null = open(os.devnull, 'w')
        model_filename = os.path.join(os.path.dirname(__file__), 'simple_qscore_model')
        self.model = badread.qscore_model.QScoreModel(model_filename, output=null)
        null.close()

    def check_min_positions(self, sequence, fragment, expected_min):
        min_position_counts = collections.defaultdict(int)
        for _ in range(self.trials):
            qscores, _, _ = badread.qscore_model.get_qscores(sequence, fragment, self.model)
            self.assertEqual(len(qscores), len(sequence))
            qscores = [badread.qscore_model.qscore_char_to_val(q) for q in qscores]
            min_position_counts[qscores.index(min(qscores))] += 1
        self.assertEqual(len(min_position_counts), len(sequence))
        self.assertEqual(max(min_position_counts, key=min_position_counts.get), expected_min)

    def test_get_qscores_1(self):
        self.check_min_positions('ACGACTACCTCAGACT',
                                 'ACGACTACGTCAGACT', 8)

    def test_get_qscores_2(self):
        self.check_min_positions('ACGACTCACGTCAGACT',
                                 'ACGACTACGTCAGACT', 6)


class TestMakeQScoreModel(unittest.TestCase):

    def setUp(self):
        self.null = open(os.devnull, 'w')
        self.ref_filename = os.path.join(os.path.dirname(__file__), 'test_alignment_ref.fasta')
        self.ref_filename_bad = os.path.join(os.path.dirname(__file__),
                                             'test_alignment_ref_bad_names.fasta')
        self.reads_filename = os.path.join(os.path.dirname(__file__), 'test_alignment_reads.fastq')
        self.reads_filename_bad = os.path.join(os.path.dirname(__file__),
                                               'test_alignment_reads_bad_names.fastq')
        self.paf_filename = os.path.join(os.path.dirname(__file__), 'test_alignment.paf')
        self.Args = collections.namedtuple('Args', ['reference', 'reads', 'alignment', 'k_size',
                                                    'max_alignments', 'max_del', 'min_occur',
                                                    'max_output'])

    def tearDown(self):
        self.null.close()

    def test_make_model_defaults(self):
        args = self.Args(reference=self.ref_filename, reads=self.reads_filename,
                         alignment=self.paf_filename, k_size=9, max_alignments=None, max_del=6,
                         min_occur=100, max_output=10000)
        with badread.misc.captured_output() as (out, err):
            badread.qscore_model.make_qscore_model(args, output=self.null, dot_interval=1)
        out = out.getvalue()
        out_lines = out.splitlines()
        self.assertEqual(len(out_lines), 6)
        self.assertTrue(out_lines[0].startswith('overall;'))
        self.assertTrue(out_lines[-1].startswith('=========;'))

    def test_make_model_k_size(self):
        args = self.Args(reference=self.ref_filename, reads=self.reads_filename,
                         alignment=self.paf_filename, k_size=5, max_alignments=None, max_del=6,
                         min_occur=100, max_output=10000)
        with badread.misc.captured_output() as (out, err):
            badread.qscore_model.make_qscore_model(args, output=self.null, dot_interval=1)
        out = out.getvalue()
        out_lines = out.splitlines()
        self.assertEqual(len(out_lines), 4)
        self.assertTrue(out_lines[0].startswith('overall;'))
        self.assertTrue(out_lines[-1].startswith('=====;'))

    def test_make_model_max_output(self):
        args = self.Args(reference=self.ref_filename, reads=self.reads_filename,
                         alignment=self.paf_filename, k_size=9, max_alignments=None, max_del=6,
                         min_occur=100, max_output=2)
        with badread.misc.captured_output() as (out, err):
            badread.qscore_model.make_qscore_model(args, output=self.null, dot_interval=1)
        out = out.getvalue()
        out_lines = out.splitlines()
        self.assertEqual(len(out_lines), 3)
        self.assertTrue(out_lines[0].startswith('overall;'))
        self.assertTrue(out_lines[1].startswith('=;'))
        self.assertTrue(out_lines[2].startswith('===;'))

    def test_make_model_bad_read_names(self):
        args = self.Args(reference=self.ref_filename, reads=self.reads_filename_bad,
                         alignment=self.paf_filename, k_size=9, max_alignments=None, max_del=6,
                         min_occur=100, max_output=2)
        with badread.misc.captured_output() as _:
            with self.assertRaises(SystemExit) as cm:
                badread.qscore_model.make_qscore_model(args, output=self.null, dot_interval=1)
        self.assertTrue('are you sure your read file and alignment file' in str(cm.exception))

    def test_make_model_bad_ref_names(self):
        args = self.Args(reference=self.ref_filename_bad, reads=self.reads_filename,
                         alignment=self.paf_filename, k_size=9, max_alignments=None, max_del=6,
                         min_occur=100, max_output=2)
        with badread.misc.captured_output() as _:
            with self.assertRaises(SystemExit) as cm:
                badread.qscore_model.make_qscore_model(args, output=self.null, dot_interval=1)
        self.assertTrue('are you sure your reference file and alignment file' in str(cm.exception))


class TestBugs(unittest.TestCase):
    """
    These tests are for bugs I found (and fixed).
    """
    def setUp(self):
        null = open(os.devnull, 'w')
        self.model = badread.qscore_model.QScoreModel('nanopore2018', output=null)
        null.close()

    def test_bug_1(self):
        """
        This case involves an alignment that ends in deletions. It caused a crash in the
        get_qscores function.
        """
        seq = 'CGGGCGCAACGCGTTCGATGCTCCACGTCAGTGAGCCTAAGCATATAAGCGAAAGGCT'
        frag = 'CGTCCGCTACGGCGGCAGTTCCCCATTCTTCCCCCGCATCGAGTGATAAACCGTAAACATGGGCGTAGACGGCATCCCCT'
        qscores, _, _ = badread.qscore_model.get_qscores(seq, frag, self.model)
        self.assertEqual(len(seq), len(qscores))

    def test_bug_2(self):
        """
        This case involves an alignment that starts with deletions. It caused a crash in the
        get_qscores function.
        """
        seq = 'CGGCGGCAGTTCCCCATTCTTCCCCCGCATCGAGTGATAAACCGTAAACATGGGCGTAGACGGCATCCCCT'
        frag = 'ATATCGGCGGCAGTTCCCCATTCTTCCCCCGCATCGAGTGATAAACCGTAAACATGGGCGTAGACGGCATCCCCT'
        qscores, _, _ = badread.qscore_model.get_qscores(seq, frag, self.model)
        self.assertEqual(len(seq), len(qscores))
