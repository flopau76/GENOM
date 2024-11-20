from kmers import encode_nucl, encode_nucl_rev, encode_kmer, stream_kmers
import unittest

class TestEncode(unittest.TestCase):
    def test_encode_nucl(self):
        expected = ['A', 'C', 'T', 'G']
        for i in range(4):
            self.assertEqual(encode_nucl(expected[i]), i)

    def test_encode_nucl_rev(self):
        expected = ['T', 'G', 'A', 'C']
        for i in range(4):
            self.assertEqual(encode_nucl_rev(expected[i]), i)

    def test_encode_kmer(self):
        seq = 'ACT'
        k = 3
        expected = (0b000110, 0b001110) # (kmer, rev_kmer)
        self.assertEqual(encode_kmer(seq, k), expected)

    def test_stream_kmers(self):
        seq = 'AGCTA'
        k = 3
        expected = [min(encode_kmer('AGC', 3)), min(encode_kmer('GCT', 3)), min(encode_kmer('CTA', 3))]
        self.assertEqual(list(stream_kmers([seq], k)), expected)

if __name__ == "__main__":
    unittest.main()