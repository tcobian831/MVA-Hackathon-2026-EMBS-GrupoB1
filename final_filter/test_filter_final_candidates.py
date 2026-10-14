"""Focused regression tests: strict limits, missing AF and allele matching."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import filter_final_candidates as app


class FilterTests(unittest.TestCase):
    def test_thresholds_missing_and_protected_control(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            fields = ['Allele', 'Consequence', 'BIOTYPE', 'SYMBOL', 'CANONICAL', 'MANE'] + ['gnomADg_' + g + '_AF' for g in app.GROUPS]
            header = '##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: ' + '|'.join(fields) + '">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\n'
            def row(pos, score, afs, window=False, chrom='chr1', ref='A', alt='G', original=None, extra=None):
                ann = '|'.join([alt, 'missense_variant', 'protein_coding', 'TEST', 'YES', ''] + afs)
                info = 'CSQ=' + ann + (',' + extra if extra else '')
                if score is not None: info += ';CADD_PHRED=' + str(score)
                if window: info += ';F1_REASON=MVA_WINDOW'
                if original: info += ';ORIGINAL=' + original
                return '\t'.join([chrom, str(pos), '.', ref, alt, '.', 'PASS', info, 'GT:DP:GQ', '0/1:30:99']) + '\n'
            empty = [''] * 6
            # 1/3 exactly on CADD cutoffs fail; 5 exactly on AF fails.
            data = [row(1, 23, empty), row(2, 23.1, empty), row(3, 18, empty, True),
                    row(4, 18.1, empty, True), row(5, 30, ['.0002'] + [''] * 5),
                    row(6, 30, ['.0003'] + [''] * 5), row(7, 30, ['.0001'] + [''] * 5),
                    row(8, None, empty), row(9, None, ['.001'] + [''] * 5),
                    row(40209701, None, ['.1'] * 6, chrom='chr15', ref='T'),
                    # Wrong allele C has a common AF; own allele G must pass.
                    row(11, 30, ['.0001'] * 6, original='chr1|11|A|C,G|2',
                        extra='|'.join(['C','missense_variant','protein_coding','TEST','YES',''] + ['.1'] * 6))]
            inp = d / 'input.vcf'; inp.write_text(header + ''.join(data))
            with patch('sys.argv', ['filter', str(inp), '--out-dir', str(d / 'out')]), contextlib.redirect_stdout(io.StringIO()):
                app.main()
            s = json.loads((d / 'out/summary.json').read_text())
            self.assertEqual(s['selected'], 5)
            self.assertEqual(s['pending'], 1)
            self.assertEqual(s['status_counts']['rejected_cadd'], 2)
            self.assertEqual(s['status_counts']['rejected_frequency'], 3)
            self.assertTrue(s['control_found'])
            # Existing output must never be overwritten.
            with patch('sys.argv', ['filter', str(inp), '--out-dir', str(d / 'out')]), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit): app.main()


if __name__ == '__main__':
    unittest.main()
