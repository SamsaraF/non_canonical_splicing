import unittest
import numpy as np
import pandas as pd


class TestSplicingExtract(unittest.TestCase):
    def test_splicing_event_count(self):
        # 20 mock cases
        for case in range(1, 21):
            splicing_site_defs = f'./data/test_data/sj_info.test_case_{case}.tsv'
            expected_output = f'./results/test_case_{case}.splicing_site.tsv'

            # ignore strand compare because NC events strand is None
            in_df = pd.read_table(splicing_site_defs, header=None, usecols=[0, 1, 2, 6], index_col=[0, 1, 2], names=['chr', 'start', 'end', 'count'])
            exp_df = pd.read_table(expected_output, header=None, index_col=[0, 1, 2], names=['chr', 'start', 'end', 'count'])

            self.assertEqual(in_df.shape, exp_df.shape)
            self.assertEqual(sorted(in_df.index.tolist()), sorted(exp_df.index.tolist()))

            all_df = pd.merge(in_df, exp_df, left_index=True, right_index=True, suffixes=('_in', '_exp'))
            all_df['diff'] = (all_df['count_in'] - all_df['count_exp']) / all_df['count_in']

            # allow count difference less than 5%
            self.assertTrue(np.all(all_df['diff'].abs() < 0.05), f'Case {case} has more than 5% difference in splicing site counts:\n{all_df}')


    def test_splicing_event_type(self):
        # 20 mock cases
        for case in range(1, 21):
            splicing_site_defs = f'./data/test_data/sj_info.test_case_{case}.tsv'
            expected_output = f'./results/test_case_{case}.splicing_site.tsv'

            in_df = pd.read_table(splicing_site_defs, header=None, usecols=[0, 1, 2, 3], index_col=[0, 1, 2], names=['chr', 'start', 'end', 'type'])
            exp_df = pd.read_table(expected_output, header=None, index_col=[0, 1, 2], names=['chr', 'start', 'end', 'type'])

            self.assertEqual(in_df.shape, exp_df.shape)
            self.assertEqual(sorted(in_df.index.tolist()), sorted(exp_df.index.tolist()))

            all_df = pd.merge(in_df, exp_df, left_index=True, right_index=True, suffixes=('_in', '_exp'))
            all_df['type_match'] = all_df['type_in'] == all_df['type_exp']

            self.assertTrue(all_df['type_match'].all(), f'Case {case} has mismatched splicing event types:\n{all_df[~all_df["type_match"]]}')



if __name__ == '__main__':
    unittest.main()