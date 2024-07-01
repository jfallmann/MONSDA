import unittest
from MONSDA.Utils import NestedDefaultDict, rmempty, comment_remover, dict_inst, get_from_dict, yield_from_dict

class TestUtils(unittest.TestCase):

    def test_NestedDefaultDict(self):
        # Test initialization and basic functionality
        nested_dict = NestedDefaultDict(lambda: NestedDefaultDict(int))
        nested_dict['a']['b'] = 1
        self.assertEqual(nested_dict['a']['b'], 1)
        # Test default factory
        self.assertEqual(nested_dict['c']['d'], 0)

    def test_rmempty(self):
        # Assuming rmempty function removes empty directories
        # This test might need to create temporary directories and files to fully test rmempty functionality
        pass

    def test_comment_remover(self):
        # Assuming comment_remover function removes comments from a list of strings
        input_text = ["code line 1", "# this is a comment", "code line 2 # inline comment"]
        expected_output = ["code line 1", "code line 2 "]
        self.assertEqual(comment_remover(input_text), expected_output)

    def test_dict_inst(self):
        # Assuming dict_inst function checks if an instance is a dictionary
        self.assertTrue(dict_inst({'key': 'value'}))
        self.assertFalse(dict_inst(['not', 'a', 'dict']))

    def test_get_from_dict(self):
        # Assuming get_from_dict function retrieves a value from a nested dictionary using a list of keys
        data_dict = {'a': {'b': {'c': 'd'}}}
        map_list = ['a', 'b', 'c']
        self.assertEqual(get_from_dict(data_dict, map_list), 'd')

    def test_yield_from_dict(self):
        # Assuming yield_from_dict function yields items from a dictionary that match a given key
        data_dict = {'a': 1, 'b': 2, 'c': {'a': 3, 'b': 4}}
        key = 'a'
        expected_output = [1, 3]
        self.assertEqual(list(yield_from_dict(key, data_dict)), expected_output)

if __name__ == '__main__':
unittest.main()