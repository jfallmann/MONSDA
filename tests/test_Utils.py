import unittest

from MONSDA.Utils import (
    NestedDefaultDict,
    add_to_innermost_key_by_list,
    check_ref,
    comment_remover,
    convertcol,
    depth,
    dict_inst,
    find_all_values_on_key,
    find_innermost_value_from_dict,
    find_key_for_value,
    get_dict_hash,
    get_from_dict,
    gethighest_dict,
    gethighest_list,
    getlowest_dict,
    getlowest_list,
    idfromfa,
    isinvalid,
    isvalid,
    keys_from_dict,
    keysets_from_dict,
    list_all_keys_of_dict,
    list_all_values_of_dict,
    makelogdir,
    makeoutdir,
    merge_dicts,
    multi_replace,
    nested_set,
    npprint,
    parseseq,
    removekey,
    rmempty,
    sub_dict,
    subset_dict,
    value_extract,
    yield_from_dict,
)


class TestUtils(unittest.TestCase):

    def test_NestedDefaultDict(self):
        nested_dict = NestedDefaultDict(lambda: NestedDefaultDict(int))
        nested_dict['a']['b'] = 1
        self.assertEqual(nested_dict["a"]["b"], 1)
        self.assertEqual(nested_dict['c']['d'], 0)

    def test_rmempty(self):
        files = ["file1.txt", "file2.txt"]
        self.assertEqual(rmempty(files), files)

    def test_comment_remover(self):
        input_text = [
            "code line 1",
            "// this is a comment",
            "code line 2 /* inline comment */",
        ]
        expected_output = ["code line 1", "code line 2 "]
        self.assertEqual(comment_remover(input_text), expected_output)

    def test_dict_inst(self):
        self.assertTrue(dict_inst({'key': 'value'}))
        self.assertFalse(dict_inst(['not', 'a', 'dict']))

    def test_get_from_dict(self):
        data_dict = {'a': {'b': {'c': 'd'}}}
        map_list = ['a', 'b', 'c']
        self.assertEqual(get_from_dict(data_dict, map_list), ["d"])

    def test_yield_from_dict(self):
        data_dict = {'a': 1, 'b': 2, 'c': {'a': 3, 'b': 4}}
        key = 'a'
        expected_output = [1, 3]
        self.assertEqual(list(yield_from_dict(key, data_dict)), expected_output)

    def test_sub_dict(self):
        data_dict = {"a": {"b": {"c": "d"}}}
        map_list = ["a", "b"]
        self.assertEqual(sub_dict(data_dict, map_list), {"c": "d"})

    def test_subset_dict(self):
        data_dict = {"a": {"b": {"c": "d"}}}
        map_list = ["a", "b"]
        self.assertEqual(subset_dict(data_dict, map_list), {"a": {"b": {"c": "d"}}})

    def test_nested_set(self):
        dic = {}
        keys = ["a", "b", "c"]
        value = "d"
        nested_set(dic, keys, value)
        self.assertEqual(dic, {"a": {"b": {"c": "d"}}})

    def test_merge_dicts(self):
        d1 = {"a": 1, "b": {"c": 2}}
        d2 = {"b": {"d": 3}, "e": 4}
        self.assertEqual(merge_dicts(d1, d2), {"a": 1, "b": {"c": 2, "d": 3}, "e": 4})

    def test_keysets_from_dict(self):
        data_dict = {"a": {"b": {"c": "d"}}}
        self.assertEqual(keysets_from_dict(data_dict), [("a", "b", "c")])

    def test_keys_from_dict(self):
        data_dict = {"a": {"b": {"c": "d"}}}
        self.assertEqual(keys_from_dict(data_dict), {0: ["a"], 1: ["b"], 2: ["c"]})

    def test_depth(self):
        data_dict = {"a": {"b": {"c": "d"}}}
        self.assertEqual(depth(data_dict), 3)

    def test_list_all_keys_of_dict(self):
        data_dict = {"a": {"b": {"c": "d"}}}
        self.assertEqual(list(list_all_keys_of_dict(data_dict)), ["a", "b", "c"])

    def test_list_all_values_of_dict(self):
        data_dict = {"a": {"b": {"c": "d"}}}
        self.assertEqual(list(list_all_values_of_dict(data_dict)), [("c", "d")])

    def test_find_all_values_on_key(self):
        data_dict = {"a": 1, "b": 2, "c": {"a": 3, "b": 4}}
        key = "a"
        expected_output = [1, 3]
        self.assertEqual(list(find_all_values_on_key(key, data_dict)), expected_output)

    def test_find_key_for_value(self):
        data_dict = {"a": 1, "b": 2, "c": {"a": 3, "b": 4}}
        val = 3
        self.assertEqual(list(find_key_for_value(val, data_dict)), ["a"])

    def test_value_extract(self):
        data_dict = {"a": 1, "b": 2, "c": {"a": 3, "b": 4}}
        key = "a"
        expected_output = [1, 3]
        self.assertEqual(list(value_extract(key, data_dict)), expected_output)

    def test_find_innermost_value_from_dict(self):
        data_dict = {"a": {"b": {"c": "d"}}}
        self.assertEqual(find_innermost_value_from_dict(data_dict), "d")

    def test_removekey(self):
        data_dict = {"a": 1, "b": 2}
        self.assertEqual(removekey(data_dict, "a"), {"b": 2})

    def test_getlowest_list(self):
        a = [5, 1, 3, 2, 4]
        n = 3
        self.assertEqual(getlowest_list(a, n), [1, 2, 3])

    def test_gethighest_list(self):
        a = [5, 1, 3, 2, 4]
        n = 3
        self.assertEqual(gethighest_list(a, n), [3, 4, 5])

    def test_getlowest_dict(self):
        a = {"a": 5, "b": 1, "c": 3, "d": 2, "e": 4}
        n = 3
        self.assertEqual(getlowest_dict(a, n), {"b": 1, "d": 2, "c": 3})

    def test_gethighest_dict(self):
        a = {"a": 5, "b": 1, "c": 3, "d": 2, "e": 4}
        n = 3
        self.assertEqual(gethighest_dict(a, n), {"a": 5, "e": 4, "c": 3})

    def test_toarray(self):
        # Assuming toarray function converts a file to np.array
        # This test might need to create a temporary file to fully test toarray functionality
        pass

    def test_convertcol(self):
        self.assertEqual(convertcol("1.23"), 1.23)
        self.assertTrue(np.isnan(convertcol("nan")))

    def test_isvalid(self):
        self.assertTrue(isvalid("1.23"))
        self.assertFalse(isvalid("nan"))

    def test_isinvalid(self):
        self.assertTrue(isinvalid("nan"))
        self.assertFalse(isinvalid("1.23"))

    def test_makeoutdir(self):
        outdir = "test_dir"
        self.assertEqual(makeoutdir(outdir), os.path.abspath(outdir))
        shutil.rmtree(outdir)

    def test_parseseq(self):
        sequence = "ATCG"
        self.assertEqual(parseseq(sequence).read(), ">Seq1:default:nochrom:(.)\nATCG")

    def test_npprint(self):
        a = np.array([1.23, 4.56])
        expected_output = "1\t1.2300000\n2\t4.5600000\n"
        with self.assertLogs(level="INFO") as log:
            npprint(a)
            self.assertIn(expected_output.strip(), log.output[0])

    def test_idfromfa(self):
        id = "cluster1:chr19.tRNA5-LysCTT(+)"
        self.assertEqual(idfromfa(id), ["cluster1", "chr19", "+"])

    def test_cluster2trna(self):
        # Assuming cluster2trna function converts cluster to tRNA
        # This test might need to create a temporary file to fully test cluster2trna functionality
        pass

    def test_check_ref(self):
        reference = "test_ref"
        with open(reference, "w") as f:
            f.write("test")
        self.assertEqual(check_ref(reference), reference)
        os.remove(reference)

    def test_multi_replace(self):
        repl = {"a": "1", "b": "2"}
        text = "a and b"
        self.assertEqual(multi_replace(repl, text), "1 and 2")

    def test_makelogdir(self):
        logdir = "test_logdir"
        self.assertEqual(makelogdir(logdir), os.path.abspath(logdir))
        shutil.rmtree(logdir)

    def test_get_dict_hash(self):
        d = {"a": 1, "b": 2}
        self.assertEqual(
            get_dict_hash(d),
            hashlib.sha256(bytes(str(sorted(d.items())), "utf-8")).hexdigest(),
        )

    def test_add_to_innermost_key_by_list(self):
        addto = {}
        toadd = "value"
        keylist = ["a", "b", "c"]
        self.assertEqual(
            add_to_innermost_key_by_list(addto, toadd, keylist),
            {"a": {"b": {"c": "value"}}},
        )


if __name__ == '__main__':
    unittest.main()
