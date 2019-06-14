import numpy
from collections import OrderedDict
from itertools import repeat

from ArgumentConfig.AnalysisType import AnalysisType
from Utilities.SafeCastUtil import SafeCastUtil


class GeneListComboUtility(object):

    @staticmethod
    def determineGeneListCombos(gene_lists, feature_names):
        gene_sets_across_files = {}
        for feature in feature_names:
            split = feature.split(".")
            if gene_sets_across_files.get(split[0]) is not None:
                gene_sets_across_files[split[0]].append(feature)
            else:
                gene_sets_across_files[split[0]] = [feature]

        numerical_permutations = GeneListComboUtility.generateNumericalPermutations(gene_lists, gene_sets_across_files)
        gene_list_keys = SafeCastUtil.safeCast(gene_lists.keys(), list)
        file_keys = SafeCastUtil.safeCast(gene_sets_across_files.keys(), list)
        gene_list_combos = []
        for perm in numerical_permutations:
            feature_strings = []
            for i in range(0, len(perm)):
                file_name = file_keys[i]
                gene_list = gene_lists[gene_list_keys[SafeCastUtil.safeCast(perm[i], int)]]
                if len(gene_list) > 0:
                    feature_strings.append([file_name + "." + gene for gene in gene_list if len(gene.strip()) > 0])
            if len(feature_strings) > 0:
                gene_list_combos.append(feature_strings)

        file_keys = SafeCastUtil.safeCast(gene_sets_across_files.keys(), list)
        gene_list_keys = SafeCastUtil.safeCast(gene_lists.keys(), list)
        expected_combo_length = (len(gene_list_keys) ** len(file_keys)) - 1
        return gene_list_combos, expected_combo_length

    @staticmethod
    def generateNumericalPermutations(gene_lists, gene_sets_across_files):
        max_depth = len(gene_lists) - 1
        num_files = len(gene_sets_across_files)
        all_arrays = []
        current_array = SafeCastUtil.safeCast(numpy.zeros(num_files, dtype=numpy.int), list)
        target_index = num_files - 1
        while target_index >= 0:
            if current_array not in all_arrays:
                clone_array = current_array[:]
                all_arrays.append(clone_array)
            if current_array[target_index] < max_depth:
                current_array[target_index] += 1
                while len(current_array) > target_index + 1 and current_array[target_index + 1] < max_depth:
                    target_index += 1
            else:
                target_index -= 1
                for subsequent_index in range(target_index, len(current_array) - 1):
                    current_array[subsequent_index + 1] = 0
        return all_arrays

    @staticmethod
    def generateFeatureSetString(feature_set, gene_lists, combine_gene_lists, analysis_type):
        feature_map = {}
        for feature_list in feature_set:
            for feature in feature_list:
                file_name = feature.split(".")[0]
                feature_name = feature.split(".")[1:][0]
                if feature_map.get(file_name):
                    feature_map[file_name].append(feature_name)
                else:
                    feature_map[file_name] = [feature_name]

        feature_set_string = ""
        num_gene_lists_deduped = len(GeneListComboUtility.fetchAllGeneListGenesDeduped(gene_lists))
        for file_key in feature_map.keys():
            if combine_gene_lists and len(feature_map[file_key]) == num_gene_lists_deduped:
                feature_set_string += (file_key + ":ALL_GENE_LISTS ")
            else:
                for gene_list_key in gene_lists.keys():
                    if len(feature_map[file_key]) == len(gene_lists[gene_list_key]):
                        feature_map[file_key].sort()
                        gene_lists[gene_list_key].sort()
                        same_list = True
                        for i in range(0, len(gene_lists[gene_list_key])):
                            if gene_lists[gene_list_key][i] != feature_map[file_key][i]:
                                same_list = False
                        if same_list:
                            feature_set_string += (file_key + ":" + gene_list_key + " ")
        if feature_set_string == "" and analysis_type is AnalysisType.SPEARMAN_NO_GENE_LISTS:
            return "all_features"  # TODO: This is a bit lazy, do it smarter.
        return feature_set_string.strip()

    @staticmethod
    def fetchAllGeneListGenesDeduped(gene_lists):
        all_genes = SafeCastUtil.safeCast(gene_lists.values(), list)
        concated_genes = SafeCastUtil.safeCast(numpy.concatenate(all_genes), list)
        dedupded_genes = list(OrderedDict(zip(concated_genes, repeat(None))))
        return dedupded_genes
