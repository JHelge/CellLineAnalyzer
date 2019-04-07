import logging
import multiprocessing
import os
import re
import mmap

from ArgumentConfig.IndividualTrainConfig import IndividualTrainConfig
from ArgumentConfig.ProcessedArguments import ProcessedArguments
from ArgumentConfig.RSENConfig import RSENConfig
from SupportedMachineLearningAlgorithms import SupportedMachineLearningAlgorithms
from Utilities.PercentageBarUtility import PercentageBarUtility
from Utilities.SafeCastUtil import SafeCastUtil


class ArgumentProcessingService(object):

    log = logging.getLogger(__name__)
    logging.basicConfig()
    log.setLevel(logging.INFO)

    ARGUMENTS_FILE = "arguments.txt"
    GENE_LISTS = "gene_list"

    UNFILLED_VALUE_PLACEHOLDER = "'0'"

    RESULTS = "results"
    IS_CLASSIFIER = "is_classifier"
    FEATURES = "features"
    FEATURE_NAMES = "featureNames"
    INNER_MONTE_CARLO_PERMUTATIONS = "inner_monte_carlo_permutations"
    OUTER_MONTE_CARLO_PERMUTATIONS = "outer_monte_carlo_permutations"
    DATA_SPLIT = "data_split"
    NUM_THREADS = "num_threads"
    ALGORITHM_CONFIGS = "algorithm_configs"
    RECORD_DIAGNOSTICS = "record_diagnostics"

    # RSEN Specific Arguments
    RSEN_P_VAL = "rsen_p_val"
    RSEN_K_VAL = "rsen_k_val"
    RSEN_COMBINE_GENE_LISTS = "rsen_combine_gene_lists"
    BINARY_CATEGORICAL_MATRIX = "binary_categorical_matrix"

    # For AnalysisType.FULL_CLA_SPECIFIC_COMBO
    SPECIFIC_COMBOS = "specific_combos"

    # For AnalysisType.SPEARMAN_NO_GENE_LISTS
    SPEARMAN_CORR = "spearman_corr"

    # For AnalysisType.INDIVIDUAL_TRAIN
    INDIVIDUAL_TRAIN_ALGORITHM = "individual_train_algorithm"
    INDIVIDUAL_TRAIN_HYPERPARAMS = "individual_train_hyperparams"
    INDIVIDUAL_TRAIN_FEATURE_GENE_LIST_COMBO = "individual_train_combo"

    def __init__(self, input_folder):
        self.input_folder = input_folder

    def handleInputFolder(self):
        directory_contents = os.listdir(self.input_folder)
        if not self.validateDirectoryContents(directory_contents):
            self.log.error("Invalid directory contents, needs a %s file.", self.ARGUMENTS_FILE)
            return None

        arguments = self.fetchArguments(self.input_folder + "/" + self.ARGUMENTS_FILE)
        results_file = arguments.get(self.RESULTS)
        is_classifier = SafeCastUtil.safeCast(arguments.get(self.IS_CLASSIFIER), int) == 1
        analyze_all = self.fetchOrReturnDefault(arguments.get(self.SPEARMAN_CORR), bool, False)

        if is_classifier is None or results_file is None:
            return None
        results_list = self.validateAndExtractResults(results_file, is_classifier)

        gene_lists = self.extractGeneLists()

        write_diagnostics = self.fetchOrReturnDefault(arguments.get(self.RECORD_DIAGNOSTICS), bool, False)
        feature_files = [file for file in os.listdir(self.input_folder) if self.fileIsFeatureFile(file, results_file)]
        if analyze_all:
            feature_map = self.createAndValidateFullFeatureMatrix(results_list, feature_files)
        else:
            feature_map = self.createAndValidateFeatureMatrix(results_list, gene_lists, write_diagnostics, feature_files)
        binary_cat_matrix = self.fetchBinaryCatMatrixIfApplicable(arguments, gene_lists, results_list, analyze_all)

        algorithm_configs = self.handleAlgorithmConfigs(arguments)
        if not feature_map or not results_list:
            return None
        inner_monte_carlo_perms = self.fetchOrReturnDefault(arguments.get(self.INNER_MONTE_CARLO_PERMUTATIONS), int, 10)
        outer_monte_carlo_perms = self.fetchOrReturnDefault(arguments.get(self.OUTER_MONTE_CARLO_PERMUTATIONS), int, 10)
        data_split = self.fetchOrReturnDefault(arguments.get(self.DATA_SPLIT), float, 0.8)
        num_threads = self.fetchOrReturnDefault(arguments.get(self.NUM_THREADS), int, multiprocessing.cpu_count())

        individual_train_config = self.createIndividualTrainConfig(arguments)
        rsen_config = self.createRSENConfig(arguments, binary_cat_matrix)
        specific_combos = self.determineSpecificCombos(arguments.get(self.SPECIFIC_COMBOS))

        return ProcessedArguments(results_list, is_classifier, feature_map, gene_lists, inner_monte_carlo_perms,
                                  outer_monte_carlo_perms, data_split, algorithm_configs, num_threads,
                                  write_diagnostics, individual_train_config, rsen_config, specific_combos, analyze_all)

    def validateDirectoryContents(self, directory_contents):
        return self.ARGUMENTS_FILE in directory_contents

    def fetchArguments(self, arguments_file):
        arguments = {}
        with open(arguments_file) as data_file:
            try:
                for line in data_file:
                    line_trimmed_split = line.strip().split("=")
                    if len(line_trimmed_split) > 1:
                        arguments[line_trimmed_split[0]] = line_trimmed_split[1]
            except ValueError as valueError:
                self.log.error(valueError)
            finally:
                self.log.debug("Closing file %s", arguments_file)
                data_file.close()
        return arguments

    def validateAndExtractResults(self, results_file, is_classifier):
        sample_list = []
        cast_type = float
        if is_classifier:
            cast_type = int
        results_path = self.input_folder + "/" + results_file
        with open(results_path) as data_file:
            try:
                for line_index, line in enumerate(data_file):
                    if len(re.findall(r'^\s*$', line)) > 0 or line_index == 0:  # header or whitespace
                        continue
                    line_trimmed_split = line.strip().split(",")
                    if len(line_trimmed_split) != 2:
                        self.log.error("Each line in %s must be 2 columns. Aborting argument processing.",
                                       results_file)
                        raise ValueError("Each line in results file must be 2 columns.")
                    cell_line = line_trimmed_split[0]
                    cell_result = SafeCastUtil.safeCast(line_trimmed_split[1], cast_type)
                    if cell_line in sample_list:
                        self.log.error("Repeated cell line name: %s. Aborting argument processing.", cell_line)
                        raise ValueError("Repeated cell line name.")
                    else:
                        sample_list.append([cell_line, cell_result])
            except ValueError as valueError:
                self.log.error(valueError)
            finally:
                self.log.debug("Closing file %s", results_file)
                data_file.close()
        return sample_list

    def extractGeneLists(self):
        gene_lists = {"null_gene_list": []}
        files = os.listdir(self.input_folder)
        for file in [f for f in files if self.GENE_LISTS in f]:
            file_path = self.input_folder + "/" + file
            with open(file_path) as gene_list_file:
                genes = gene_list_file.read().strip().split(",")
                genes_deduped = []
                [genes_deduped.append(g.strip()) for g in genes if g not in genes_deduped and len(g.strip()) > 0]
                if len(genes_deduped) > 0:
                    gene_lists[file.split(".csv")[0]] = genes_deduped
                else:
                    self.log.warning("No genes found in gene list %s, will not process.", file)

        return gene_lists

    def createAndValidateFullFeatureMatrix(self, results_list, feature_files):
        feature_matrix = {self.FEATURE_NAMES: []}
        for cell_line in results_list:
            feature_matrix[cell_line[0]] = []
        for file in feature_files:
            self.log.info("Fetching all features for file %s", file)
            file_name = file.split(".")[0]
            features_path = self.input_folder + "/" + file
            num_lines = self.countLines(features_path)
            with open(features_path) as feature_file:
                num_features = 0
                valid_indices = []
                try:
                    for line_index, line in enumerate(feature_file):
                        line_split = line.split(",")
                        percent_done, percentage_bar = PercentageBarUtility.calculateAndCreatePercentageBar(line_index,
                                                                                                            num_lines)
                        self.log.info("Total progress extracting rows for %s: %s%% done:\n %s",
                                      file, percent_done, percentage_bar)
                        if line_index == 0:
                            valid_indices, feature_names = self.fetchUniqueFeatureNamesAndIndices(line_split, file_name)
                            feature_matrix[self.FEATURE_NAMES] += feature_names
                            num_features = len(line_split)
                            continue
                        for i in range(0, num_features):
                            if i not in valid_indices:
                                continue
                            feature = line_split[i] if line_split[i] is not None else self.UNFILLED_VALUE_PLACEHOLDER
                            feature_as_float = SafeCastUtil.safeCast(feature, float)
                            if feature_as_float is not None:
                                feature_matrix[results_list[line_index - 1][0]].append(feature_as_float)
                            else:
                                feature_matrix[results_list[line_index - 1][0]].append(feature.strip())

                except ValueError as valueError:
                    self.log.error(valueError)
                finally:
                    self.log.debug("Closing file %s", feature_file)
                    feature_file.close()
        return feature_matrix

    def countLines(self, file):
        open_file = open(file, "r+")
        buf = mmap.mmap(open_file.fileno(), 0)
        lines = 0
        readline = buf.readline
        while readline():
            lines += 1
        return lines

    def fetchUniqueFeatureNamesAndIndices(self, line_split, file_name):
        unvalidated_features = [file_name + "." + name.strip() for name in line_split if len(name.strip()) > 0]
        valid_indices = []
        valid_features = []
        for i in range(0, len(unvalidated_features)):
            if unvalidated_features.count(unvalidated_features[i]) == 1:
                valid_indices.append(i)

        for i in range(0, len(unvalidated_features)):
            if i in valid_indices:
                valid_features.append(unvalidated_features[i])
        return valid_indices, valid_features

    def createAndValidateFeatureMatrix(self, results_list, gene_lists, write_diagnostics, feature_files):
        incomplete_features = []
        for file in feature_files:
            features_path = self.input_folder + "/" + file
            validated_features, num_features = self.validateGeneLists(features_path, file, gene_lists)
            incomplete_features.append([file, validated_features, num_features])

        if write_diagnostics:
            self.writeDiagnostics(incomplete_features)

        feature_matrix = {self.FEATURE_NAMES: []}
        for file in feature_files:
            features_path = self.input_folder + "/" + file
            self.extractFeatureMatrix(feature_matrix, features_path, file, gene_lists, results_list)
        return feature_matrix

    def validateGeneLists(self, features_path, file, gene_lists):
        features_missing_from_files = {}
        num_features = 0
        with open(features_path) as feature_file:
            try:
                for line_index, line in enumerate(feature_file):
                    if line_index == 0:
                        feature_names = line.split(",")
                        num_features = len(feature_names)
                        features_missing_from_files = self.validateAndTrimGeneList(feature_names, gene_lists, file)
                    break
            except ValueError as valueError:
                self.log.error(valueError)
                return features_missing_from_files
            finally:
                self.log.debug("Closing file %s", feature_file)
                feature_file.close()
        return features_missing_from_files, num_features

    def validateAndTrimGeneList(self, feature_list, gene_lists, file):
        unused_features = {}
        for key in gene_lists.keys():
            for gene in gene_lists[key]:
                if gene not in [feature.strip() for feature in feature_list]:
                    index = gene_lists[key].index(gene)
                    if unused_features.get(key) is None:
                        unused_features[key] = [[gene, index]]
                    else:
                        unused_features[key].append([gene, (index + len(unused_features[key]))])
                    self.log.warning("Incomplete dataset: gene %s from gene list %s not found in file %s. "
                                     "Will not process this gene in this file.", gene, key, file)
        return unused_features

    def writeDiagnostics(self, features_removed):
        with open(self.input_folder + "/Diagnostics.txt", "w") as diagnostics_file:
            try:
                diagnostics_file.write("### Diagnostics ###\n")
                for feature_file in features_removed:
                    diagnostics_file.write("\nFeatures from gene list(s) not available in " + feature_file[0] + ":\n")
                    for gene_list in feature_file[1].keys():
                        num_genes_missing = len(feature_file[1][gene_list])
                        percent_genes_missing = round((num_genes_missing / feature_file[2]) * 100, 2)
                        diagnostics_file.write("\t" + SafeCastUtil.safeCast(num_genes_missing, str) + " (" +
                                                      SafeCastUtil.safeCast(percent_genes_missing, str) + " %" + ")"
                                                    " features not present in " + gene_list + ".csv:\n")
                        for gene in feature_file[1][gene_list]:
                            diagnostics_file.write("\t\t" + gene[0] + " at index "
                                                   + SafeCastUtil.safeCast(gene[1], str) + "\n")
                diagnostics_file.write("\n\n######################\n\n")
            except ValueError as error:
                self.log.error("Error writing to file %s. %s", diagnostics_file, error)
            finally:
                diagnostics_file.close()

    def extractFeatureMatrix(self, feature_matrix, features_path, file, gene_lists, results_list):
        self.log.info("Extracting important features for %s.", file)
        gene_list_features = []
        for gene_list in gene_lists.values():
            for gene_list_feature in gene_list:
                if gene_list_feature not in gene_list_features:
                    gene_list_features.append(gene_list_feature)

        with open(features_path) as feature_file:
            try:
                important_feature_indices = []
                for line_index, line in enumerate(feature_file):
                    if line_index == 0:
                        feature_names = line.split(",")
                        for gene_list_feature in gene_list_features:
                            important_index = None
                            feature_name = self.determineFeatureName(gene_list_feature, file)
                            for i in range(0, len(feature_names)):
                                if feature_names[i].strip() == gene_list_feature.strip():
                                    important_index = i
                            if feature_name not in feature_matrix[self.FEATURE_NAMES]:
                                feature_matrix[self.FEATURE_NAMES].append(feature_name)
                            important_feature_indices.append(important_index)
                    else:
                        features = self.extractCastedFeatures(line, important_feature_indices)
                        cell_line = results_list[line_index - 1]
                        if not cell_line[0] in feature_matrix:
                            feature_matrix[cell_line[0]] = features
                        else:
                            feature_matrix[cell_line[0]] = feature_matrix[cell_line[0]] + features
                        if line_index > len(results_list):
                            self.log.error("Invalid line count for %s", file)
                            raise ValueError("Invalid line count for" + file + ". Must be " +
                                             SafeCastUtil.safeCast(file, str) + "lines long.")
            except ValueError as valueError:
                self.log.error("Please verify results file is the same number of rows as all feature files.")
                self.log.error(valueError)
                return None
            finally:
                feature_file.close()
                self.log.debug("Closing file %s", feature_file)

    def fileIsFeatureFile(self, file, results_file):
        algorithm_files = [algo + ".csv" for algo in SupportedMachineLearningAlgorithms.fetchAlgorithms()]

        return file != results_file and file != self.ARGUMENTS_FILE and self.GENE_LISTS not in file and\
               file not in algorithm_files and ".csv" in file.lower()

    def determineFeatureName(self, feature_name, file):
        return SafeCastUtil.safeCast(file.split(".")[0] + "." + feature_name.strip(), str)

    def extractCastedFeatures(self, line, important_feature_indices):
        important_features = []
        feature_values = line.strip().split(",")
        for index in important_feature_indices:
            if index is None:
                # TODO: Verify that this is acceptable, it works for one hot encoding and should never vary in any model
                important_features.append(self.UNFILLED_VALUE_PLACEHOLDER)
            else:
                if SafeCastUtil.safeCast(feature_values[index], float) is not None:
                    important_features.append(SafeCastUtil.safeCast(feature_values[index].strip(), float))
                else:
                    important_features.append(SafeCastUtil.safeCast(feature_values[index].strip(), str))
        return important_features

    def handleAlgorithmConfigs(self, arguments):
        algos = SupportedMachineLearningAlgorithms.fetchAlgorithms()
        configs = {}
        default_inner_perms = self.fetchOrReturnDefault(arguments.get(self.INNER_MONTE_CARLO_PERMUTATIONS), int, 10)
        default_outer_perms = self.fetchOrReturnDefault(arguments.get(self.OUTER_MONTE_CARLO_PERMUTATIONS), int, 10)

        for algo in algos:
            algo_config = arguments.get(algo)
            if algo_config is None:
                configs[algo] = [True, default_inner_perms, default_outer_perms]
            else:
                config_split = [param.strip() for param in algo_config.split(",")]
                if len(config_split) >= 3:
                    configs[algo] = [config_split[0] == 'True',
                                     SafeCastUtil.safeCast(config_split[1], int),
                                     SafeCastUtil.safeCast(config_split[2], int)]
        return configs

    def createRSENConfig(self, arguments, binary_cat_matrix):
        rsen_p_val = self.fetchOrReturnDefault(arguments.get(self.RSEN_P_VAL), float, 0.0)
        rsen_k_val = self.fetchOrReturnDefault(arguments.get(self.RSEN_P_VAL), float, 0.1)
        rsen_combine_gene_lists = self.fetchOrReturnDefault(arguments.get(self.RSEN_COMBINE_GENE_LISTS), bool, False)
        rsen_config = RSENConfig(binary_cat_matrix, rsen_p_val, rsen_k_val, rsen_combine_gene_lists)
        return rsen_config

    def createIndividualTrainConfig(self, arguments):
        individual_train_algorithm = self.fetchOrReturnDefault(arguments.get(self.INDIVIDUAL_TRAIN_ALGORITHM), str,
                                                               None)
        individual_train_hyperparams = self.fetchOrReturnDefault(arguments.get(self.INDIVIDUAL_TRAIN_HYPERPARAMS), str,
                                                                 "")
        individual_train_feature_gene_list_combo = self.fetchOrReturnDefault(
            arguments.get(self.INDIVIDUAL_TRAIN_FEATURE_GENE_LIST_COMBO),
            str, None)
        individual_train_config = IndividualTrainConfig(individual_train_algorithm, individual_train_hyperparams,
                                                        individual_train_feature_gene_list_combo)
        return individual_train_config

    def fetchBinaryCatMatrixIfApplicable(self, arguments, gene_lists, results_list, analyze_all):
        binary_matrix_file = arguments.get(ArgumentProcessingService.BINARY_CATEGORICAL_MATRIX)
        if binary_matrix_file is not None:
            if analyze_all:
                return self.createAndValidateFullFeatureMatrix(results_list, [binary_matrix_file])
            return self.createAndValidateFeatureMatrix(results_list, gene_lists, False, [binary_matrix_file])
        else:
            return None

    def fetchOrReturnDefault(self, field, to_type, default):
        if field:
            if field.lower() == 'false' and to_type is bool:
                return False
            return SafeCastUtil.safeCast(field, to_type)
        else:
            return default

    def determineSpecificCombos(self, combos):
        if combos is None:
            return []
        return [combo.strip().replace("\"", "") for combo in combos.split(",")]

