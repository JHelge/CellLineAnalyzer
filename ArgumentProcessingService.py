import multiprocessing
import os
import logging
import re

from Utilities.SafeCastUtil import SafeCastUtil
from SupportedMachineLearningAlgorithms import SupportedMachineLearningAlgorithms


class ArgumentProcessingService(object):

    log = logging.getLogger(__name__)
    logging.basicConfig()
    log.setLevel(logging.INFO)

    ARGUMENTS_FILE = "arguments.txt"
    GENE_LISTS = "gene_list"

    RESULTS = "results"
    IS_CLASSIFIER = "is_classifier"
    FEATURES = "features"
    FEATURE_NAMES = "featureNames"
    INNER_MONTE_CARLO_PERMUTATIONS = "inner_monte_carlo_permutations"
    OUTER_MONTE_CARLO_PERMUTATIONS = "outer_monte_carlo_permutations"
    DATA_SPLIT = "data_split"
    NUM_THREADS = "num_threads"
    SKIP_RF = "skip_rf"
    SKIP_LINEAR_SVM = "skip_linear_svm"
    SKIP_RBF_SVM = "skip_rbf_svm"
    SKIP_ELASTIC_NET = "skip_elastic_net"
    SKIP_LINEAR_REGRESSION = "skip_linear_regression"
    RECORD_DIAGNOSTICS = "record_diagnostics"

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
        if is_classifier is not None and results_file is not None:
            results_list = self.validateAndExtractResults(results_file, is_classifier)
            gene_lists = self.extractGeneList()
            write_diagnostics = self.fetchOrReturnDefault(arguments.get(self.RECORD_DIAGNOSTICS), bool, False)
            feature_map = self.createAndValidateFeatureMatrix(results_list, results_file, gene_lists, write_diagnostics)
            if feature_map and gene_lists and results_list:
                return {
                    self.RESULTS: results_list,
                    self.IS_CLASSIFIER: is_classifier,
                    self.FEATURES: feature_map,
                    self.GENE_LISTS: gene_lists,
                    self.INNER_MONTE_CARLO_PERMUTATIONS: self.fetchOrReturnDefault(arguments.get(self.INNER_MONTE_CARLO_PERMUTATIONS), int, 10),
                    self.OUTER_MONTE_CARLO_PERMUTATIONS: self.fetchOrReturnDefault(arguments.get(self.OUTER_MONTE_CARLO_PERMUTATIONS), int, 10),
                    self.DATA_SPLIT: self.fetchOrReturnDefault(arguments.get(self.DATA_SPLIT), float, 0.8),
                    self.SKIP_RF: self.fetchOrReturnDefault(arguments.get(self.SKIP_RF), bool, False),
                    self.SKIP_LINEAR_SVM: self.fetchOrReturnDefault(arguments.get(self.SKIP_LINEAR_SVM), bool, False),
                    self.SKIP_RBF_SVM: self.fetchOrReturnDefault(arguments.get(self.SKIP_RBF_SVM), bool, False),
                    self.SKIP_ELASTIC_NET: self.fetchOrReturnDefault(arguments.get(self.SKIP_ELASTIC_NET), bool, False),
                    self.SKIP_LINEAR_REGRESSION: self.fetchOrReturnDefault(arguments.get(self.SKIP_LINEAR_REGRESSION),
                                                                           bool, False),
                    self.NUM_THREADS: self.fetchOrReturnDefault(arguments.get(self.NUM_THREADS), int,
                                                                multiprocessing.cpu_count()),
                    self.RECORD_DIAGNOSTICS: write_diagnostics,
                    self.INDIVIDUAL_TRAIN_ALGORITHM: self.fetchOrReturnDefault(arguments.get(self.INDIVIDUAL_TRAIN_ALGORITHM), str, None),
                    self.INDIVIDUAL_TRAIN_HYPERPARAMS: self.fetchOrReturnDefault(arguments.get(self.INDIVIDUAL_TRAIN_HYPERPARAMS), str, ""),
                    self.INDIVIDUAL_TRAIN_FEATURE_GENE_LIST_COMBO: self.fetchOrReturnDefault(arguments.get(self.INDIVIDUAL_TRAIN_FEATURE_GENE_LIST_COMBO), str, None)
                }
            else:
                return None
        else:
            return None

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

    def extractGeneList(self):
        gene_lists = {"null_gene_list": []}
        files = os.listdir(self.input_folder)
        for file in [f for f in files if self.GENE_LISTS in f]:
            file_path = self.input_folder + "/" + file
            with open(file_path) as gene_list_file:
                gene_lists[file.split(".csv")[0]] = gene_list_file.read().strip().split(",")

        return gene_lists

    def createAndValidateFeatureMatrix(self, results_list, results_file, gene_lists, write_diagnostics):
        files = os.listdir(self.input_folder)
        incomplete_features = []
        for file in [file for file in files if self.fileIsFeatureFile(file, results_file)]:
            features_path = self.input_folder + "/" + file
            incomplete_features.append([file, self.validateGeneLists(features_path, file, gene_lists)])
        if write_diagnostics:
            self.writeDiagnostics(incomplete_features)

        feature_matrix = {self.FEATURE_NAMES: []}
        for file in [file for file in files if self.fileIsFeatureFile(file, results_file)]:
            features_path = self.input_folder + "/" + file
            self.extractFeatureMatrix(feature_matrix, features_path, file, gene_lists, results_list)
        return feature_matrix

    def validateGeneLists(self, features_path, file, gene_lists):
        features_missing_from_files = {}
        with open(features_path) as feature_file:
            try:
                for line_index, line in enumerate(feature_file):
                    if line_index == 0:
                        feature_names = line.split(",")
                        features_missing_from_files = self.validateAndTrimGeneList(feature_names, gene_lists, file)
                    break
            except ValueError as valueError:
                self.log.error(valueError)
                return None
            finally:
                self.log.debug("Closing file %s", feature_file)
        return features_missing_from_files

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
                    diagnostics_file.write("\nFeatures from gene lists not available in " + feature_file[0] + ":\n")
                    for gene_list in feature_file[1].keys():
                        diagnostics_file.write("\tFeatures needed from gene list " + gene_list + ":\n")
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
                                if feature_names[i] == gene_list_feature:
                                    important_index = i
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
        rf_analysis = SupportedMachineLearningAlgorithms.RANDOM_FOREST + ".csv"
        linear_svm_analysis = SupportedMachineLearningAlgorithms.LINEAR_SVM + ".csv"
        rbf_svm_analysis = SupportedMachineLearningAlgorithms.RADIAL_BASIS_FUNCTION_SVM + ".csv"
        elastic_net_analysis = SupportedMachineLearningAlgorithms.ELASTIC_NET + ".csv"
        return file != results_file and file != self.ARGUMENTS_FILE and self.GENE_LISTS not in file and\
               file != rf_analysis and file != rbf_svm_analysis and file != linear_svm_analysis and\
               file != elastic_net_analysis and ".csv" in file.lower()

    def determineFeatureName(self, feature_name, file):
        return SafeCastUtil.safeCast(file.split(".")[0] + "." + feature_name.strip(), str)

    def extractCastedFeatures(self, line, important_feature_indices):
        important_features = []
        feature_values = line.split(",")
        for index in important_feature_indices:
            if index is None:
                # TODO: Verify that this is acceptable, it works for one hot encoding and should never vary in any model
                important_features.append(0)
            else:
                if SafeCastUtil.safeCast(feature_values[index], float) is not None:
                    important_features.append(SafeCastUtil.safeCast(feature_values[index].strip(), float))
                else:
                    important_features.append(SafeCastUtil.safeCast(feature_values[index].strip(), str))
        return important_features

    def fetchOrReturnDefault(self, field, to_type, default):
        if field:
            if field.lower() == 'false' and to_type is bool:
                return False
            return SafeCastUtil.safeCast(field, to_type)
        else:
            return default
