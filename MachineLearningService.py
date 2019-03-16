import csv
import gc
import logging
import multiprocessing
import numpy
import os
import threading
import psutil

from joblib import Parallel, delayed

from ArgumentConfig.AnalysisType import AnalysisType
from ArgumentConfig.ProcessedArguments import ProcessedArguments
from HTMLWritingService import HTMLWritingService
from SupportedMachineLearningAlgorithms import SupportedMachineLearningAlgorithms
from ArgumentProcessingService import ArgumentProcessingService
from DataFormattingService import DataFormattingService
from Trainers.RandomForestTrainer import RandomForestTrainer
from Trainers.LinearSVMTrainer import LinearSVMTrainer
from Trainers.RadialBasisFunctionSVMTrainer import RadialBasisFunctionSVMTrainer
from Trainers.ElasticNetTrainer import ElasticNetTrainer
from Trainers.RidgeRegressionTrainer import RidgeRegressionTrainer
from Trainers.LassoRegressionTrainer import LassoRegressionTrainer
from Trainers.RandomSubsetElasticNetTrainer import RandomSubsetElasticNetTrainer
from Trainers.AbstractModelTrainer import AbstractModelTrainer
from Utilities.GeneListComboUtility import GeneListComboUtility
from Utilities.ModelTrainerFactory import ModelTrainerFactory
from Utilities.SafeCastUtil import SafeCastUtil


class MachineLearningService(object):

    log = logging.getLogger(__name__)
    logging.basicConfig()
    log.setLevel(logging.INFO)

    MAXIMUM_FEATURES_RECORDED = 20
    DELIMITER = " --- "

    #TODO: consider extracting these to a helper class.
    SCORE_AND_HYPERPARAM_PHRASE = "score and optimal hyperparams for outer perm "
    FEATURE_FILE_GENE_LIST_COMBO = "feature file: gene list combo"
    R_SQUARED_SCORE = "R^2 score"
    PERCENT_ACCURATE_PREDICTIONS = "percentage accurate predictions"

    def __init__(self, data):
        self.inputs = data

    def analyze(self, input_folder):
        gene_list_combos = self.determineGeneListCombos()

        is_classifier = self.inputs.is_classifier
        analysis_type = self.inputs.analysisType()

        if analysis_type is AnalysisType.INDIVIDUAL_TRAIN:
            self.analyzeIndividualGeneListCombo(gene_list_combos, input_folder, is_classifier)
        elif analysis_type is AnalysisType.FULL_CLA_SPECIFIC_COMBO:
            self.analyzeGeneListCombos(self.determineSpecificCombos(gene_list_combos), input_folder, is_classifier)
        else:
            self.analyzeGeneListCombos(gene_list_combos, input_folder, is_classifier)

    def determineGeneListCombos(self):
        feature_names = self.inputs.features.get(ArgumentProcessingService.FEATURE_NAMES)
        if self.inputs.analysisType() is AnalysisType.SPEARMAN_NO_GENE_LISTS:
            return [[feature_names]]

        gene_lists = self.inputs.gene_lists
        combos, expected_length = GeneListComboUtility.determineGeneListCombos(gene_lists, feature_names)
        if len(combos) != expected_length:
            self.log.warning("Unexpected number of combos detected, should be %s but instead created %s.\n%s",
                             expected_length, len(combos), combos)
        return combos

    def analyzeIndividualGeneListCombo(self, gene_list_combos, input_folder, is_classifier):
        config = self.inputs.individual_train_config
        target_combo = config.combo
        target_algorithm = config.algorithm
        hyperparams = config.hyperparams.split(",")
        casted_params = [SafeCastUtil.safeCast(param, float) for param in hyperparams]

        rsen_config = self.inputs.rsen_config

        outer_monte_carlo_loops = self.inputs.outer_monte_carlo_permutations
        for gene_list_combo in gene_list_combos:
            plain_text_name = self.generateFeatureSetString(gene_list_combo)
            if plain_text_name == target_combo:
                trainer = ModelTrainerFactory.createTrainerFromTargetAlgorithm(is_classifier, target_algorithm, rsen_config)
                for permutation in range(0, outer_monte_carlo_loops):
                    results = self.inputs.results
                    formatted_data = self.formatData(self.inputs, True, True)
                    training_matrix = self.trimMatrixByFeatureSet(DataFormattingService.TRAINING_MATRIX,
                                                                  gene_list_combo, formatted_data)
                    testing_matrix = self.trimMatrixByFeatureSet(DataFormattingService.TESTING_MATRIX, gene_list_combo,
                                                                 formatted_data)
                    features, relevant_results = trainer.populateFeaturesAndResultsByCellLine(training_matrix, results)
                    feature_names = training_matrix.get(ArgumentProcessingService.FEATURE_NAMES)
                    model = trainer.buildModel(relevant_results, features, casted_params, feature_names)
                    model_score = trainer.fetchPredictionsAndScore(model, testing_matrix, results)
                    score = model_score[0]
                    accuracy = model_score[1]
                    importances = trainer.fetchFeatureImportances(model, gene_list_combo)
                    for key in importances.keys():
                        importances[key] = [importances[key]]
                    ordered_importances = self.averageAndSortImportances(importances, 1)
                    ordered_phrases = self.averageAndSortImportantRSENPhrases(
                                                            trainer.fetchModelPhrases(model, gene_list_combo), trainer)

                    numbered_combo = target_combo + " RUN " + SafeCastUtil.safeCast(permutation, str)
                    self.log.debug("Final score and accuracy of individual analysis for feature gene combo %s "
                                   "using algorithm %s: %s, %s", numbered_combo, target_algorithm, score, accuracy)
                    score_and_hyperparam = [self.generateScoreAndHyperParam(score, hyperparams, trainer)]
                    line = self.generateLine(accuracy, numbered_combo, ordered_importances, ordered_phrases, score,
                                             score_and_hyperparam)
                    self.writeToCSVInLock(line, input_folder, target_algorithm, outer_monte_carlo_loops, 1)
                return
        self.log.info("Gene list feature file %s combo not found in current dataset.", target_combo)
        return

    def generateLine(self, accuracy, combo, ordered_importances, ordered_phrases, score, score_and_hyperparam):
        return numpy.concatenate([[combo, score, accuracy], score_and_hyperparam,
                                  ordered_importances, ordered_phrases])

    def generateScoreAndHyperParam(self, score, hyperparam, trainer):
        hyperparam_string = ""
        keys = SafeCastUtil.safeCast(trainer.hyperparameters.keys(), list)
        for i in range(0, len(keys)):
            hyperparam_string += (keys[i] + ": " + SafeCastUtil.safeCast(hyperparam[i], str))
            if i < len(keys) - 1:
                hyperparam_string += ", "

        return SafeCastUtil.safeCast(score, str) + self.DELIMITER + hyperparam_string

    def shouldTrainAlgorithm(self, algorithm):
        configs = self.inputs.algorithm_configs
        return configs is not None and configs.get(algorithm) is not None and configs.get(algorithm)[0]

    def determineSpecificCombos(self, all_combos):
        specific_combos = self.inputs.specific_combos
        selected_combos = {}
        for specific_combo in specific_combos:
            for combo in all_combos:
                combo_string = self.generateFeatureSetString(combo)
                if specific_combo == combo_string and selected_combos.get(combo_string) is None:
                    selected_combos[combo_string] = combo
                else:
                    this_sorted_combo = sorted(combo_string.split(" "))
                    specific_sorted_combo = sorted(specific_combo.split(" "))
                    if this_sorted_combo == specific_sorted_combo and selected_combos.get(combo_string) is None:
                        selected_combos[combo_string] = combo
        selected_combo_names = SafeCastUtil.safeCast(selected_combos.keys(), list)
        if len(selected_combo_names) < len(specific_combos):
            self.log.warning("Not all specified combos were available in this data folder.\n"
                             "Specified combos: %s\n Selected combos: %s", specific_combos, selected_combo_names)
        else:
            self.log.info("Only running analysis on following combos:\n %s", selected_combo_names)
        return SafeCastUtil.safeCast(selected_combos.values(), list)

    def analyzeGeneListCombos(self, gene_list_combos, input_folder, is_classifier):
        enet = SupportedMachineLearningAlgorithms.ELASTIC_NET
        rig_reg = SupportedMachineLearningAlgorithms.RIDGE_REGRESSION
        las_reg = SupportedMachineLearningAlgorithms.LASSO_REGRESSION
        lin_svm = SupportedMachineLearningAlgorithms.LINEAR_SVM
        rbf = SupportedMachineLearningAlgorithms.RADIAL_BASIS_FUNCTION_SVM
        rf = SupportedMachineLearningAlgorithms.RANDOM_FOREST
        rsen = SupportedMachineLearningAlgorithms.RANDOM_SUBSET_ELASTIC_NET

        #TODO: Extract all these to the new ModelTrainerFactory class.
        if not is_classifier and self.shouldTrainAlgorithm(enet):
            elasticnet_trainer = ElasticNetTrainer(is_classifier)
            elasticnet_trainer.logTrainingMessage(self.monteCarloPermsByAlgorithm(enet, True),
                                                  self.monteCarloPermsByAlgorithm(enet, False),
                                                  len(gene_list_combos))
            self.handleParallellization(gene_list_combos, input_folder, elasticnet_trainer)

        if not is_classifier and self.shouldTrainAlgorithm(rig_reg):
            ridge_regression_trainer = RidgeRegressionTrainer(is_classifier)
            ridge_regression_trainer.logTrainingMessage(self.monteCarloPermsByAlgorithm(rig_reg, True),
                                                        self.monteCarloPermsByAlgorithm(rig_reg, False),
                                                        len(gene_list_combos))
            self.handleParallellization(gene_list_combos, input_folder, ridge_regression_trainer)

        if not is_classifier and self.shouldTrainAlgorithm(las_reg):
            lasso_regression_trainer = LassoRegressionTrainer(is_classifier)
            lasso_regression_trainer.logTrainingMessage(self.monteCarloPermsByAlgorithm(rig_reg, True),
                                                        self.monteCarloPermsByAlgorithm(rig_reg, False),
                                                        len(gene_list_combos))
            self.handleParallellization(gene_list_combos, input_folder, lasso_regression_trainer)

        if self.shouldTrainAlgorithm(lin_svm):
            linear_svm_trainer = LinearSVMTrainer(is_classifier)
            linear_svm_trainer.logTrainingMessage(self.monteCarloPermsByAlgorithm(lin_svm, True),
                                                  self.monteCarloPermsByAlgorithm(lin_svm, False),
                                                  len(gene_list_combos))
            self.handleParallellization(gene_list_combos, input_folder, linear_svm_trainer)

        if self.shouldTrainAlgorithm(rbf):
            rbf_svm_trainer = RadialBasisFunctionSVMTrainer(is_classifier)
            rbf_svm_trainer.logTrainingMessage(self.monteCarloPermsByAlgorithm(rbf, True),
                                               self.monteCarloPermsByAlgorithm(rbf, False), len(gene_list_combos))
            self.handleParallellization(gene_list_combos, input_folder, rbf_svm_trainer)

        if self.shouldTrainAlgorithm(rf):
            rf_trainer = RandomForestTrainer(is_classifier)
            rf_trainer.logTrainingMessage(self.monteCarloPermsByAlgorithm(rf, True),
                                          self.monteCarloPermsByAlgorithm(rf, False), len(gene_list_combos))
            self.handleParallellization(gene_list_combos, input_folder, rf_trainer)

        rsen_config = self.inputs.rsen_config
        binary_cat_matrix = rsen_config.binary_cat_matrix
        if self.shouldTrainAlgorithm(rsen) and not is_classifier and binary_cat_matrix is not None:
            p_val = rsen_config.p_val
            k_val = rsen_config.k_val
            rsen_trainer = RandomSubsetElasticNetTrainer(is_classifier, binary_cat_matrix, p_val, k_val)
            rsen_trainer.logTrainingMessage(self.monteCarloPermsByAlgorithm(rsen, True),
                                            self.monteCarloPermsByAlgorithm(rsen, False),
                                            len(gene_list_combos))
            self.handleParallellization(gene_list_combos, input_folder, rsen_trainer)

    def monteCarloPermsByAlgorithm(self, algorithm, outer):
        monte_carlo_config = self.inputs.algorithm_configs.get(algorithm)
        return monte_carlo_config[1] if outer else monte_carlo_config[2]

    def handleParallellization(self, gene_list_combos, input_folder, trainer):
        max_nodes = multiprocessing.cpu_count()
        requested_threads = self.inputs.num_threads
        nodes_to_use = numpy.amin([requested_threads, max_nodes])

        valid_combos = self.fetchValidGeneListCombos(input_folder, gene_list_combos, trainer)

        Parallel(n_jobs=nodes_to_use)(delayed(self.runMonteCarloSelection)(feature_set, trainer, input_folder,
                                                                           len(valid_combos))
                                      for feature_set in valid_combos)
        self.logMemoryUsageAndGarbageCollect()

    def fetchValidGeneListCombos(self, input_folder, gene_list_combos, trainer):
        valid_combos = [feature_set for feature_set in gene_list_combos if trainer.shouldProcessFeatureSet(feature_set)]

        rsen_config = self.inputs.rsen_config
        if trainer.algorithm == SupportedMachineLearningAlgorithms.RANDOM_SUBSET_ELASTIC_NET and \
                rsen_config.combine_gene_lists:
            all_genes = GeneListComboUtility.fetchAllGeneListGenesDeduped(self.inputs.gene_lists)
            # TODO: Can fail if "." in feature name.
            bin_cat_matrix = rsen_config.binary_cat_matrix.get(ArgumentProcessingService.FEATURE_NAMES)[0].split(".")[0]
            full_gene_list = [bin_cat_matrix + "." + gene for gene in all_genes if len(gene.strip()) > 0]

            new_combos = []
            for combo in valid_combos:
                new_combo = []
                for feature_set in combo:
                    if bin_cat_matrix in feature_set[0]:
                        new_combo.append(full_gene_list)
                    else:
                        new_combo.append(feature_set)
                if new_combo not in new_combos:
                    new_combos.append(new_combo)
            return self.trimAnalyzedCombos(input_folder, new_combos, trainer)
        else:
            return self.trimAnalyzedCombos(input_folder, valid_combos, trainer)

    def trimAnalyzedCombos(self, input_folder, valid_combos, trainer):
        file_name = trainer.algorithm + ".csv"
        if file_name not in os.listdir(input_folder):
            return valid_combos

        existing_combo_strings = []
        with open(input_folder + "/" + file_name) as analyzed_file:
            try:
                for line_index, line in enumerate(analyzed_file):
                    if line_index == 0:
                        continue
                    existing_combo_strings.append(line.strip().split(",")[0])
            except ValueError as error:
                self.log.error("Error reading existing combos from analysis file: %s", analyzed_file, error)
            finally:
                analyzed_file.close()

        trimmed_combos = []
        for combo in valid_combos:
            if self.generateFeatureSetString(combo) not in existing_combo_strings:
                trimmed_combos.append(combo)
        return trimmed_combos

    def runMonteCarloSelection(self, feature_set, trainer, input_folder, num_combos):
        scores = []
        accuracies = []
        importances = {}
        feature_set_as_string = self.generateFeatureSetString(feature_set)
        outer_perms = self.monteCarloPermsByAlgorithm(trainer.algorithm, True)
        important_rsen_phrases = {}
        scores_and_hyperparams = []

        for i in range(1, outer_perms + 1):
            formatted_data = self.formatData(self.inputs, True, True)
            training_matrix = self.trimMatrixByFeatureSet(DataFormattingService.TRAINING_MATRIX, feature_set,
                                                          formatted_data)
            testing_matrix = self.trimMatrixByFeatureSet(DataFormattingService.TESTING_MATRIX, feature_set,
                                                         formatted_data)

            self.log.debug("Computing outer Monte Carlo Permutation %s for %s.", i, feature_set_as_string)

            optimal_hyperparams = self.determineOptimalHyperparameters(feature_set, formatted_data, trainer)
            record_diagnostics = self.inputs.record_diagnostics
            trainer.logIfBestHyperparamsOnRangeThreshold(optimal_hyperparams, record_diagnostics, input_folder)
            trainer.logOptimalHyperParams(optimal_hyperparams, self.generateFeatureSetString(feature_set),
                                          record_diagnostics, input_folder)

            prediction_data = self.fetchOuterPermutationModelScore(feature_set, trainer,
                                                                   optimal_hyperparams, testing_matrix,
                                                                   training_matrix)
            scores.append(prediction_data[0])
            accuracies.append(prediction_data[1])
            for importance in prediction_data[2].keys():
                if importances.get(importance) is not None:
                    importances[importance].append(prediction_data[2].get(importance))
                else:
                    importances[importance] = [prediction_data[2].get(importance)]
            if len(prediction_data) == 4 and \
                    trainer.algorithm == SupportedMachineLearningAlgorithms.RANDOM_SUBSET_ELASTIC_NET:
                for phrase in prediction_data[3].keys():
                    if important_rsen_phrases.get(phrase) is not None:
                        important_rsen_phrases[phrase].append(prediction_data[3].get(phrase))
                    else:
                        important_rsen_phrases[phrase] = [prediction_data[3].get(phrase)]
            scores_and_hyperparams.append(self.generateScoreAndHyperParam(prediction_data[0], optimal_hyperparams,
                                                                          trainer))

            self.logMemoryUsageAndGarbageCollect()

        average_score = numpy.mean(scores)
        average_accuracy = numpy.mean(accuracies)
        self.log.debug("Average score and accuracy of all Monte Carlo runs for %s: %s, %s",
                       feature_set_as_string, average_score, average_accuracy)
        ordered_importances = self.averageAndSortImportances(importances, outer_perms)

        ordered_phrases = self.averageAndSortImportantRSENPhrases(important_rsen_phrases, trainer)

        line = self.generateLine(average_accuracy, feature_set_as_string, ordered_importances, ordered_phrases,
                                 average_score, scores_and_hyperparams)
        self.writeToCSVInLock(line, input_folder, trainer.algorithm, num_combos, outer_perms)
        self.saveOutputToTxtFile(scores, accuracies, feature_set_as_string, input_folder, trainer.algorithm)

    def generateFeatureSetString(self, feature_set):
        return GeneListComboUtility.generateFeatureSetString(feature_set, self.inputs.gene_lists,
                                                             self.inputs.rsen_config.combine_gene_lists,
                                                             self.inputs.analysisType)

    def fetchOuterPermutationModelScore(self, feature_set, trainer, optimal_hyperparams, testing_matrix,
                                        training_matrix):
        # TODO: Handle hyperparams with n
        results = self.inputs.results
        features, relevant_results = trainer.populateFeaturesAndResultsByCellLine(training_matrix, results)
        feature_names = training_matrix.get(ArgumentProcessingService.FEATURE_NAMES)
        model = trainer.buildModel(relevant_results, features, optimal_hyperparams, feature_names)
        score, accuracy = trainer.fetchPredictionsAndScore(model, testing_matrix, results)
        return [score, accuracy, trainer.fetchFeatureImportances(model, feature_set),
                trainer.fetchModelPhrases(model, feature_set)]

    def averageAndSortImportances(self, importances, outer_loops):
        for key in importances.keys():
            if len(importances[key]) != outer_loops:
                self.log.warning("Different amount of importances for feature %s than expected. Should be %s but is "
                                 "instead %s.", key, outer_loops, len(importances[key]))
        ordered = []
        [ordered.append({"feature": key, "importance": numpy.sum(importances[key]) / outer_loops}) for key in
         importances.keys()]
        ordered = sorted(ordered, key=lambda k: k["importance"], reverse=True)
        trimmed = ordered[:self.MAXIMUM_FEATURES_RECORDED]

        final_imps = []
        for i in range(0, self.MAXIMUM_FEATURES_RECORDED):
            if i < len(trimmed):
                summary = trimmed[i].get("feature") + self.DELIMITER + \
                          SafeCastUtil.safeCast(trimmed[i].get("importance"), str)
                final_imps.append(summary)
            else:
                final_imps.append("")
        return final_imps

    def averageAndSortImportantRSENPhrases(self, important_rsen_phrases, trainer):
        if trainer.algorithm == SupportedMachineLearningAlgorithms.RANDOM_SUBSET_ELASTIC_NET:
            ordered_phrases = []
            [ordered_phrases.append({"phrase": key, "score": numpy.average(important_rsen_phrases[key])}) for key
             in important_rsen_phrases.keys()]
            ordered_phrases = sorted(ordered_phrases, key=lambda k: k["score"], reverse=True)
            trimmed = ordered_phrases[:self.MAXIMUM_FEATURES_RECORDED]
            final_phrases = []
            for i in range(0, self.MAXIMUM_FEATURES_RECORDED):
                if i < len(trimmed):
                    summary = trimmed[i].get("phrase") + self.DELIMITER + \
                              SafeCastUtil.safeCast(trimmed[i].get("score"), str)
                    final_phrases.append(summary)
                else:
                    final_phrases.append("")
            return final_phrases
        else:
            return []

    def determineInnerHyperparameters(self, feature_set, formatted_data, trainer):
        inner_model_hyperparams = {}
        inner_perms = self.monteCarloPermsByAlgorithm(trainer.algorithm, False)
        for j in range(1, inner_perms + 1):
            self.logMemoryUsageAndGarbageCollect()
            formatted_inputs = self.reformatInputsByTrainingMatrix(
                formatted_data.get(DataFormattingService.TRAINING_MATRIX),
                formatted_data.get(ArgumentProcessingService.FEATURE_NAMES))
            further_formatted_data = self.formatData(formatted_inputs, False, False)
            inner_validation_matrix = self.trimMatrixByFeatureSet(DataFormattingService.TESTING_MATRIX, feature_set,
                                                                  further_formatted_data)
            inner_train_matrix = self.trimMatrixByFeatureSet(DataFormattingService.TRAINING_MATRIX, feature_set,
                                                             further_formatted_data)
            model_data = trainer.hyperparameterize(inner_train_matrix, inner_validation_matrix, self.inputs.results)
            for data in model_data.keys():
                if inner_model_hyperparams.get(data) is not None:
                    inner_model_hyperparams[data].append(model_data[data])
                else:
                    inner_model_hyperparams[data] = [model_data[data]]
        return inner_model_hyperparams

    def logMemoryUsageAndGarbageCollect(self):
        processes = psutil.Process()
        procs = [processes] + processes.children(recursive=True)
        for process in procs:
            rss = process.memory_info().rss
            memory_usage_mb = numpy.round(rss / 1e6, 2)
            self.log.debug("Memory usage for PID %s: %s: MB", process.pid, memory_usage_mb)
        gc.collect()

    def formatData(self, inputs, should_scale, should_one_hot_encode):
        data_formatting_service = DataFormattingService(inputs)
        return data_formatting_service.formatData(should_scale, should_one_hot_encode)

    def reformatInputsByTrainingMatrix(self, training_matrix, feature_names):
        real_inputs = self.inputs

        features = {}
        results = []
        features[ArgumentProcessingService.FEATURE_NAMES] = feature_names
        for training_cell in training_matrix.keys():
            for input_cell in real_inputs.features.keys():
                if training_cell is input_cell:
                    features[training_cell] = training_matrix.get(training_cell)
                    for result in real_inputs.results:
                        if result[0] is training_cell:
                            results.append(result)
                            break
                    break
        return ProcessedArguments(results, real_inputs.is_classifier, features, real_inputs.gene_lists,
                                  real_inputs.inner_monte_carlo_permutations,
                                  real_inputs.outer_monte_carlo_permutations, real_inputs.data_split,
                                  real_inputs.algorithm_configs, real_inputs.num_threads,
                                  real_inputs.record_diagnostics,
                                  real_inputs.individual_train_config, real_inputs.rsen_config,
                                  real_inputs.specific_combos, False)

    def determineOptimalHyperparameters(self, feature_set, formatted_data, trainer):
        inner_model_hyperparams = self.determineInnerHyperparameters(feature_set, formatted_data, trainer)
        highest_average = trainer.DEFAULT_MIN_SCORE
        best_hyperparam = None
        for hyperparam_set in inner_model_hyperparams.keys():
            if hyperparam_set == AbstractModelTrainer.ADDITIONAL_DATA:
                continue
            average = numpy.average([results[0] for results in inner_model_hyperparams[hyperparam_set]])  # raw score
            if average > highest_average:
                best_hyperparam = SafeCastUtil.safeCast(hyperparam_set, list)
                highest_average = average
        additional_data = inner_model_hyperparams.get(AbstractModelTrainer.ADDITIONAL_DATA)
        if additional_data:
            best_hyperparam.append(additional_data)
        return best_hyperparam

    def trimMatrixByFeatureSet(self, matrix_type, gene_lists, formatted_inputs):
        full_matrix = formatted_inputs.get(matrix_type)
        trimmed_matrix = {
            ArgumentProcessingService.FEATURE_NAMES: []
        }

        important_indices = []
        feature_names = formatted_inputs.get(ArgumentProcessingService.FEATURE_NAMES)
        for i in range(0, len(feature_names)):
            for gene_list in gene_lists:
                for gene in gene_list:
                    if gene == feature_names[i]:
                        important_indices.append(i)
                        trimmed_matrix[ArgumentProcessingService.FEATURE_NAMES].append(gene)

        for cell_line in full_matrix.keys():
            new_cell_line_features = []
            for j in range(0, len(full_matrix[cell_line])):
                if j in important_indices:
                    new_cell_line_features.append(full_matrix[cell_line][j])
            trimmed_matrix[cell_line] = new_cell_line_features
        return trimmed_matrix

    def writeToCSVInLock(self, line, input_folder, ml_algorithm, num_combos, outer_perms):
        lock = threading.Lock()
        lock.acquire(True)
        self.lockThreadMessage()

        file_name = ml_algorithm + ".csv"
        write_action = "w"
        if file_name in os.listdir(input_folder):
            write_action = "a"
        with open(input_folder + "/" + file_name, write_action, newline='') as csv_file:
            try:
                writer = csv.writer(csv_file)
                if write_action == "w":
                    writer.writerow(self.getCSVFileHeader(self.inputs.is_classifier, ml_algorithm, outer_perms))
                writer.writerow(line)
            except ValueError as error:
                self.log.error("Error writing to file %s. %s", file_name, error)
            finally:
                csv_file.close()

        total_lines = 0
        with open(input_folder + "/" + file_name) as csv_file:
            try:
                reader = csv.reader(csv_file)
                total_lines += (len(SafeCastUtil.safeCast(reader, list)) - 1)
            except ValueError as error:
                self.log.error("Error reading lines from file %s. %s", file_name, error)
            finally:
                csv_file.close()
                self.logPercentDone(total_lines, num_combos, ml_algorithm)

        self.unlockThreadMessage()
        lock.release()

    @staticmethod
    def getCSVFileHeader(is_classifier, ml_algorithm, outer_perms):
        header = [MachineLearningService.FEATURE_FILE_GENE_LIST_COMBO]
        if is_classifier:
            header.append(MachineLearningService.PERCENT_ACCURATE_PREDICTIONS)
            header.append("accuracy score")
        else:
            header.append(MachineLearningService.R_SQUARED_SCORE)
            header.append("mean squared error")
        for i in range(1, outer_perms + 1):
            header.append(MachineLearningService.SCORE_AND_HYPERPARAM_PHRASE + SafeCastUtil.safeCast(i, str))
        if ml_algorithm == SupportedMachineLearningAlgorithms.RADIAL_BASIS_FUNCTION_SVM:
            return header
        feature_analysis = " most important feature"
        for i in range(1, MachineLearningService.MAXIMUM_FEATURES_RECORDED + 1):
            suffix = MachineLearningService.generateNumericalSuffix(i)
            header.append(SafeCastUtil.safeCast(i, str) + suffix + feature_analysis)
        if ml_algorithm == SupportedMachineLearningAlgorithms.RANDOM_SUBSET_ELASTIC_NET:
            phrase_analysis = " most significant boolean phrase"
            for i in range(1, MachineLearningService.MAXIMUM_FEATURES_RECORDED + 1):
                suffix = MachineLearningService.generateNumericalSuffix(i)
                header.append(SafeCastUtil.safeCast(i, str) + suffix + phrase_analysis)
        return header

    @staticmethod
    def generateNumericalSuffix(i):
        if i == 1 and i != 11:
            return "st"
        elif i == 2 and i != 12:
            return "nd"
        elif i == 3 and i != 13:
            return "rd"
        else:
            return "th"

    def logPercentDone(self, total_lines, num_combos, ml_algorithm):
        percent_done = numpy.round((total_lines / num_combos) * 100, 1)
        percentage_bar = "["
        for i in range(0, 100):
            if i < percent_done:
                percentage_bar += "="
            elif i >= percent_done:
                if i > 0 and (i - 1) < percent_done:
                    percentage_bar += ">"
                else:
                    percentage_bar += " "
        percentage_bar += "]"
        self.log.info("Total progress for %s: %s%% done: %s", ml_algorithm, percent_done, percentage_bar)

    def saveOutputToTxtFile(self, scores, accuracies, feature_set_as_string, input_folder, algorithm):
        lock = threading.Lock()
        self.lockThreadMessage()
        lock.acquire(True)

        file_name = HTMLWritingService.RECORD_FILE
        write_action = "w"
        if file_name in os.listdir(input_folder):
            write_action = "a"
        with open(input_folder + "/" + file_name, write_action) as output_file:
            try:
                output_file.write(algorithm + MachineLearningService.DELIMITER + feature_set_as_string +
                                  MachineLearningService.DELIMITER + SafeCastUtil.safeCast(scores, str) +
                                  MachineLearningService.DELIMITER + SafeCastUtil.safeCast(accuracies, str)
                                  + "\n")
            except ValueError as error:
                self.log.error("Error saving output of %s analysis to memory: %s", algorithm, error)
            finally:
                output_file.close()
                self.unlockThreadMessage()
                lock.release()

    def lockThreadMessage(self):
        self.log.debug("Locking current thread %s.", threading.current_thread())

    def unlockThreadMessage(self):
        self.log.debug("Releasing current thread %s.", threading.current_thread())
