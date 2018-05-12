from abc import ABC, abstractmethod
import logging
from Utilities.SafeCastUtil import SafeCastUtil
import numpy
from sklearn.metrics import mean_squared_error
from sklearn.metrics import accuracy_score

class AbstractModelTrainer(ABC):

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    DEFAULT_MIN_SCORE = -10

    @abstractmethod
    def __init__(self, algorithm, hyperparameters, is_classifier):
        self.algorithm = algorithm
        self.hyperparameters = hyperparameters
        self.is_classifier = is_classifier

    @abstractmethod
    def hyperparameterize(self, training_matrix, testing_matrix, results):
        pass

    @abstractmethod
    def train(self, results, features, hyperparams):
        pass

    @abstractmethod
    def setModelDataDictionary(self, model_data, hyperparam_set, current_model_score):
        pass

    @abstractmethod
    def logOptimalHyperParams(self, hyperparams, feature_set_as_string):
        pass

    def trainingMessage(self, outer_monte_carlo_perms, inner_monte_carlo_perms, num_gene_list_combos):
        num_models = self.determineNumModelsToCreate(outer_monte_carlo_perms, inner_monte_carlo_perms, num_gene_list_combos)
        self.log.info("Running permutations on %s different combinations of features. Requires creation of %s "
                      "different %s models.", SafeCastUtil.safeCast(num_gene_list_combos, str),
                      num_models, self.algorithm)

    def determineNumModelsToCreate(self, outer_monte_carlo_perms, inner_monte_carlo_perms, num_gene_list_combos):
        num_models = outer_monte_carlo_perms * inner_monte_carlo_perms * num_gene_list_combos
        for hyperparam_set in self.hyperparameters.values():
            num_models *= len(hyperparam_set)
        return num_models

    def loopThroughHyperparams(self, hyperparams, training_matrix, testing_matrix, results):
        features, relevant_results = self.populateFeaturesAndResultsByCellLine(training_matrix, results)

        model_data = {}
        for hyperparam_set in self.fetchAllHyperparamPermutations(hyperparams):
            model = self.train(relevant_results, features, hyperparam_set)
            current_model_score = self.fetchPredictionsAndScore(model, testing_matrix, results)
            self.setModelDataDictionary(model_data, hyperparam_set, current_model_score)
        return model_data

    def fetchAllHyperparamPermutations(self, hyperparams):
        all_perms = []
        hyperparam_keys = SafeCastUtil.safeCast(hyperparams.keys(), list)
        zero_filled_indices = SafeCastUtil.safeCast(numpy.zeros(len(hyperparam_keys)), list)
        target_index = len(zero_filled_indices) - 1
        current_perm = zero_filled_indices[:]
        while target_index >= 0:
            current_hyperparams = []
            for i in range(0, len(current_perm)):
                current_hyperparams.append(hyperparams[hyperparam_keys[i]][SafeCastUtil.safeCast(current_perm[i], int)])
            if current_hyperparams not in all_perms:
                clone_array = current_hyperparams[:]
                all_perms.append(clone_array)

            if current_perm[target_index] < len(hyperparams[hyperparam_keys[target_index]]) - 1:
                current_perm[target_index] += 1
                while len(current_perm) > target_index + 1 and current_perm[target_index + 1] < len(hyperparams[hyperparam_keys[target_index]]):
                    target_index += 1
            else:
                target_index -= 1
                for subsequent_index in range(target_index, len(current_perm) - 1):
                    current_perm[subsequent_index + 1] = 0
        return all_perms

    def fetchPredictionsAndScore(self, model, testing_matrix, results):
        if model is None:
            return self.DEFAULT_MIN_SCORE
        features, relevant_results = self.populateFeaturesAndResultsByCellLine(testing_matrix, results)
        predictions = model.predict(features)
        score = model.score(features, relevant_results)
        if self.is_classifier:
            accuracy = accuracy_score(relevant_results, predictions)
        else:
            accuracy = mean_squared_error(relevant_results, predictions)
        del model
        return score, accuracy

    def populateFeaturesAndResultsByCellLine(self, matrix, results):
        features = []
        relevant_results = []
        for cell in matrix.keys():
            features.append(matrix[cell])
            for result in results:
                if result[0] == cell:
                    relevant_results.append(result[1])
        return features, relevant_results
