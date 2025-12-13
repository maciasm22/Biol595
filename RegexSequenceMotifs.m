%% Regular Expression Motif Finder with Cross-Validation
% Finds minimal regular expressions that distinguish aligned sequences
% from randomized controls using user-defined equivalence classes

clear all; close all; clc;

%% Configuration
config = struct();

% Equivalence classes (amino acid groupings)
config.equivClasses = {
    {'I','L','M','V'},  % Hydrophobic aliphatic
    {'D','E','N','Q'},  % Acidic/Amide
    {'H','K','R'},      % Basic
    {'S','T'},          % Hydroxyl
    {'A','G'},          % Small
    {'F','W','Y'},      % Aromatic
    {'C'},              % Cysteine
    {'P'}               % Proline
};

% Alternative: no equivalence (each amino acid separate)
config.noEquivClasses = arrayfun(@(x) {char(x)}, 'ACDEFGHIKLMNPQRSTVWY', 'UniformOutput', false);

% Parameters
config.trainingFraction = 0.8;
config.numFolds = 5;  % For cross-validation
config.minRegexLength = 3;
config.maxRegexLength = 15;
config.populationSize = 50;
config.numGenerations = 100;
config.mutationRate = 0.15;
config.crossoverRate = 0.7;

%% Menu
fprintf('\n=== Regular Expression Motif Finder ===\n');
fprintf('1. Run single experiment\n');
fprintf('2. Test effect of training set size\n');
fprintf('3. Compare equivalence classes vs. no equivalence\n');
fprintf('4. Full analysis (all tests)\n');
choice = input('Select option (1-4): ');

% Generate or load test alignment
fprintf('\nGenerating test multiple sequence alignment...\n');
alignment = generateTestAlignment(30, 50);  % 30 sequences, 50 positions

if choice == 1
    fprintf('\n--- Single Experiment ---\n');
    results = runSingleExperiment(alignment, config);
    displayResults(results, config);
    
elseif choice == 2
    fprintf('\n--- Training Set Size Analysis ---\n');
    fractions = [0.3, 0.5, 0.7, 0.8, 0.9];
    trainingSizeResults = testTrainingSizeEffect(alignment, config, fractions);
    plotTrainingSizeEffect(trainingSizeResults);
    
elseif choice == 3
    fprintf('\n--- Equivalence Class Comparison ---\n');
    comparisonResults = compareEquivalenceClasses(alignment, config);
    displayComparison(comparisonResults);
    
elseif choice == 4
    fprintf('\n--- Full Analysis ---\n');
    
    % Single experiment
    fprintf('\n1. Running single experiment with equivalence classes...\n');
    results1 = runSingleExperiment(alignment, config);
    displayResults(results1, config);
    
    % Training size effect
    fprintf('\n2. Testing training set size effect...\n');
    fractions = [0.3, 0.5, 0.7, 0.8, 0.9];
    trainingSizeResults = testTrainingSizeEffect(alignment, config, fractions);
    plotTrainingSizeEffect(trainingSizeResults);
    
    % Equivalence comparison
    fprintf('\n3. Comparing equivalence classes...\n');
    comparisonResults = compareEquivalenceClasses(alignment, config);
    displayComparison(comparisonResults);
end

%% Function: Generate Test Alignment with conserved motifs
function alignment = generateTestAlignment(numSeqs, seqLength)
    % Create artificial alignment with conserved patterns
    aminoAcids = 'ACDEFGHIKLMNPQRSTVWY';
    alignment = cell(numSeqs, 1);
    
    % Define conserved motifs at specific positions
    motif1_pos = 5:10;
    motif1 = 'ILMVDE';  % Somewhat conserved
    
    motif2_pos = 20:25;
    motif2 = 'HKRSTY';  % Conserved region
    
    motif3_pos = 35:40;
    motif3 = 'FWYAGP';  % Another conserved region
    
    for i = 1:numSeqs
        seq = '';
        
        for pos = 1:seqLength
            if ismember(pos, motif1_pos)
                % High conservation
                idx = pos - motif1_pos(1) + 1;
                if rand() < 0.8
                    seq = [seq, motif1(idx)];
                else
                    seq = [seq, aminoAcids(randi(length(aminoAcids)))];
                end
                
            elseif ismember(pos, motif2_pos)
                % Very high conservation
                idx = pos - motif2_pos(1) + 1;
                if rand() < 0.95
                    seq = [seq, motif2(idx)];
                else
                    seq = [seq, aminoAcids(randi(length(aminoAcids)))];
                end
                
            elseif ismember(pos, motif3_pos)
                % Moderate conservation
                idx = pos - motif3_pos(1) + 1;
                if rand() < 0.7
                    seq = [seq, motif3(idx)];
                else
                    seq = [seq, aminoAcids(randi(length(aminoAcids)))];
                end
            else
                % Random position
                seq = [seq, aminoAcids(randi(length(aminoAcids)))];
            end
        end
        
        alignment{i} = seq;
    end
end

%% Function: Scramble sequences (preserve composition)
function scrambled = scrambleSequences(sequences)
    scrambled = cell(size(sequences));
    for i = 1:length(sequences)
        seq = sequences{i};
        scrambled{i} = seq(randperm(length(seq)));
    end
end

%% Function: Split data into training and testing
function [trainSeqs, testSeqs] = splitData(sequences, fraction)
    numTrain = round(length(sequences) * fraction);
    indices = randperm(length(sequences));
    
    trainIdx = indices(1:numTrain);
    testIdx = indices(numTrain+1:end);
    
    trainSeqs = sequences(trainIdx);
    testSeqs = sequences(testIdx);
end

%% Function: Convert equivalence classes to regex character class
function charClass = makeCharClass(letters, equivClasses)
    % Find which equivalence class contains these letters
    for i = 1:length(equivClasses)
        if all(ismember(letters, equivClasses{i}))
            % All letters in this class
            charClass = ['[' cell2mat(equivClasses{i}) ']'];
            return;
        end
    end
    
    % If not in single class, use literal characters
    if length(letters) == 1
        charClass = letters{1};
    else
        charClass = ['[' strjoin(letters, '') ']'];
    end
end

%% Function: Generate random regex pattern
function pattern = generateRandomRegex(minLen, maxLen, equivClasses)
    len = randi([minLen, maxLen]);
    pattern = '';
    
    for i = 1:len
        % Randomly select an equivalence class
        classIdx = randi(length(equivClasses));
        selectedClass = equivClasses{classIdx};
        
        if length(selectedClass) == 1
            pattern = [pattern, selectedClass{1}];
        else
            pattern = [pattern, '[' strjoin(selectedClass, '') ']'];
        end
        
        % Occasionally add quantifiers
        if rand() < 0.2
            quantifiers = {'+', '*', '?', '{2,3}'};
            pattern = [pattern, quantifiers{randi(length(quantifiers))}];
        end
    end
end

%% Function: Evaluate regex fitness
function [fitness, precision, recall, f1] = evaluateRegex(pattern, posSeqs, negSeqs)
    try
        % Test on positive sequences (real alignment)
        truePos = 0;
        for i = 1:length(posSeqs)
            if ~isempty(regexp(posSeqs{i}, pattern, 'once'))
                truePos = truePos + 1;
            end
        end
        
        % Test on negative sequences (scrambled)
        falsePos = 0;
        for i = 1:length(negSeqs)
            if ~isempty(regexp(negSeqs{i}, pattern, 'once'))
                falsePos = falsePos + 1;
            end
        end
        
        % Calculate metrics
        trueNeg = length(negSeqs) - falsePos;
        falseNeg = length(posSeqs) - truePos;
        
        if truePos + falsePos > 0
            precision = truePos / (truePos + falsePos);
        else
            precision = 0;
        end
        
        if truePos + falseNeg > 0
            recall = truePos / (truePos + falseNeg);
        else
            recall = 0;
        end
        
        if precision + recall > 0
            f1 = 2 * (precision * recall) / (precision + recall);
        else
            f1 = 0;
        end
        
        % Fitness favors high F1 and shorter patterns
        fitness = f1 - 0.01 * length(pattern);
        
    catch
        % Invalid regex
        fitness = 0;
        precision = 0;
        recall = 0;
        f1 = 0;
    end
end

%% Function: Genetic Algorithm for regex optimization
function [bestPattern, bestFitness, bestMetrics] = optimizeRegex(trainPos, trainNeg, config)
    % Initialize population
    population = cell(config.populationSize, 1);
    for i = 1:config.populationSize
        population{i} = generateRandomRegex(config.minRegexLength, ...
                                            config.maxRegexLength, ...
                                            config.equivClasses);
    end
    
    bestPattern = '';
    bestFitness = -inf;
    bestMetrics = struct('precision', 0, 'recall', 0, 'f1', 0);
    
    for gen = 1:config.numGenerations
        % Evaluate fitness
        fitness = zeros(config.populationSize, 1);
        metrics = cell(config.populationSize, 1);
        
        for i = 1:config.populationSize
            [fitness(i), p, r, f] = evaluateRegex(population{i}, trainPos, trainNeg);
            metrics{i} = struct('precision', p, 'recall', r, 'f1', f);
        end
        
        % Track best
        [maxFit, maxIdx] = max(fitness);
        if maxFit > bestFitness
            bestFitness = maxFit;
            bestPattern = population{maxIdx};
            bestMetrics = metrics{maxIdx};
        end
        
        % Selection (tournament)
        newPopulation = cell(config.populationSize, 1);
        for i = 1:config.populationSize
            % Tournament selection
            idx1 = randi(config.populationSize);
            idx2 = randi(config.populationSize);
            if fitness(idx1) > fitness(idx2)
                newPopulation{i} = population{idx1};
            else
                newPopulation{i} = population{idx2};
            end
        end
        
        % Crossover and mutation
        for i = 1:2:config.populationSize-1
            if rand() < config.crossoverRate
                % Single point crossover
                p1 = newPopulation{i};
                p2 = newPopulation{i+1};
                if length(p1) > 1 && length(p2) > 1
                    point = randi(min(length(p1), length(p2)));
                    newPopulation{i} = [p1(1:point), p2(point+1:end)];
                    newPopulation{i+1} = [p2(1:point), p1(point+1:end)];
                end
            end
            
            % Mutation
            if rand() < config.mutationRate
                newPopulation{i} = generateRandomRegex(config.minRegexLength, ...
                                                       config.maxRegexLength, ...
                                                       config.equivClasses);
            end
        end
        
        population = newPopulation;
        
        if mod(gen, 20) == 0
            fprintf('  Generation %d: Best F1 = %.3f, Pattern = %s\n', ...
                    gen, bestMetrics.f1, bestPattern);
        end
    end
end

%% Function: Cross-validation
function cvResults = crossValidate(alignment, config)
    numSeqs = length(alignment);
    foldSize = floor(numSeqs / config.numFolds);
    
    allPrecisions = zeros(config.numFolds, 1);
    allRecalls = zeros(config.numFolds, 1);
    allF1s = zeros(config.numFolds, 1);
    allPatterns = cell(config.numFolds, 1);
    
    fprintf('Running %d-fold cross-validation...\n', config.numFolds);
    
    for fold = 1:config.numFolds
        fprintf('  Fold %d/%d...', fold, config.numFolds);
        
        % Split into training and validation
        testIdx = ((fold-1)*foldSize + 1):min(fold*foldSize, numSeqs);
        trainIdx = setdiff(1:numSeqs, testIdx);
        
        trainSeqs = alignment(trainIdx);
        valSeqs = alignment(testIdx);
        
        % Generate scrambled negatives
        trainNeg = scrambleSequences(trainSeqs);
        valNeg = scrambleSequences(valSeqs);
        
        % Optimize regex on training fold
        [pattern, ~, ~] = optimizeRegex(trainSeqs, trainNeg, config);
        
        % Evaluate on validation fold
        [~, prec, rec, f1] = evaluateRegex(pattern, valSeqs, valNeg);
        
        allPrecisions(fold) = prec;
        allRecalls(fold) = rec;
        allF1s(fold) = f1;
        allPatterns{fold} = pattern;
        
        fprintf(' F1=%.3f\n', f1);
    end
    
    cvResults = struct();
    cvResults.precision = mean(allPrecisions);
    cvResults.recall = mean(allRecalls);
    cvResults.f1 = mean(allF1s);
    cvResults.precisionStd = std(allPrecisions);
    cvResults.recallStd = std(allRecalls);
    cvResults.f1Std = std(allF1s);
    cvResults.patterns = allPatterns;
end

%% Function: Run single experiment
function results = runSingleExperiment(alignment, config)
    % Split data
    [trainSeqs, testSeqs] = splitData(alignment, config.trainingFraction);
    
    % Generate scrambled sequences
    trainNeg = scrambleSequences(trainSeqs);
    testNeg = scrambleSequences(testSeqs);
    
    fprintf('Training sequences: %d, Test sequences: %d\n', ...
            length(trainSeqs), length(testSeqs));
    
    % Optimize regex
    fprintf('Optimizing regular expression...\n');
    [bestPattern, ~, trainMetrics] = optimizeRegex(trainSeqs, trainNeg, config);
    
    % Test on held-out test set
    [~, testPrec, testRec, testF1] = evaluateRegex(bestPattern, testSeqs, testNeg);
    
    % Cross-validation
    cvResults = crossValidate(alignment, config);
    
    results = struct();
    results.pattern = bestPattern;
    results.patternLength = length(bestPattern);
    results.trainMetrics = trainMetrics;
    results.testPrecision = testPrec;
    results.testRecall = testRec;
    results.testF1 = testF1;
    results.cvResults = cvResults;
end

%% Function: Test effect of training set size
function results = testTrainingSizeEffect(alignment, config, fractions)
    numReplicates = 10;
    numFractions = length(fractions);
    
    avgPrecision = zeros(numFractions, 1);
    avgRecall = zeros(numFractions, 1);
    avgF1 = zeros(numFractions, 1);
    avgLength = zeros(numFractions, 1);
    stdPrecision = zeros(numFractions, 1);
    stdRecall = zeros(numFractions, 1);
    stdF1 = zeros(numFractions, 1);
    
    for f = 1:numFractions
        frac = fractions(f);
        fprintf('\nTesting training fraction %.2f (%d replicates):\n', frac, numReplicates);
        
        precisions = zeros(numReplicates, 1);
        recalls = zeros(numReplicates, 1);
        f1s = zeros(numReplicates, 1);
        lengths = zeros(numReplicates, 1);
        
        for rep = 1:numReplicates
            fprintf('  Replicate %d/%d...', rep, numReplicates);
            
            tempConfig = config;
            tempConfig.trainingFraction = frac;
            tempConfig.numGenerations = 50;  % Faster for multiple runs
            
            result = runSingleExperiment(alignment, tempConfig);
            
            precisions(rep) = result.testPrecision;
            recalls(rep) = result.testRecall;
            f1s(rep) = result.testF1;
            lengths(rep) = result.patternLength;
            
            fprintf(' F1=%.3f\n', result.testF1);
        end
        
        avgPrecision(f) = mean(precisions);
        avgRecall(f) = mean(recalls);
        avgF1(f) = mean(f1s);
        avgLength(f) = mean(lengths);
        stdPrecision(f) = std(precisions);
        stdRecall(f) = std(recalls);
        stdF1(f) = std(f1s);
    end
    
    results = struct();
    results.fractions = fractions;
    results.precision = avgPrecision;
    results.recall = avgRecall;
    results.f1 = avgF1;
    results.length = avgLength;
    results.precisionStd = stdPrecision;
    results.recallStd = stdRecall;
    results.f1Std = stdF1;
end

%% Function: Compare equivalence classes vs no equivalence
function results = compareEquivalenceClasses(alignment, config)
    numReplicates = 10;
    
    fprintf('\nTesting WITH equivalence classes...\n');
    withEquiv = struct();
    withEquiv.precisions = zeros(numReplicates, 1);
    withEquiv.recalls = zeros(numReplicates, 1);
    withEquiv.f1s = zeros(numReplicates, 1);
    withEquiv.lengths = zeros(numReplicates, 1);
    
    for rep = 1:numReplicates
        fprintf('  Replicate %d/%d...', rep, numReplicates);
        result = runSingleExperiment(alignment, config);
        withEquiv.precisions(rep) = result.testPrecision;
        withEquiv.recalls(rep) = result.testRecall;
        withEquiv.f1s(rep) = result.testF1;
        withEquiv.lengths(rep) = result.patternLength;
        fprintf(' F1=%.3f\n', result.testF1);
    end
    
    fprintf('\nTesting WITHOUT equivalence classes (individual amino acids)...\n');
    noEquivConfig = config;
    noEquivConfig.equivClasses = config.noEquivClasses;
    
    noEquiv = struct();
    noEquiv.precisions = zeros(numReplicates, 1);
    noEquiv.recalls = zeros(numReplicates, 1);
    noEquiv.f1s = zeros(numReplicates, 1);
    noEquiv.lengths = zeros(numReplicates, 1);
    
    for rep = 1:numReplicates
        fprintf('  Replicate %d/%d...', rep, numReplicates);
        result = runSingleExperiment(alignment, noEquivConfig);
        noEquiv.precisions(rep) = result.testPrecision;
        noEquiv.recalls(rep) = result.testRecall;
        noEquiv.f1s(rep) = result.testF1;
        noEquiv.lengths(rep) = result.patternLength;
        fprintf(' F1=%.3f\n', result.testF1);
    end
    
    results = struct();
    results.withEquiv = withEquiv;
    results.noEquiv = noEquiv;
end

%% Function: Display results
function displayResults(results, config)
    fprintf('\n=== RESULTS ===\n');
    fprintf('Best Regular Expression: %s\n', results.pattern);
    fprintf('Pattern Length: %d characters\n', results.patternLength);
    fprintf('\nTraining Set Performance:\n');
    fprintf('  Precision: %.3f\n', results.trainMetrics.precision);
    fprintf('  Recall:    %.3f\n', results.trainMetrics.recall);
    fprintf('  F1 Score:  %.3f\n', results.trainMetrics.f1);
    fprintf('\nTest Set Performance:\n');
    fprintf('  Precision: %.3f\n', results.testPrecision);
    fprintf('  Recall:    %.3f\n', results.testRecall);
    fprintf('  F1 Score:  %.3f\n', results.testF1);
    fprintf('\nCross-Validation Results (%d folds):\n', config.numFolds);
    fprintf('  Precision: %.3f ± %.3f\n', results.cvResults.precision, results.cvResults.precisionStd);
    fprintf('  Recall:    %.3f ± %.3f\n', results.cvResults.recall, results.cvResults.recallStd);
    fprintf('  F1 Score:  %.3f ± %.3f\n', results.cvResults.f1, results.cvResults.f1Std);
end

%% Function: Plot training size effect
function plotTrainingSizeEffect(results)
    figure('Position', [100, 100, 1400, 500]);
    
    subplot(1, 3, 1);
    errorbar(results.fractions, results.precision, results.precisionStd, 'o-', 'LineWidth', 2);
    hold on;
    errorbar(results.fractions, results.recall, results.recallStd, 's-', 'LineWidth', 2);
    errorbar(results.fractions, results.f1, results.f1Std, '^-', 'LineWidth', 2);
    xlabel('Training Set Fraction');
    ylabel('Score');
    title('Performance vs Training Set Size');
    legend('Precision', 'Recall', 'F1', 'Location', 'best');
    grid on;
    ylim([0, 1.1]);
    
    subplot(1, 3, 2);
    plot(results.fractions, results.length, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Training Set Fraction');
    ylabel('Pattern Length (characters)');
    title('Regex Length vs Training Set Size');
    grid on;
    
    subplot(1, 3, 3);
    bar(results.fractions, [results.precision, results.recall, results.f1]);
    xlabel('Training Set Fraction');
    ylabel('Score');
    title('Performance Metrics by Training Size');
    legend('Precision', 'Recall', 'F1', 'Location', 'best');
    grid on;
    ylim([0, 1.1]);
    
    sgtitle('Effect of Training Set Size on Performance (10 replicates per condition)');
    
    % Print summary
    fprintf('\n=== TRAINING SET SIZE ANALYSIS ===\n');
    fprintf('Effect of training fraction on performance:\n');
    for i = 1:length(results.fractions)
        fprintf('  %.0f%% training: F1=%.3f±%.3f, Length=%d\n', ...
                results.fractions(i)*100, results.f1(i), results.f1Std(i), ...
                round(results.length(i)));
    end
end

%% Function: Display equivalence class comparison
function displayComparison(results)
    figure('Position', [100, 100, 1000, 600]);
    
    % Compare metrics
    subplot(2, 2, 1);
    withData = [mean(results.withEquiv.precisions), mean(results.withEquiv.recalls), mean(results.withEquiv.f1s)];
    noData = [mean(results.noEquiv.precisions), mean(results.noEquiv.recalls), mean(results.noEquiv.f1s)];
    
    b = bar([withData; noData]);
    set(gca, 'XTickLabel', {'With Equivalence', 'No Equivalence'});
    ylabel('Score');
    title('Performance Comparison');
    legend('Precision', 'Recall', 'F1', 'Location', 'best');
    grid on;
    ylim([0, 1.1]);
    
    % Box plots
    subplot(2, 2, 2);
    boxplot([results.withEquiv.f1s; results.noEquiv.f1s], ...
            [ones(size(results.withEquiv.f1s)); 2*ones(size(results.noEquiv.f1s))], ...
            'Labels', {'With Equiv', 'No Equiv'});
    ylabel('F1 Score');
    title('F1 Score Distribution');
    grid on;
    
    subplot(2, 2, 3);
    boxplot([results.withEquiv.lengths; results.noEquiv.lengths], ...
            [ones(size(results.withEquiv.lengths)); 2*ones(size(results.noEquiv.lengths))], ...
            'Labels', {'With Equiv', 'No Equiv'});
    ylabel('Pattern Length');
    title('Regex Length Distribution');
    grid on;
    
    subplot(2, 2, 4);
    data = [mean(results.withEquiv.f1s), mean(results.withEquiv.lengths);
            mean(results.noEquiv.f1s), mean(results.noEquiv.lengths)];
    bar(data);
    set(gca, 'XTickLabel', {'With Equivalence', 'No Equivalence'});
    legend('F1 Score', 'Pattern Length', 'Location', 'best');
    title('Summary Statistics');
    grid on;
    
    sgtitle('Equivalence Classes vs Individual Amino Acids (10 replicates each)');
    
    % Statistical summary
    fprintf('\n=== EQUIVALENCE CLASS COMPARISON ===\n');
    fprintf('\nWith Equivalence Classes:\n');
    fprintf('  Precision: %.3f ± %.3f\n', mean(results.withEquiv.precisions), std(results.withEquiv.precisions));
    fprintf('  Recall:    %.3f ± %.3f\n', mean(results.withEquiv.recalls), std(results.withEquiv.recalls));
    fprintf('  F1 Score:  %.3f ± %.3f\n', mean(results.withEquiv.f1s), std(results.withEquiv.f1s));
    fprintf('  Avg Length: %.1f ± %.1f\n', mean(results.withEquiv.lengths), std(results.withEquiv.lengths));
    
    fprintf('\nWithout Equivalence Classes:\n');
    fprintf('  Precision: %.3f ± %.3f\n', mean(results.noEquiv.precisions), std(results.noEquiv.precisions));
    fprintf('  Recall:    %.3f ± %.3f\n', mean(results.noEquiv.recalls), std(results.noEquiv.recalls));
    fprintf('  F1 Score:  %.3f ± %.3f\n', mean(results.noEquiv.f1s), std(results.noEquiv.f1s));
    fprintf('  Avg Length: %.1f ± %.1f\n', mean(results.noEquiv.lengths), std(results.noEquiv.lengths));
    
    % Statistical test
    [~, pValue] = ttest2(results.withEquiv.f1s, results.noEquiv.f1s);
    fprintf('\nStatistical Test (t-test on F1 scores):\n');
    fprintf('  p-value: %.4f\n', pValue);
    if pValue < 0.05
        fprintf('  Result: SIGNIFICANT difference (p < 0.05)\n');
    else
        fprintf('  Result: No significant difference (p >= 0.05)\n');
    end
    
    fprintf('\n=== INTERPRETATION ===\n');
    fprintf('Equivalence classes [ILMV], [DENQ], [HKR], [ST], [AG], [FWY], C, P:\n');
    if mean(results.withEquiv.f1s) > mean(results.noEquiv.f1s)
        fprintf('  ✓ BETTER performance (%.1f%% improvement in F1)\n', ...
                100*(mean(results.withEquiv.f1s) - mean(results.noEquiv.f1s))/mean(results.noEquiv.f1s));
    else
        fprintf('  ✗ WORSE performance\n');
    end
    
    if mean(results.withEquiv.lengths) < mean(results.noEquiv.lengths)
        fprintf('  ✓ SHORTER patterns (%.1f%% reduction)\n', ...
                100*(mean(results.noEquiv.lengths) - mean(results.withEquiv.lengths))/mean(results.noEquiv.lengths));
    else
        fprintf('  ✗ LONGER patterns\n');
    end
end