%% DNA Assembly Simulator with De Bruijn Graph
% Main script for genome assembly simulation and analysis

clear all; close all; clc;

%% Parameters
params = struct();
params.genomeLength = 500;
params.baseFreqs = [0.25, 0.25, 0.25, 0.25]; % A, T, G, C
params.repeatLength = 20;
params.numRepeats = 3;
params.readLength = 50;
params.coverage = 10;
params.errorRate = 0.01;
params.kmerLength = 15;

%% Menu Selection
fprintf('\n=== DNA Assembly Simulator ===\n');
fprintf('Option 1. Run single assembly\n');
fprintf('Option 2. Run statistical analysis\n');
fprintf('Option 3. Run both\n');
choice = input('Select option (1-3): ');

if choice == 1 || choice == 3
    fprintf('\n--- Running Single Assembly ---\n');
    results = runSingleAssembly(params);
    displayResults(results, params);
    visualizeGraph(results);
end

if choice == 2 || choice == 3
    fprintf('\n--- Running Statistical Analysis ---\n');
    stats = runStatisticalAnalysis(params);
    plotStatistics(stats, params);
end

%% Function: Run Single Assembly
function results = runSingleAssembly(params)
    % Generate genome
    genome = generateDNA(params.genomeLength, params.baseFreqs);
    genome = insertRepeats(genome, params.repeatLength, params.numRepeats);
    
    % Generate reads
    reads = generateReads(genome, params.readLength, params.coverage, params.errorRate);
    
    % Extract k-mers
    kmers = extractKmers(reads, params.kmerLength);
    
    % Build De Bruijn graph
    [graph, edges] = buildDeBruijnGraph(kmers);
    
    % Assemble contigs
    contigs = assembleContigs(graph);
    
    % Calculate Lander-Waterman prediction
    predicted = landerWaterman(params.genomeLength, params.readLength, params.coverage);
    
    % Package results
    results = struct();
    results.genome = genome;
    results.numReads = length(reads);
    results.numKmers = length(kmers);
    results.contigs = contigs;
    results.numContigs = length(contigs);
    results.predicted = predicted;
    results.graph = graph;
    results.edges = edges;
end

%% Function: Generate Random DNA Sequence
function sequence = generateDNA(length, baseFreqs)
    bases = 'ATGC';
    cumFreqs = cumsum(baseFreqs);
    
    sequence = '';
    for i = 1:length
        r = rand();
        idx = find(r < cumFreqs, 1, 'first');
        sequence = [sequence, bases(idx)];
    end
end

%% Function: Insert Repetitive Sequences
function sequence = insertRepeats(sequence, repeatLength, numRepeats)
    if numRepeats < 2 || length(sequence) < repeatLength * numRepeats
        return;
    end
    
    % Generate a random repeat unit
    repeat = generateDNA(repeatLength, [0.25, 0.25, 0.25, 0.25]);
    repeatSeq = repmat(repeat, 1, numRepeats);
    
    % Insert at random position
    maxPos = length(sequence) - length(repeatSeq);
    if maxPos > 0
        insertPos = randi(maxPos);
        sequence = [sequence(1:insertPos-1), repeatSeq, ...
                   sequence(insertPos+length(repeatSeq):end)];
    end
end

%% Function: Generate Sequencing Reads
function reads = generateReads(sequence, readLength, coverage, errorRate)
    numReads = ceil((length(sequence) * coverage) / readLength);
    reads = cell(numReads, 1);
    bases = 'ATGC';
    
    for i = 1:numReads
        % Random starting position
        maxStart = length(sequence) - readLength + 1;
        if maxStart < 1
            continue;
        end
        startPos = randi(maxStart);
        read = sequence(startPos:startPos+readLength-1);
        
        % Add sequencing errors
        for j = 1:length(read)
            if rand() < errorRate
                % Replace with different base
                currentBase = read(j);
                otherBases = bases(bases ~= currentBase);
                read(j) = otherBases(randi(3));
            end
        end
        
        reads{i} = read;
    end
end

%% Function: Extract K-mers from Reads
function kmers = extractKmers(reads, k)
    kmerSet = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    
    for i = 1:length(reads)
        read = reads{i};
        for j = 1:length(read)-k+1
            kmer = read(j:j+k-1);
            kmerSet(kmer) = true;
        end
    end
    
    kmers = keys(kmerSet);
end

%% Function: Build De Bruijn Graph
function [graph, edges] = buildDeBruijnGraph(kmers)
    graph = containers.Map('KeyType', 'char', 'ValueType', 'any');
    edges = struct('from', {}, 'to', {}, 'label', {});
    edgeIdx = 1;
    
    for i = 1:length(kmers)
        kmer = kmers{i};
        prefix = kmer(1:end-1);
        suffix = kmer(2:end);
        
        % Add edge to graph
        if isKey(graph, prefix)
            neighbors = graph(prefix);
            neighbors{end+1} = suffix;
            graph(prefix) = neighbors;
        else
            graph(prefix) = {suffix};
        end
        
        % Store edge information
        edges(edgeIdx).from = prefix;
        edges(edgeIdx).to = suffix;
        edges(edgeIdx).label = kmer;
        edgeIdx = edgeIdx + 1;
    end
end

%% Function: Assemble Contigs from De Bruijn Graph
function contigs = assembleContigs(graph)
    contigs = {};
    visited = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    
    % Calculate in-degrees and out-degrees
    [inDegree, outDegree] = calculateDegrees(graph);
    
    % Find starting nodes
    allNodes = keys(graph);
    startNodes = {};
    
    for i = 1:length(allNodes)
        node = allNodes{i};
        inDeg = 0;
        if isKey(inDegree, node)
            inDeg = inDegree(node);
        end
        
        outDeg = length(graph(node));
        
        % Start from nodes with special properties
        if outDeg > 0 && (inDeg == 0 || inDeg ~= 1 || outDeg ~= 1)
            startNodes{end+1} = node;
        end
    end
    
    % Also include unvisited nodes
    for i = 1:length(allNodes)
        node = allNodes{i};
        if ~ismember(node, startNodes)
            startNodes{end+1} = node;
        end
    end
    
    % Traverse from each starting node
    for i = 1:length(startNodes)
        start = startNodes{i};
        if isKey(visited, start)
            continue;
        end
        
        current = start;
        path = {current};
        visited(current) = true;
        
        % Follow linear path
        while isKey(graph, current)
            neighbors = graph(current);
            
            % Break at branch points
            if length(neighbors) > 1
                break;
            end
            
            next = neighbors{1};
            
            % Break if next has multiple incoming edges
            if isKey(inDegree, next) && inDegree(next) > 1
                break;
            end
            
            % Break at cycles
            if isKey(visited, next)
                break;
            end
            
            path{end+1} = next;
            visited(next) = true;
            current = next;
        end
        
        % Reconstruct sequence from path
        if length(path) > 0
            contigSeq = path{1};
            for j = 2:length(path)
                contigSeq = [contigSeq, path{j}(end)];
            end
            
            contigs{end+1} = struct('sequence', contigSeq, 'path', {path});
        end
    end
end

%% Function: Calculate Node Degrees
function [inDegree, outDegree] = calculateDegrees(graph)
    inDegree = containers.Map('KeyType', 'char', 'ValueType', 'double');
    outDegree = containers.Map('KeyType', 'char', 'ValueType', 'double');
    
    allNodes = keys(graph);
    
    for i = 1:length(allNodes)
        node = allNodes{i};
        neighbors = graph(node);
        outDegree(node) = length(neighbors);
        
        for j = 1:length(neighbors)
            neighbor = neighbors{j};
            if isKey(inDegree, neighbor)
                inDegree(neighbor) = inDegree(neighbor) + 1;
            else
                inDegree(neighbor) = 1;
            end
        end
    end
end

%% Function: Lander-Waterman Equation
function predicted = landerWaterman(genomeLength, readLength, coverage)
    lambda = coverage;
    predicted = max(1, round((genomeLength / readLength) * exp(-lambda)));
end

%% Function: Display Results
function displayResults(results, params)
    fprintf('\n=== Assembly Results ===\n');
    fprintf('Genome Length: %d bp\n', length(results.genome));
    fprintf('Reads Generated: %d\n', results.numReads);
    fprintf('Unique K-mers: %d\n', results.numKmers);
    fprintf('Contigs Assembled: %d\n', results.numContigs);
    fprintf('Lander-Waterman Predicted: %d\n', results.predicted);
    fprintf('Prediction Error: %d\n', abs(results.numContigs - results.predicted));
    
    fprintf('\n=== Contig Details ===\n');
    for i = 1:min(10, length(results.contigs))
        contig = results.contigs{i};
        fprintf('Contig %d: %d bp - %s...\n', i, length(contig.sequence), ...
                contig.sequence(1:min(50, length(contig.sequence))));
    end
    
    if length(results.contigs) > 10
        fprintf('... and %d more contigs\n', length(results.contigs) - 10);
    end
end

%% Function: Visualize De Bruijn Graph
function visualizeGraph(results)
    figure('Position', [100, 100, 1200, 600]);
    
    % Create subplot for graph visualization
    subplot(1, 2, 1);
    
    % Get unique nodes (limit for visualization)
    allNodes = {};
    for i = 1:min(100, length(results.edges))
        allNodes{end+1} = results.edges(i).from;
        allNodes{end+1} = results.edges(i).to;
    end
    uniqueNodes = unique(allNodes);
    
    % Create positions in circular layout
    numNodes = length(uniqueNodes);
    theta = linspace(0, 2*pi, numNodes+1);
    theta = theta(1:end-1);
    radius = 1;
    x = radius * cos(theta);
    y = radius * sin(theta);
    
    nodePos = containers.Map(uniqueNodes, num2cell([x; y]', 2));
    
    hold on;
    
    % Draw edges
    for i = 1:min(200, length(results.edges))
        edge = results.edges(i);
        if isKey(nodePos, edge.from) && isKey(nodePos, edge.to)
            fromPos = nodePos(edge.from);
            toPos = nodePos(edge.to);
            plot([fromPos(1), toPos(1)], [fromPos(2), toPos(2)], ...
                 'Color', [0.6, 0.8, 1], 'LineWidth', 0.5);
        end
    end
    
    % Draw contig paths
    colors = lines(min(10, length(results.contigs)));
    for i = 1:min(10, length(results.contigs))
        contig = results.contigs{i};
        path = contig.path;
        
        for j = 1:length(path)-1
            if isKey(nodePos, path{j}) && isKey(nodePos, path{j+1})
                fromPos = nodePos(path{j});
                toPos = nodePos(path{j+1});
                plot([fromPos(1), toPos(1)], [fromPos(2), toPos(2)], ...
                     'Color', colors(i, :), 'LineWidth', 2);
            end
        end
    end
    
    % Draw nodes
    scatter(x, y, 50, [0.2, 0.4, 0.8], 'filled');
    
    axis equal;
    axis off;
    title('De Bruijn Graph Structure');
    legend('Edges', 'Contig Paths', 'Nodes', 'Location', 'best');
    
    % Contig length distribution
    subplot(1, 2, 2);
    contigLengths = zeros(1, length(results.contigs));
    for i = 1:length(results.contigs)
        contigLengths(i) = length(results.contigs{i}.sequence);
    end
    
    histogram(contigLengths, 20, 'FaceColor', [0.3, 0.6, 0.9]);
    xlabel('Contig Length (bp)');
    ylabel('Frequency');
    title('Contig Length Distribution');
    grid on;
    
    sgtitle(sprintf('Assembly Visualization: %d Contigs (LW Predicted: %d)', ...
            results.numContigs, results.predicted));
end

%% Function: Run Statistical Analysis
function stats = runStatisticalAnalysis(params)
    numRuns = 20;
    coverages = [5, 10, 15, 20, 25];
    errorRates = [0, 0.01, 0.05, 0.1];
    
    % Coverage analysis
    fprintf('Running coverage analysis...\n');
    coverageData = zeros(length(coverages), 4); % cov, observed, predicted, stdDev
    
    for c = 1:length(coverages)
        cov = coverages(c);
        fprintf('  Coverage %dx: ', cov);
        
        contigCounts = zeros(numRuns, 1);
        for run = 1:numRuns
            if mod(run, 5) == 0
                fprintf('.');
            end
            
            tempParams = params;
            tempParams.coverage = cov;
            result = runSingleAssembly(tempParams);
            contigCounts(run) = result.numContigs;
        end
        
        coverageData(c, 1) = cov;
        coverageData(c, 2) = mean(contigCounts);
        coverageData(c, 3) = landerWaterman(params.genomeLength, params.readLength, cov);
        coverageData(c, 4) = std(contigCounts);
        
        fprintf(' Done (avg: %.1f contigs)\n', coverageData(c, 2));
    end
    
    % Error rate analysis
    fprintf('Running error rate analysis...\n');
    errorData = zeros(length(errorRates), 4); % error, observed, predicted, ratio
    
    for e = 1:length(errorRates)
        err = errorRates(e);
        fprintf('  Error rate %.2f%%: ', err * 100);
        
        contigCounts = zeros(numRuns, 1);
        for run = 1:numRuns
            if mod(run, 5) == 0
                fprintf('.');
            end
            
            tempParams = params;
            tempParams.errorRate = err;
            result = runSingleAssembly(tempParams);
            contigCounts(run) = result.numContigs;
        end
        
        errorData(e, 1) = err;
        errorData(e, 2) = mean(contigCounts);
        errorData(e, 3) = landerWaterman(params.genomeLength, params.readLength, params.coverage);
        errorData(e, 4) = errorData(e, 2) / errorData(e, 3);
        
        fprintf(' Done (avg: %.1f contigs, ratio: %.2f)\n', ...
                errorData(e, 2), errorData(e, 4));
    end
    
    stats = struct();
    stats.coverageData = coverageData;
    stats.errorData = errorData;
end

%% Function: Plot Statistical Results
function plotStatistics(stats, params)
    figure('Position', [100, 100, 1400, 500]);
    
    % Coverage analysis
    subplot(1, 3, 1);
    coverageData = stats.coverageData;
    errorbar(coverageData(:, 1), coverageData(:, 2), coverageData(:, 4), ...
             'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Observed');
    hold on;
    plot(coverageData(:, 1), coverageData(:, 3), 's--', ...
         'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'LW Predicted');
    xlabel('Coverage');
    ylabel('Number of Contigs');
    title('Coverage vs. Contig Count');
    legend('Location', 'best');
    grid on;
    
    % Error rate analysis - bar chart
    subplot(1, 3, 2);
    errorData = stats.errorData;
    x = 1:size(errorData, 1);
    bar(x, [errorData(:, 2), errorData(:, 3)]);
    set(gca, 'XTickLabel', arrayfun(@(e) sprintf('%.0f%%', e*100), ...
        errorData(:, 1), 'UniformOutput', false));
    xlabel('Error Rate');
    ylabel('Number of Contigs');
    title('Error Rate Impact on Assembly');
    legend('Observed', 'LW Predicted', 'Location', 'best');
    grid on;
    
    % Prediction accuracy
    subplot(1, 3, 3);
    plot(errorData(:, 1) * 100, errorData(:, 4), 'o-', ...
         'LineWidth', 2, 'MarkerSize', 10, 'Color', [0.8, 0.2, 0.2]);
    hold on;
    yline(1, '--k', 'LineWidth', 1.5);
    xlabel('Error Rate (%)');
    ylabel('Observed / Predicted Ratio');
    title('LW Prediction Accuracy vs. Error Rate');
    grid on;
    ylim([0.5, max(errorData(:, 4)) * 1.2]);
    
    sgtitle('Statistical Analysis Results (20 runs per condition)');
    
    % Print analysis summary
    fprintf('\n=== ANALYSIS SUMMARY ===\n\n');
    fprintf('1. LANDER-WATERMAN EQUATION PERFORMANCE:\n');
    fprintf('   The LW equation predicts contig counts based on:\n');
    fprintf('   - Genome length, read length, and coverage\n');
    fprintf('   - Assumes random sequencing without errors\n');
    fprintf('   - Works best at low coverage where gaps dominate\n\n');
    
    fprintf('2. COVERAGE EFFECTS:\n');
    fprintf('   As coverage increases:\n');
    fprintf('   - Fewer contigs (better assembly)\n');
    fprintf('   - LW prediction becomes more accurate\n');
    fprintf('   - At high coverage, repeats become limiting factor\n\n');
    
    fprintf('3. ERROR RATE IMPACT:\n');
    fprintf('   Observed/Predicted ratios at different error rates:\n');
    for i = 1:size(errorData, 1)
        if errorData(i, 4) > 1
            comparison = 'more';
        else
            comparison = 'fewer';
        end
        fprintf('   %.0f%% error: %.2f (%.0f%% %s than predicted)\n', ...
                errorData(i, 1) * 100, errorData(i, 4), ...
                abs(1 - errorData(i, 4)) * 100, comparison);
    end
    fprintf('\n   Higher errors create false k-mers → graph fragmentation\n');
    fprintf('   → More contigs than LW predicts\n\n');
    
    fprintf('4. REPETITIVE SEQUENCES:\n');
    fprintf('   Tandem repeats cause:\n');
    fprintf('   - Ambiguous branches in De Bruijn graph\n');
    fprintf('   - Contig breaks at ambiguous junctions\n');
    fprintf('   - More contigs than predicted for random sequences\n\n');
    
    fprintf('5. CONCLUSION:\n');
    fprintf('   LW equation provides good baseline but:\n');
    fprintf('   - Underestimates with sequencing errors\n');
    fprintf('   - Underestimates with repetitive sequences\n');
    fprintf('   - Most accurate for high-quality, random sequences\n');
end