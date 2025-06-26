clear
clc;
t1=tic;
lines = upper(readlines("text.txt"));
fullText = strjoin(lines, ' ');
 
 
sentences = strtrim(split(fullText, '.'));
sentences(sentences == "") = [];  
cleanText = regexprep(fullText, '[\d.,;:!?()"-]', '');  
cleanText = lower(cleanText); 
 
 
stopwords = ["a", "of", "on", "i", "for", "with", "the", "at", "from", "in", "to"];
 
words = split(cleanText);
words = words(~ismember(words, stopwords));
words(words == "") = [];  
 
wordCountMap = containers.Map();
for i = 1:length(words)
    word = words{i}; 
    if isempty(word)
        continue;
    end
    if isKey(wordCountMap, word)  
        wordCountMap(word) = wordCountMap(word) + 1;  
    else
        wordCountMap(word) = 1; 
    end
end
 
 
sentenceScores = zeros(1, length(sentences));
 
for i = 1:length(sentences)
 
    sentence = lower(sentences{i});
    sentenceWords = split(sentence);
    sentenceWords = sentenceWords(~ismember(sentenceWords, stopwords)); 
    sentenceWords(sentenceWords == "") = [];  
    
   
    score = 0;
    for j = 1:length(sentenceWords)
        word = sentenceWords{j};
        if isKey(wordCountMap, word) 
            score = score + wordCountMap(word); 
        end
    end
 
    sentenceScores(i) = score;  
end
 
[sortedScores, sortedIndices] = sort(sentenceScores, 'descend');
numSentences = input('Πόσες προτάσεις θέλεις να περιλαμβάνει η περίληψη: ');
numSentences = min(numSentences, length(sentences));
topSentences = sentences(sortedIndices(1:numSentences));
summary = strjoin(topSentences, '. ');
disp('Η περίληψη του κειμένου είναι:');
disp(summary);
 t2=toc(t1);
 disp("the time needed is :" +t2);
wordFrequency = cell2mat(values(wordCountMap)); 
figure;
hist(wordFrequency, 10);  
xlabel('Συχνότητα Λέξεων');
ylabel('Αριθμός Λέξεων');
title('Ιστόγραμμα Συχνοτήτων Λέξεων');
grid on;