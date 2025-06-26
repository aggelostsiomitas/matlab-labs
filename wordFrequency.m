function wordCount = wordFrequency(text)
text = lower(text);  
    text = regexprep(text, '[\d.,;:!?()"-]', '');  
    
    words = strsplit(text);
    wordCount = {};  
    
    for i = 1:length(words)
        word = words{i}; 
        if isempty(word)
            continue;
        end
        idx = find(strcmp(wordCount(:, 1), word)); 
        if isempty(idx)  
            wordCount = [wordCount; {word, 1}];  
        else
            wordCount{idx, 2} = wordCount{idx, 2} + 1;  
        end
    end
end