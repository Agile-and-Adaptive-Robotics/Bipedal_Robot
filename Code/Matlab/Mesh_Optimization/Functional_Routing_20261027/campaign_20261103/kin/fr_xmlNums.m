function v = fr_xmlNums(node)
%FR_XMLNUMS Parse a whitespace-separated number list from XML text.
parts = strsplit(strtrim(char(node.getTextContent())));
v = str2double(parts);
end
