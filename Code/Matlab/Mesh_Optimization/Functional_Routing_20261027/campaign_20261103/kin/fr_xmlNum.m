function v = fr_xmlNum(node)
%FR_XMLNUM Parse a single number from an XML element's text content.
v = str2double(strtrim(char(node.getTextContent())));
end
