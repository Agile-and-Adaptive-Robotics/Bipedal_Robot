function n = fr_xmlChild(node, tag)
%FR_XMLCHILD First direct ELEMENT child of node with the given tag.
n = [];
if isempty(node); return; end
kids = node.getChildNodes();
for i = 0:kids.getLength()-1
    c = kids.item(i);
    if c.getNodeType() == c.ELEMENT_NODE && strcmp(char(c.getNodeName()), tag)
        n = c;
        return
    end
end
end
