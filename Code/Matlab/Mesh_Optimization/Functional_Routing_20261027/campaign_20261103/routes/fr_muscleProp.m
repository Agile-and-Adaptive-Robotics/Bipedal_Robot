function v = fr_muscleProp(osimPath, muscleName, propTag)
%FR_MUSCLEPROP Read a scalar property (e.g. max_isometric_force) of one
% muscle straight from the .osim XML.
doc = xmlread(osimPath);
root = doc.getDocumentElement();
muscles = root.getElementsByTagName('Thelen2003Muscle');
for i = 0:muscles.getLength()-1
    m = muscles.item(i);
    if strcmp(char(m.getAttribute('name')), muscleName)
        el = m.getElementsByTagName(propTag);
        if el.getLength() == 0
            error('fr_muscleProp:NoProp', ...
                'muscle %s has no %s in %s', muscleName, propTag, osimPath);
        end
        v = str2double(strtrim(char(el.item(0).getTextContent())));
        return
    end
end
error('fr_muscleProp:NoMuscle', 'muscle %s not found', muscleName);
end
