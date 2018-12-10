% conversion of hd EEG from netstation into mat file
% initial txt file comes from the Netstation acquisition,
% only first lines were removed manually as well as last lines
% informations saved are:
% Electrod 'Name', 'type' and 3D coordinates ('xcoord','ycoord','zcoord')

filename = 'C:\Users\Doro\Github\EEGtoolswc\branches\brancheDC\FASST\CSG\coordinates256.txt';
fileID = fopen(filename);
C = textscan(fileID,'%s');
fclose(fileID);
C = strrep(C{1}, '<sensor>', '');
C = strrep(C, '<', '');
C = strrep(C, '>', ' ');
C = strrep(C,'name /name ','');
C = strrep(C,'/number ','');
C = strrep(C,'/sensor ','');
C = strrep(C,'number','E');
C = strrep(C,'/x ','');
C = strrep(C,'/y ','');
C = strrep(C,'/z ','');
C = strrep(C,'/type ','');

ik = 1;
data = cell(1,2);
for iC = 1 : max(size(C))
    if ~strcmp(C{iC},'')
        data = textscan(C{iC},'%s %f');
        switch char(data{1})
            case 'E'
                Name{ik} = [char(data{1}) num2str(data{2})];
            case 'type'
                types(ik) = data{2};
            case 'x'
                xcoord(ik) = data{2};
            case 'y'
                ycoord(ik) = data{2};
            case 'z'
                zcoord(ik) = data{2};
                ik = ik+1;
        end

    end
end

save('C:\Users\Doro\Github\EEGtoolswc\branches\brancheDC\FASST\CSG\NetStationMap256.mat','Name','types','xcoord','ycoord','zcoord')