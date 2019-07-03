
function [N,n,x,n_obs,n_para,lPl] = read_groops_normals( file_path )

tmp = importdata([file_path '-normals.txt']);
us = tmp.data;
for i=1:size(us,1)
    us(i,:) = circshift(us(i,:),[0 i-1]);
end
us(isnan(us)) = 0;
N = us + triu(us, 1)';

tmp = importdata([file_path '-normals.rightHandSide.txt']);
n = tmp.data;

try
    tmp = importdata([file_path '-x.txt']);
    x = tmp.data;
catch
    x = nan(size(n));
end

tmp = xmlread( [file_path '-normals.info.xml']);
element = tmp.getElementsByTagName('observationCount');
n_obs = str2double(element.item(0).getTextContent);
element = tmp.getElementsByTagName('parameterCount');
n_para = str2double(element.item(0).getTextContent);
element = tmp.getElementsByTagName('lPl');
lPl = str2double(element.item(0).getTextContent);

end