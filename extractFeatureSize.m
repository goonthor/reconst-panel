function fsize = extractFeatureSize(handles)                              % 
fsizex1 = str2double(get(handles.featuresizex1,'String'));                %
fsizex2 = str2double(get(handles.featuresizex2,'String'));                %
fsizex3 = str2double(get(handles.featuresizex3,'String'));                %
fsizex4 = str2double(get(handles.featuresizex4,'String'));                %
fsize = [fsizex1 fsizex2 fsizex3 fsizex4];                                %
fsize(isnan(fsize))=[];
if numel(fsize) ~= handles.pntdim
    fsize=zeros(handles.pntdim,1);
end