function newvar=loadname(fullfilename,name)

var=load(fullfilename);
structname=who('-file',fullfilename); %WHO -FILE FILENAME lists the variables in the specified .MAT file
aux=sprintf('var.%s',structname{1});  % STR = SPRINTF(FORMAT, A, ...) applies the FORMAT to all elements of
                                      %   array A and any additional array arguments in column order, and returns
                                      %   the results to string STR.
newvar=eval(aux);                     %where s is a string, causes MATLAB to execute
                                      %   the string as an expression or statement.