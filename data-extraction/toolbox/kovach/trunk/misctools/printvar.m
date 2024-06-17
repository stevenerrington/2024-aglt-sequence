function printvar(var,fid,newline,depth)

%  printvar(var,fid)
%
% Prints the contents of variable 'var' to a file or to std out if fid is
% unspecified.

if nargin < 3
    newline = true;
end
if nargin < 4
    depth = 0;
end
if nargin < 2 || isempty(fid)
    fid = 1;
end

switch class(var)
   case 'char',
       fprintf(fid,'''%s''',var);
   case 'cell'
       fprintf(fid,'{');

       for j = 1:length(var)
           printvar(var{j},fid,false,depth+1)
           if j < length(var)
               fprintf(fid,',');
           end
       end
       fprintf(fid,'}');
   case 'struct' 
        if depth>0
            fprintf(fid,['...\n',repmat('\t',1,depth)]);
        end
        fprintf(fid,'struct( ');
        fldn= fieldnames(var);
        for k = 1:length(fldn)
            fprintf(fid,'''%s'' , ',fldn{k});
            %if length(var)>1
                printvar({var.(fldn{k})},fid,false,depth+1)
%             else
%                 printvar(fid,var.(fldn{k}),false,depth+2)
%             end
            if k < length(fldn)
               fprintf(fid,[',...\n',repmat('\t',1,depth+2)]);
            end
        end
        fprintf(fid,['...\n',repmat('\t',1,depth+1),')']);

   otherwise
       if isscalar(var)
         fprintf(fid,'%s',sprintf(' %g ',var));
       else
           fprintf(fid,'[%s]',sprintf(' %g ',var));
       end   
end
if newline
    fprintf(fid,';\n');
end