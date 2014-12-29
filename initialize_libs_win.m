
% startup % setup cvx

cur_file_path=mfilename('fullpath');
[path_str, file_name, ext] = fileparts(cur_file_path);
[parent_path_str, cur_dir_name, ext] = fileparts(path_str);
if ~(ismcc || isdeployed)
    % addpath(strcat(path_str,'/queuingmva'));
    %             addpath(strcat(path_str));
    addpath(strcat(path_str,'/opera'));
    addpath(strcat(path_str,'/xml_io_tools'));
    addpath(strcat(path_str,'/estimation_myself/code'));
    addpath(strcat(path_str,'/estimation_myself/data'));
    %             addpath(strcat(path_str,'/policy'));
    %             addpath(strcat(path_str,'/utils'));
    %             addpath(strcat(path_str,'/control'));
    %             addpath(strcat(path_str,'/figure'));
    %             addpath(strcat(path_str,'/recent'));
    %             addpath(strcat(path_str,'/utils/graphViz4Matlab'));
    %             addpath(strcat(path_str,'/utils/GraphViz2Mat1.2'));
    %
    %
    % %             obj.include_jar('model.jar','/home/zigorat/workspace2/model/dist/model.jar');
    % %             obj.include_jar('opera.jar','/home/zigorat/workspace2/Opera/dist/opera.jar');
    %            javaaddpath({strcat(path_str,'/opera/model.jar'),...
    %                strcat(path_str,'/opera/opera.jar')});
    disp(strcat(path_str,'\opera\operaAll-new.jar') );
    
    javaaddpath({strcat(path_str,'\opera\operaAll-new.jar') });
    
    %javac -target 1.5 -source 1.5  sed\*.java
    %jar cvf sed.jar sed/*.class
    %  javaaddpath({strcat(path_str,'\opera\classes') });
    javaaddpath({strcat(path_str,'\opera\sed.jar') });
    
    % %             javaaddpath({'/home/zigorat/workspace2/model/dist/model.jar',...
    % %                  '/home/zigorat/workspace2/Opera/dist/opera.jar'});
    %
    % cvx_path = 'C:\mf\tool\cvx'; %'/home/zigbie/my/tools';  %/cse/home/hamoun/


elseif true
    addpath(strcat(ctfroot, '/cvx'));
    addpath(strcat(ctfroot, '/cvx/structures'));
    addpath(strcat(ctfroot, '/cvx/lib'));
    addpath(strcat(ctfroot, '/cvx/functions')) ;
    addpath(strcat(ctfroot, '/cvx/commands'));
    addpath(strcat(ctfroot, '/cvx/builtins')) ;
    addpath(strcat(ctfroot, '/cvx/shims')) ;
    strcat(ctfroot, '/cvx/shims')
    try
        cvx_setup
        cvx_startup
        % run('C:\mf\tool\cvx\cvx_startup.m');
    catch me
        me
    end
else
    cvx_sdpt3
    cvx_sedumi
end
