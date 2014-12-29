            
            % startup % setup cvx

            cur_file_path=mfilename('fullpath');
            [path_str, file_name, ext] = fileparts(cur_file_path);
            [parent_path_str, cur_dir_name, ext] = fileparts(path_str);            
         
             % addpath(strcat(path_str,'/queuingmva'));
             
%             addpath(strcat(path_str));
             addpath(strcat(path_str,'/opera'));
             addpath(strcat(path_str,'/xml_io_tools'));
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
             javaaddpath({strcat(path_str,'/opera/operaAll.jar') }); 
             
% %             javaaddpath({'/home/zigorat/workspace2/model/dist/model.jar',...
% %                  '/home/zigorat/workspace2/Opera/dist/opera.jar'});
% 
% cvx_path = '/Users/hamoun/mf/lib'; %'/home/zigbie/my/tools';  %/cse/home/hamoun/
% addpath(strcat(cvx_path, '/cvx'));
% addpath(strcat(cvx_path, '/cvx/structures'));
%    addpath(strcat(cvx_path, '/cvx/lib'));
%     addpath(strcat(cvx_path, '/cvx/functions')) ;
%     addpath(strcat(cvx_path, '/cvx/commands')); 
%     addpath(strcat(cvx_path, '/cvx/builtins')) ;
%     
    