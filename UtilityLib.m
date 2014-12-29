classdef UtilityLib
   properties

        color = { 'r' 'g' 'b' 'c'	'm'	'y'	'k'	'w'	};
        lnstyle = {'-','--','-.',':'};
        marker = {'s','o','h','.','x','+','d','^','v','>','<','*'};
   end
    methods(Static)
        function print_figure(figure,width,height,filename)
             set(gcf, 'PaperUnits', 'inches');
             set(gcf, 'PaperSize', [8.5 11]);
             % papersize = get(gcf, 'PaperSize')
%              width = 5.5;         % Initialize a variable for width.
%              height = 3;          % Initialize a variable for height.
             % left = (papersize(1)- width)/2
             % bottom = (papersize(2)- height)/2
             myfiguresize = [0, 0, width, height];
             set(gcf, 'PaperPosition', myfiguresize);
             set(gcf, 'PaperPositionMode', 'auto');
             print('-depsc',figure,filename); %add depsc to be colored %,'-tiff' %-deps -djpeg -depsc
        end
         
        function test
           handle = figure;
           plot(sin(1:100));
           UtilityLib.print_figure(handle,9,7,'figures\workload-interarr');
        end
        
    end
end






