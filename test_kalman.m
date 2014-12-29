
import opera.*;
import opera.KalmanFilter.*; 
% http://www.mathworks.com/matlabcentral/answers/12422-macosx-encoding-problem
%  
cur_file_path=mfilename('fullpath');   [path_str, file_name, ext] = fileparts(cur_file_path);


theModel =  OperaModel(); % OperaModel
theModel.setModel(strcat(path_str,'/opera/output/Simple DB Operations.model.pxl'));

kalmanConfig = KalmanConfiguration(); %KalmanConfiguration
kalmanConfig.withConfigFile(strcat(path_str,'/opera/output/mykalman.config'))...
    .withSetting(KalmanConfiguration.ITERATIONS_MAX, '1')...
    .withSetting(KalmanConfiguration.FILE_MODEL_RESULTS, strcat(path_str,'/opera/output/Simple DB Operations.results.xml'))...
    .withModel(theModel)...
     .withSetting(KalmanConfiguration.MODEL_UPDATE, 'true')...
    .withSetting(KalmanConfiguration.FILE_TRACE, strcat(path_str,'/opera/output/Simple DB Operations.trace.txt'));

theEstimator =  KalmanEstimator(kalmanConfig); %KalmanEstimator
theEstimator.EstimateModelParameters().toString()

theModel.SaveModelToXmlFile(strcat(path_str,'/opera/output/testtestl.pxl'));
import javax.xml.transform.*;
a = TransformerFactory.newInstance()
 xmlTransformer = a.newTransformer();
 xmlTransformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, 'no');
	    	xmlTransformer.setOutputProperty(OutputKeys.INDENT, 'yes');
	    	xmlTransformer.setOutputProperty(OutputKeys.ENCODING, 'UTF-8');
	    	xmlTransformer.setOutputProperty('{http://xml.apache.org/xslt}indent-amount', '2');