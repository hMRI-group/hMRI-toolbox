% % Unit tests class definition template
% %   
% 
% classdef hmri_function_template_test < matlab.unittest.TestCase
%     properties (TestParameter)
%         % Augment TestParameter with parameters over which tests will run,
%         % as well as parameters needed by the test functions.
%         specifyParamAsList = {1,10,100};
% 
%     end
%     properties
%         % Parameters independent of TestParameter
%     end
%     
%     methods (Test)
%         
%         %% Test Functions
%         function testSingleContrast(testCase,specifyParamAsList)
%             
%             % Check every element of output is equal to reference to within a 
%             % pre-defined threshold, defined as object property. 
%             % %TODO -rationale for threshold.
%             assertEqual(testCase,output,reference,'AbsTol',tolerance)
%         end
%         
%        
%  
%         function testNonStructInput(testCase)
%             
%             % Check that an error is thrown if the input isn't a struct
%             assertError(testCase, @() hmri_function(zeros([32,56,8])), 'hmri:structError');
%             
%         end
%         
%  
%     end
%     
%     %% Test Setup and Teardown Functions
%     
%     methods(TestMethodSetup)
%         % These methods are run before each test
%         
%         function methodName(testCase) %#ok<MANU>
%         end
%         
%     end
% 
%     methods(TestClassSetup)
%         % These methods are run when instantiating the class
%         
%         function methodName2(testCase) %#ok<MANU>
%         end
% 
%     end
%  
%     methods(TestMethodTeardown)
%         
%         function deleteTempData(testCase)
%             % Destructor 
%             delete(testCase.src_file)
%         end
%     end
%     
%     methods(TestClassTeardown)
% 
%     end
% end
