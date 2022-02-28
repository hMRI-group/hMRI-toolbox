classdef hmri_test_utils
    
    % Class with general utility methods for unit testing the hMRI toolbox.
    
    methods(Static)
        % Ernst equation
        function S=ernst(alpha,TR,R1)
            
            S=sin(alpha).*(1-exp(-TR.*R1))./(1-cos(alpha).*exp(-TR.*R1));
            
        end
        
        % TODO: warn about activating legacy RNG
        function seedRandomNumberGenerator
            % set randomness method and seed for reproducibility
            %
            % it would be best to use a separate random stream for each test
            % for simplicity we set the global stream
            % this will seed the RNG separately on every worker
            rng(319, 'twister');
        end
        
    end
end