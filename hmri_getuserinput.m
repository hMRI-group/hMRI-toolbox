    function output = getuserinput(mainstr)

prompt = {'do you want to input customized parameters: type 0 OR 1' };
dlgtitle = 'User customized variable?';
dims = [1 35];
definput = {'false'};
answer1 = inputdlg(prompt,dlgtitle,dims,definput);


    if str2double(answer1{1})
prompt = {'Enter number of parameters you want to input (nonzero entry)' };
dlgtitle = 'Number of user customized parameters';
dims = [1 35];
definput = {'0'};
answer2 = inputdlg(prompt,dlgtitle,dims,definput);
        if str2double(answer2{1})~=0
            for i=1:str2double(answer2{1})
prompt = {'Parameter name:', 'Parameter value:', 'Parameter type: (char or num)' };
dlgtitle = 'User customized variable';
dims = [1 35];
definput = {'name','value', 'type'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
disp(answer)
                if strcmp(answer{3},'char')
          
mainstr= setfield(mainstr, answer{1}, answer{2});
                end
                if strcmp(answer{3}, 'num')
                   
    mainstr= setfield(mainstr, answer{1}, str2num(answer{2}));
                end
            end
        else

        end
    end
output = mainstr;

    end