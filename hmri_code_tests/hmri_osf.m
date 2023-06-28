classdef hmri_osf < handle
    % HMRI_OSF Reading access to osf storage nodes.
    %   This class provides an interface to interact with osf storage nodes
    %   for listing and downloading files and directories.
    %
    %   Creation:
    %       % Access the qMRI BIDS example data at https://osf.io/k4bs5/
    %       osf = hmri_osf()
    %       % Set output folder for downloaded content
    %       osf.target_root = '/tmp' % Set output folder
    %   
    %   Listing:
    %       % Get directory content as map of files and subdirectories.
    %       % Values of the map are structs with information about entries
    %       content = osf.ls()
    %       % Get a tree view of all files and directories
    %       osf.print_tree()
    %       % Get all names of entries inside the 'ds-mp2rage' directory
    %       content('ds-mp2rage').files.keys
    %
    %   Downloading:
    %       osf.download('<path to file/directory>')
    %       % Will download the given file or directory recursively.
    %       % Files that already exist, and are not older than the source
    %       % will not be downloaded again.
    %       osf.download('ds-mp2rage/sub-1/anat')
    %       % will create a directory 'ds-mp2rage/sub-1/anat' and fill that
    %       % with the content of the respective directory on osf.
    %   
    %   Tested with Matlab 2018a-2023a

    
    properties
        root_node = "k4bs5";       % Root node ID of the OSF project
        token = "";                % Access token for the OSF project
        target_root = string(pwd); % Target root directory for downloads
        tree_cache = [];
    end

    properties (Access=private)
        API_URL = "https://api.osf.io/v2"; % % Base URL for OSF API
        
    end


    methods (Access = private)
        function check_makedir(~,target)
            % Check if the target directory exists and create it if not
            if(isfolder(target)==0)
                mkdir(target)
            end
        end

        function req=gen_request(obj)
            % % Generate an HTTP request with appropriate headers
            if strlength(obj.token)
                hdr=matlab.net.http.HeaderField("Authorization", "Bearer " + string(obj.token));
            else
                hdr=[];
            end
            req=matlab.net.http.RequestMessage('GET',hdr);
        end

        function json = get_json(obj, url)
            % Get JSON data from a given URL
            try
                options = matlab.net.http.HTTPOptions;
                options.CertificateFilename = ''; % disable certificate check as workaround for older Matlab
                result = send(obj.gen_request(),url,options);
            catch ME
                error("Error accessing %s: %s",url,ME.message);
            end
            if result.StatusCode == matlab.net.http.StatusCode.OK
                json = jsondecode(char(result.Body));
            else
                error("Accessing %s resulted in %s",url,result.StatusCode);
            end
        end
        
        function files = node_ls(obj,node)
            % List all files and folders within a given node
            files = containers.Map;
            files_url = node.relationships.files.links.related.href;
            while(files_url)
                result = obj.get_json(files_url);
                data = result.data;
                for f = 1:length(data)
                    kind = string(data(f).attributes.kind);
                    if kind == "file"
                        entry.fullpath = string(data(f).attributes.materialized_path);
                        entry.date = datetime(data(f).attributes.date_modified,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSSSSS");
                        entry.size = data(f).attributes.size;
                        entry.download = data(f).links.download;
                        entry.isfile = true;
                        files(char(data(f).attributes.name)) = entry;
                    else
                        entry.fullpath = string(data(f).attributes.materialized_path);
                        entry.isfile = false;
                        entry.files = obj.node_ls(data(f));
                        files(char(data(f).attributes.name)) = entry;
                    end
                end
                files_url=result.links.next;
            end
        end
        
        function download_entry(obj,entry,target)
            % Download a single entry (file or directory) to the specified target
            import matlab.net.*
            import matlab.net.http.*
            if(entry.isfile)
                if(isfile(target{1}))
                    file = dir(target{1});
                    if(file.date > entry.date)
                        fprintf(2,"The existing file %s, is up-to-date. Won't download again\n", target);
                        return
                    end
                end
                fprintf("Downloading %s from %s\n",target,entry.download);
                options = matlab.net.http.HTTPOptions;
                options.CertificateFilename = ''; % disable certificate check as workaround for older Matlab
                send(obj.gen_request(),entry.download,options,io.FileConsumer(target{1}));
            else
                obj.check_makedir(target)
                for k = string(entry.files.keys)
                    obj.download_entry(entry.files(char(k)),target+"/"+k)
                end
            end
        end
        
        function download_r(obj, files, lh, rh)
            % Recursive download function
            target_dir = join([obj.target_root,lh],"/");
            obj.check_makedir(target_dir);
            if(numel(rh)==1) %leaf
                target = join([target_dir,rh],"/");
                obj.download_entry(files(char(rh(1))),target);
            else % branch
                new_lh = [lh,rh(1)];
                new_rh = rh(2:end);
                if(files.isKey(char(rh(1)))==1)
                    download_r(obj,files(char(rh(1))).files,new_lh,new_rh);
                else
                    error(join([lh,new_lh],"/")+" does not exist");
                end
            end
        end

        function indent = dash_indent(~, x)
            % Helper method
            pipe = char(9474);
            if x == 1
                indent = strcat(pipe, "   ");
            else
                indent = "    ";
            end
        end

        function tree_helper(obj, files, indents)
            % Visits files recursively and prints the file tree

            % Characters we need for beautiful output
            pipe_dash = char(9500);
            pipe_end = char(9492);
            dash = char(9472);

            % Each line starts with a certain indent the represents the
            % current depth of the entry in the filesystem. Also, we have
            % vertical pipes at specific places that indicate the parent
            % directory
            if isempty(indents)
                line_start = "";
            else
                line_start = join(arrayfun(@(x) obj.dash_indent(x), indents), "");
            end

            % Loop through each entry and print the file/directory with the
            % correct indent and line marks. For each directory, the
            % tree_helper function is called recursively.
            keys = files.keys();
            for i = 1:length(keys)
                entry = files(keys{i});
                if i == length(keys)
                    next_indent = 0;
                    entry_dash = strcat(pipe_end, dash);
                else
                    next_indent = 1;
                    entry_dash = strcat(pipe_dash, dash);
                end
                if entry.isfile == 1
                    fprintf("%s %s", strcat(line_start, entry_dash), keys{i});
                    fprintf(" (%s)\n", entry.download);
                else
                    fprintf("%s %s\n", strcat(line_start, entry_dash), keys{i});
                    tree_helper(obj, entry.files, [indents, next_indent]);
                end
            end
        end

    end


    methods
        function files = ls(obj, force_update)
            % List files and directories of the root node
            if ~exist('force_update','var')
                force_update = false;
            end

            % Download the tree only if it is not cached or we want to
            % enforce a redownload
            if force_update || isempty(obj.tree_cache)  
                node = obj.get_json(join([obj.API_URL,"nodes",obj.root_node,"files"],"/")).data;
                obj.tree_cache=obj.node_ls(node);
            end
            files=obj.tree_cache;
        end

        function print_tree(obj)
            % Prints out the directory structure like the Unix tree command
            tree_helper(obj, obj.ls(), [])
        end

        function download(obj,filename)
            % Download the specified file or directory
            obj.download_r(obj.ls(),[],split(string(filename),"/"));
        end
    end
end

