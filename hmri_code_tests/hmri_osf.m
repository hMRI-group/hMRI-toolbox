classdef hmri_osf
    % HMRI_OSF Reading access to osf storage nodes.
    %   This class provides an interface to interact with osf storage nodes
    %   for listing and downloading files and directories.
    %
    %   Creation:
    %       % Access the hMRI Toolbox example data at https://osf.io/k4bs5/
    %       osf = hmri_osf()
    %       % Set output folder for downloaded content
    %       osf.target_root = "/tmp" % Set output folder
    %   
    %   Listing:
    %       % Get directory content as map of files and subdirectories.
    %       % Values of the map are structs with information about entries
    %       content = osf.ls()
    %       % Get all names of entries inside the "ds-mp2rage" directory
    %       content("ds-mp2rage").files.keys
    %
    %   Downloading:
    %       osf.download("<path to file/directory>")
    %       % Will download the given file or directory recursively.
    %       % Files that already exist, and are not older than the source
    %       % will not be downloaded again.
    %       osf.download('ds-mp2rage/sub-1/anat')
    %       % will create a diretory 'ds-mp2rage/sub-1/anat' and fill that
    %       % with the content of the respective directory on osf.
    %   
    %   Tested with Matlab 2019b-2023a 
    %   (2019a failes because of missing certs)

    
    properties
        root_node = "k4bs5";       % Root node ID of the OSF project
        token = "";                % Access token for the OSF project
        target_root = string(pwd); % Target root directory for downloads
    end

    properties (Access=private)
        API_URL = "https://api.osf.io/v2"; % % Base URL for OSF API
    end


    methods (Access = private)
        function check_makedir(obj,target)
            % Check if the target directory exists and create it if not
            if(isfolder(target)==0)
                mkdir(target)
            end
        end

        function req=gen_request(obj)
            % % Generate an HTTP request with appropriate headers
            if strlength(obj.token)
                hdr=[matlab.net.http.HeaderField("Authorization", "Bearer " + string(obj.token))];
            else
                hdr=[];
            end
            req=matlab.net.http.RequestMessage('GET',hdr);
        end

        function json = get_json(obj, url)
            % Get JSON data from a given URL
            result = send(obj.gen_request(),url);
            json = jsondecode(string(result.Body));
        end
        
        function files = node_ls(obj,node)
            % List all files and folders within a given node
            files = containers.Map;
            files_url = node.relationships.files.links.related.href;
            data = obj.get_json(files_url).data;
            for f = 1:length(data)
                kind = string(data(f).attributes.kind);
                if kind == "file"
                    entry.fullpath = string(data(f).attributes.materialized_path);
                    entry.date = datetime(data(f).attributes.date_modified,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss.SSSSSS");
                    entry.size = data(f).attributes.size;
                    entry.download = data(f).links.download;
                    entry.isfile = true;
                    files(string(data(f).attributes.name)) = entry;
                else
                    entry.fullpath = string(data(f).attributes.materialized_path);
                    entry.isfile = false;
                    entry.files = obj.node_ls(data(f));
                    files(string(data(f).attributes.name)) = entry;
                end
            end
        end
        
        function download_entry(obj,entry,target)
            % Download a single entry (file or directory) to the specified target
            fprintf("Ping 1")
            import matlab.net.*
            import matlab.net.http.*
            fprintf("Ping 2")
            if(entry.isfile)
                if(isfile(target{1}))
                    file = dir(target{1});
                    if(file.date > entry.date)
                        fprintf(2,"The existing file %s, is up-to-date. Won't download again\n", target);
                        return
                    end
                end
                fprintf("Downloading %s from %s\n",target,entry.download);
                send(obj.gen_request(),entry.download,[],io.FileConsumer(target{1}));
            else
                obj.check_makedir(target)
                for k = string(entry.files.keys)
                    obj.download_entry(entry.files(k),target+"/"+k)
                end
            end
        end
        
        function download_r(obj, files, lh, rh)
            % Recursive download function
            target_dir = join([obj.target_root,lh],'/');
            obj.check_makedir(target_dir);
            if(numel(rh)==1) %leaf
                target = join([target_dir,rh],'/');
                obj.download_entry(files(rh(1)),target);
            else % branch
                new_lh = [lh,rh(1)];
                new_rh = rh(2:end);
                if(files.isKey(rh(1))==1)
                    download_r(obj,files(rh(1)).files,new_lh,new_rh);
                else
                    error(join([lh,new_lh],"/")+" does not exist");
                end
            end
        end
        
    end

    methods
        function files = ls(obj)
            % List files and directories of the root node
            node = obj.get_json(join([obj.API_URL,"nodes",obj.root_node,"files"],"/")).data;
            files=obj.node_ls(node);
        end

        function download(obj,filename)
            % Download the specified file or directory
            obj.download_r(obj.ls(),[],split(string(filename),"/"));
        end
    end
end

