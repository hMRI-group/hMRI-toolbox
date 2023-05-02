classdef hmri_osf
    %Reading access to osf storage nodes.
    %   Creation:
    %   hmri_osf(<node_id>,<access_token>)
    %   eg: osf = hmri_osf("k4bs5","")
    %   for https://osf.io/k4bs5/
    %   
    %   listing: 
    %   osf.ls() will return associative list of directory entries (files or 
    %   sudirectories) and structs with information to those entries.
    %   
    %
    %   downloading:
    %   osf.download(<path to file/directory>)
    %   Will download the given file or directory recursively.
    %   Files that already exist, and are not older than the source will
    %   not be downloaded again.
    %   eg: osf.download('ds-mp2rage/sub-1/anat') will create a diretory 
    %   'ds-mp2rage/sub-1/anat' and fill that with the content of the
    %   respective directory on osf.
    %   
    %   Tested with Matlab 2019b-2023a 
    %   (2019a failes because of missing certs)

    
    properties
        API_URL = "https://api.osf.io/v2";
        target_root = string(pwd);
        root_node;
        token;
    end
   
    methods (Access = private)
        function check_makedir(obj,target)
            if(isfolder(target)==0)
                mkdir(target)
            end
        end

        function req=gen_request(obj)
            if strlength(obj.token)
                hdr=[matlab.net.http.HeaderField("Authorization", "Bearer " + string(obj.token))];
            else
                hdr=[];
            end
            req=matlab.net.http.RequestMessage('GET',hdr);
        end

        function json = get_json(obj, url)
            result = send(obj.gen_request(),url);
            json = jsondecode(string(result.Body));
        end
        
        function files = node_ls(obj,node)
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
                send(obj.gen_request(),entry.download,[],io.FileConsumer(target{1}));
            else
                obj.check_makedir(target)
                for k = string(entry.files.keys)
                    obj.download_entry(entry.files(k),target+"/"+k)
                end
            end
        end
        
        function download_r(obj, files, lh, rh)
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
            node = obj.get_json(join([obj.API_URL,"nodes",obj.root_node,"files"],"/")).data;
            files=obj.node_ls(node);
        end
        function obj = hmri_osf(root_node,token)
            obj.root_node = string(root_node);
            obj.token = token;
        end
        function download(obj,filename)
            obj.download_r(obj.ls(),[],split(string(filename),"/"));
        end

    end
end

