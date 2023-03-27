classdef hmri_webdav
    %WebDav access to public ownCloud instances.
    %   Creation:
    %   dv = hmri_webdav(<server_url>,<ownCloud token>)
    %   eg: dv=hmri_webdav('https://owncloud.gwdg.de','iv2TOQwGy4FGDDZ')
    %   for https://owncloud.gwdg.de/s/iv2TOQwGy4FGDDZ
    %   
    %   listing: 
    %   dv.ls() for the root or (obj.ls('path')
    %   that will return associative list of directory entries (files or 
    %   sudirectories) and structs with information to those entries.
    %   
    %   string(dv.ls().keys) gives a simple list of the filenames.
    %   dv.ls();ans('README.txt') gives metadata for the file README.txt
    %   
    %
    %   downloading:
    %   obj.download(<path to file/directory>)
    %   Will download the given file or directory recursivly.
    %   Files that already exist, and are not older than the source will
    %   not be downloaded again.
    
    properties
        API_PATH = "remote.php/dav/public-files"
        target_root = string(pwd);
        token;
        SERVER_URI;
    end
    properties (Access = private)
        content_xml = matlab.net.http.field.ContentTypeField('application/xml; charset=UTF-8');
        depth_hdr = matlab.net.http.HeaderField('Depth', '1');
    end
    
    methods
        function obj = hmri_webdav(SERVER_URI,token)
            obj.SERVER_URI = string(SERVER_URI);
            obj.token = string(token);
        end

        function references = ls(obj,varargin)
            path = join([obj.token,varargin],'/');
            URL = join([obj.SERVER_URI,obj.API_PATH,path],'/');
            references = obj.ls_url(URL);
        end
        
        function references = ls_url(obj,URL)
            import matlab.net.*
            import matlab.net.http.*
            import matlab.io.xml.dom.*
            
            REQUEST_BODY = '<?xml version="1.0"?><d:propfind xmlns:d="DAV:" xmlns:oc="http://owncloud.org/ns"></d:propfind>';
            
            r = RequestMessage('PROPFIND',[obj.content_xml,obj.depth_hdr],REQUEST_BODY);
            result = send(r,URL);
            xml = parseString(Parser,string(result.Body));
            elements = xml.getElementsByTagName('d:response');
            references = containers.Map;
            for n=1:elements.Length-1 %first entry is directory itself (and yes it starts with "0". Yay consistency!!)
                element = elements.item(n);
                
                last_modified = element.getElementsByTagName('d:getlastmodified');
                contentlength = element.getElementsByTagName('d:getcontentlength');
                contenttype = element.getElementsByTagName('d:getcontenttype');

                if(last_modified.Length)
                    datestr = last_modified.TextContent(1:end-4);
                    timezone = last_modified.TextContent(end-2:end);
                    entry.date = datetime(datestr,'TimeZone',timezone,'InputFormat','eee, d MMM yyyy H:mm:ss');
                end
                if(contentlength.Length)
                    entry.size = str2double(contentlength.TextContent);
                end
                if(contenttype.Length)
                    entry.ctype = contenttype.TextContent;
                end

                entry.isdir = element.getElementsByTagName('d:collection').Length;

                ref = element.getElementsByTagName('d:href').TextContent;
                entry.url = strcat(obj.SERVER_URI,strip(ref,'right','/'));

                slash = strfind(ref,'/');
                filename = extractAfter(ref,slash(end));
                references(filename) = entry;
            end
        end
        
        function download_entry(obj,entry,target)
            import matlab.net.*
            import matlab.net.http.*
            if(isfile(target{1}))
                filedate = dir(target{1}).date;
                if(filedate > entry.date)
                    fprintf(2,"The existing file %s, is us-to-date. Won't download again\n", target);
                    return
                end
            end
            if(entry.isdir)
                if(isfolder(target{1})==0)
                    mkdir(target{1})
                end

                files = obj.ls_url(entry.url);
                for file = files.keys
                    download_entry(obj,files(file{1}),join([target,file],'/'));
                end
            else
                fprintf("Downloading %s from %s\n",target,obj.SERVER_URI);
                send(RequestMessage('GET'),entry.url,[],io.FileConsumer(target{1}));
            end
        end
        
        function download_r(obj, lh, rh)
            target_dir = join([obj.target_root,lh],'/');
            if(isfolder(target_dir)==0)
                mkdir(target_dir{1})
            end
            if(numel(rh)==1)
                target = join([target_dir,rh],'/');
                if(numel(lh)>0)
                    filemap = obj.ls(join(lh,'/'));
                else
                    filemap = obj.ls();
                end
                obj.download_entry(filemap(rh{1}),target);
            else
                new_lh = [lh,rh(1)];
                new_rh = rh(2:end);
                newmap = obj.ls(join(new_lh,'/'));
                if(newmap.isKey(rh{2})==1)
                    download_r(obj,new_lh,new_rh);
                else
                    error(join([rh{2},"is","not","in",new_lh]));
                end
            end
        end
        
        function download(obj,filename)
            obj.download_r([],split(filename,'/'));
        end

    end
end

