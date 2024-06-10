classdef protocol_images < handle 
    
    properties
        label = '';
        url = '';
    end
    properties (Dependent = true)
        image = [];
    end
    properties (Hidden = true)
        imgdat = [];
    end
    methods
      function me=protocol_images(url)
            me.url = url;
      end
      
      function img = get.image(me)
          
          if isempty(me.imgdat)
             me.imgdat = me.retrieve(); 
          end
          img = me.imgdat;
      end
      
      function X = retrieve(me)

        imdat = readurl_binary(me.url);
        [~,ext] = strtok(imdat.contentType,'/');
        tempfn = sprintf('%s.%s',tempname,ext(2:end));
        fid = fopen(tempfn,'w');
        fwrite(fid,imdat.data);
        fclose(fid);
        X = imread(tempfn);
        delete(tempfn);
      end
      
    end
end
    