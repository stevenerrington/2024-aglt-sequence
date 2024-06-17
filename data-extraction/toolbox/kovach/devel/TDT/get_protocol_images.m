
  function out = get_protocol_images(blkdat)

    urlbase = 'http://cruprotocol.neurosurgery.uiowa.edu/';
    if isstruct(blkdat)
        prurl = blkdat.url;
    end

    html = urlread(prurl);

    scr = regexp(html,'ShowImage([^}]+{.*?var s =([^\n]*)','tokens','once');
    imgexp = regexprep(scr,'"\s*\+[^\+]*\+*\s*"*','%s');
    imgexp=regexprep(imgexp,'\s*"','');
    imgurl = @(imgid)sprintf(['%s',imgexp{1}],urlbase,'final',num2str(blkdat.record_id),blkdat.subject_id,imgid);

    ddc = regexp(html,'<select name="ddlFinalGrid"(.*?)</select','match');
    ddc = regexp(ddc{1},'<option value="(.*?)">(.*?)</option>','tokens');

    for i = 1:length(ddc)

        url2 = imgurl(ddc{i}{1});
        html2 = urlread(url2);

        imloc = regexp(html2,'<img.*?src="(.*?)"','tokens','once');
        imlocurl = [urlbase,imloc{1}];
        
         out.gridmaps(i) = protocol_images(imlocurl);
        out.gridmaps(i).label = ddc{i}{2};
        out.gridmaps(i).url = imlocurl; 
    end
    out.blkdat = blkdat;
     

    