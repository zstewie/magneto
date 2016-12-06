% DROPBOX  Return root path of local Dropbox folder
%    PATH = DROPBOX() returns the path.
%    PATH = DROPBOX(RELATIVE) returns the full path to a file or folder
%    relative to the root Dropbox path.
%
% DROPBOX keeps track of the location of your Dropbox folder across
% different machines and operating systems. On first-time use, it asks the
% user to locate the local Dropbox folder through a GUI window. The m-file
% then updates itself to recognize different machines through their
% hostname, so that only one copy of the file is needed (that you can keep
% in your Dropbox). DROPBOX does not rely on MATLAB preferences so it is
% robust across MATLAB versions and after fresh installs.
%
% Because the file edits itself, user changes are not recommended; the name
% of the file can be changed, however.

% DO NOT EDIT THIS FILE %

function p = dropbox(varargin)

switch strtrim(evalc('system(''hostname'');'))
% Start hostname cases
    case 'Zacharys-MacBook-Pro-8.local'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a225363.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a237535.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22544b.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a2250fe.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a225218.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a238134.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'Zacharys-MacBook-Pro-9.local'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a228165.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22875f.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a238257.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a238419.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a2257b0.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22504b.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22581f.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a2384b3.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN52eigu.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22264e.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22540e.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a2265f1.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN52eir3.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22521f.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22647c.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a2259c4.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a23807b.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN52eide.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a238507.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN52eioj.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN525j9j.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a225956.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a226941.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN52ei13.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a237a67.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22671e.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a225063.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'rescomp-15-228877.stanford.edu'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN52eit6.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22534b.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22b9e7.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a225509.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a2250c5.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a225741.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN52eimo.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22564a.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a2255b9.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a221968.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22533a.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a2384b7.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a225177.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a225642.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a226213.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22567c.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22404a.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a22546a.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a225735.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
    case 'DN0a225632.SUNet'
        p = fullfile('/Users/Zach/StanfordDropbox/Dropbox',varargin{:});
%##%
% End hostname cases
    otherwise
        dbdr = uigetdir(pwd,'Please locate your Dropbox folder.');
        if ~dbdr
            error('Dropbox:dropbox:folderNotFound',...
                'Your Dropbox folder could not be located.')
        end
        tempfile = tempname();
        thisfile = which(mfilename);
        fi = fopen(thisfile,'r');
        fo = fopen(tempfile,'wt');
        if fo<0
            error('Dropbox:dropbox:tempnameFailed',...
                'Could not write temporary file to %s.',tempfile);
        end
        while ~feof(fi)
            buffer = fgetl(fi);
            if strcmp(buffer,'%##%')
                thishost = strtrim(evalc('system(''hostname'');'));
                fprintf(fo,'    case ''%s''\n',thishost);
                fprintf(fo,'        p = fullfile(''%s'',varargin{:});\n',dbdr);
            end
            fprintf(fo,'%s\n',buffer);
        end
        fclose(fi);
        fclose(fo);
        [success errmsg errid] = movefile(tempfile,thisfile,'f');
        if ~success
            error(errid,errmsg);
        end
        fprintf('Host ''%s'' added to %s.m!\n',thishost,mfilename);
end
