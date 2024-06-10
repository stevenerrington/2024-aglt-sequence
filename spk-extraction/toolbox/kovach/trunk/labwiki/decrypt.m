function [plaintext,ciphertext,header,filename] = decrypt(ciphertext,header,password)

%    plaintext = decrypt(ciphertext,header)
%
% Decrypts ciphertext encrypted by ENCRYPT.
%
%    plaintext = decrypt(filename)
%
% Reads data from file.
%
%    plaintext = decrypt
%
% Prompts for file with encrypted data.
%
%
% See also ENCRYPT

%C. Kovach 2016



import javax.crypto.*
import javax.crypto.spec.*
import java.security.*       


if nargin < 2 || isempty(header) % assume file if nargin < 2

    if nargin < 1 || isempty(ciphertext)
        fprintf('\nLocate encrypted key file.\n')
        [fname,pth] = uigetfile({'*.enc'},'Locate file with encrypted data.');
        filename = fullfile(pth,fname);
    else
        
        filename = ciphertext;
    end
    
    fid = fopen(filename,'r');
    
    %%% Read the header
    dglen = fread(fid,1,'int8=>int8');
    
    dg = char(fread(fid,dglen,'int8=>double')+128)';
    
    mdg = MessageDigest.getInstance(dg);
    check = fread(fid,mdg.getDigestLength,'int8=>int8');
    
    sltlen = fread(fid,1,'int8=>int8');
    salt = fread(fid,sltlen,'int8=>int8');
    
    ivlen  = fread(fid,1,'int8=>int8');
    iv = fread(fid,ivlen,'int8=>int8');
    
    if nargout > 2
        header = [int8(length(dg));int8(double(dg)-128)';check;int8(length(salt));salt;int8(length(iv));iv];    
    end
    %%% Done reading header; the rest is ciphertext. 
    ciphertext = fread(fid,'int8=>int8');

    fclose(fid);
else
     dglen = double(header(1));
     indx = 1;
     dg = char(double(header((1:dglen)+indx))+128);
      mdg = MessageDigest.getInstance(dg);
     indx = indx+dglen;
     check = header((1:mdg.getDigestLength)+indx);
     indx = indx+mdg.getDigestLength;
     sltlen = double(header(indx+1));
     indx = indx+1;
     salt = header((1:sltlen)+indx);
     indx = indx+sltlen;
     ivlen = double(header(indx+1));
     indx = indx+1;
     iv = header((1:ivlen)+indx);
end
 
if nargin < 3 || isempty(password)
    fprintf('\nEnter password to decrypt input: ')
    try
        [~,password] = system('read'); 
    catch
        password = getpw(' ');
    end
end

factory = SecretKeyFactory.getInstance(java.lang.String('PBKDF2WithHmacSHA1'));
spec = PBEKeySpec(password, salt, 65536, 128);
tmp = factory.generateSecret(spec);
secret = SecretKeySpec(tmp.getEncoded(), 'AES');

cipher = Cipher.getInstance('AES/CBC/PKCS5Padding');

cipher.init(Cipher.DECRYPT_MODE, secret, IvParameterSpec(iv));
try
    plaintext = char(double(cipher.doFinal(ciphertext))+128)';

catch
    error(sprintf('\nDecryption returned an error\nprobably because the password is incorrect.\n'))
end
 
if ~isequal(check,mdg.digest(uint8(plaintext)))
    fprintf('Checksum failed. Wrong password maybe?')
    plaintext = '';
    return
end


