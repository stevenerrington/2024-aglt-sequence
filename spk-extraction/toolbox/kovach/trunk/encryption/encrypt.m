function [ciphertext,header] = encrypt(plaintext,filename,password)

% Convert plaintext to password-encrypted ciphertext and, optionally write
% to file. This uses encryption tools provided by java. Plaintext must be a
% character array, while the ciphertext and header output are signed 8-bit
% integers.

%C. Kovach 2016


import javax.crypto.*
import javax.crypto.spec.*
import java.security.*       


out1 = 'x';
out2 = 'y';
while ~isequal(out1,out2) && (nargin < 3 || isempty(password))
    try
        fprintf('\nEnter password for encryption: ')
        [~,out1] = system('read'); 
        fprintf('\nEnter password again: ')
        [~,out2] = system('read');
    catch
           out1 = getpw(' ') ;
           out2 = getpw('Enter password again:') ;
    end 
    if ~strcmp(out1,out2)
        fprintf('\n\nPasswords don''t match!\n')
    end
   
    password = out1;
end


digest = 'SHA'; % Algorithm to generate checksum

rng = SecureRandom(); 
mdg = MessageDigest.getInstance(digest);
salt = int8([rng.nextInt(2^8) rng.nextInt(2^8)]-128)';
check = mdg.digest(uint8(plaintext));

factory = SecretKeyFactory.getInstance('PBKDF2WithHmacSHA1');
spec = PBEKeySpec(password, salt, 65536, 128);
tmp = factory.generateSecret(spec);
secret = SecretKeySpec(tmp.getEncoded(), 'AES');

cipher = Cipher.getInstance('AES/CBC/PKCS5Padding');
cipher.init(Cipher.ENCRYPT_MODE, secret);
iv = cipher.getIV;
ciphertext = cipher.doFinal(int8(double(plaintext)-128));
header = [int8(length(digest));int8(double(digest)-128)';check;int8(length(salt));salt;int8(length(iv));iv];    

if nargin >1 && ~isempty(filename)

    fid = fopen(filename,'w');
    fwrite(fid,[header;ciphertext],'int8');
end

