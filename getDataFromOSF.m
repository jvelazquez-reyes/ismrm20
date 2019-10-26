function getDataFromOSF

disp('Downloading dataset. DOI 10.17605/OSF.IO/A4V8H');

url = 'https://osf.io/4fztd/download/';

try 
    websave('data.zip',url);
catch 
    [~, SUCCESS, MESSAGE] = urlwrite(url,'data.zip');
    if ~SUCCESS, error(MESSAGE); end    
end
disp('Data has been downloaded ...');

unzip('data.zip'); 


end 