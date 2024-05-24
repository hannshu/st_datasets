import os
from tqdm import tqdm
import sys
import json
from urllib.request import Request, ProxyHandler, build_opener
 

def check_file_location(path, url):
    # data not in local, download from url
    if (False == os.path.exists(path)):
        proxy = check_proxy()
        # download and save file to local path
        try:
            if (False == proxy):
                opener = build_opener()
            else:
                # set proxy server
                handler = ProxyHandler(proxy)
                opener = build_opener(handler)

            # Request to download the data
            request = Request(url)
            response = opener.open(request)

            # build file folder
            os.makedirs(os.path.dirname(path), exist_ok=True)
            file_size = int(response.headers['Content-Length'])

            # download file
            data = []
            with tqdm(total=file_size, unit='B', unit_scale=True, unit_divisor=1024, desc='>>> INFO: Download dataset') as bar:
                for _ in range(file_size):
                    chunk = response.read(1024)
                    if (not chunk):
                        break
                    data.append(chunk)
                    bar.update(len(chunk))

            with open(path, "wb") as file:
                file.write(b''.join(data))

            # Close the response
            response.close()
        except Exception as e:
            # If an exception occurs, delete the partially downloaded file
            print(f">>> ERROR: Download interrupted: {e}")
            sys.exit()
    else:
        print('>>> INFO: Use local data.')

    return path


def check_proxy():
    proxy_file_path = os.path.join(os.environ['HOME'], 'st_datasets_proxy_setting.json')
    if (False == os.path.exists(proxy_file_path)):
        generate_proxy_file(proxy_file_path)
        raise FileNotFoundError(f'This is the first time you use st_datasets, please confirm proxy setting at {proxy_file_path}! If you do not need to set the proxy service, please directly rerun your script.')
    else:
        try:
            with open(proxy_file_path, 'r') as file:
                content = json.load(file)
            file.close()
            if (False == content['proxy']):
                proxy = False
            else:
                proxy = {
                    "http": content['http_proxy'], 
                    "https": content['http_proxy']
                }
        except Exception:
            print(f">>> ERROR: proxy setting file error, we will delete this file and generate a new one!")
            generate_proxy_file(proxy_file_path)
            sys.exit()
    
    return proxy


def generate_proxy_file(json_path):
    data = {
        "proxy": False,
        "http_proxy": "http://hostname:port",
        "https_proxy": "https://hostname:port",
        "_instruction": "The data contained within the st_datasets is stored in remote database. If your native network cannot access Hugging Face services and the StomicsDB ftp database, you should utilize proxy services to download those data or download the datasets from Baidu Netdisk we provided in README file and place it in your user directory. If you require the use of proxy services, please change 'proxy' to `true`, and modify the 'proxy,' 'http_proxy,' and 'https_proxy' according to the format specified above."
    }
    with open(json_path, 'w') as file:
        json.dump(data, file, indent=4)
