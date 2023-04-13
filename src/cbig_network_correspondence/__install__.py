from git import Repo
from os import path

def main():
    print("Download data")
    Repo.clone_from('https://github.com/rubykong/cbig_network_correspondence_data', path.dirname(path.abspath(__file__)))
    
if __name__ == '__main__':
    main()