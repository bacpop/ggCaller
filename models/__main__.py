import os
import tarfile
import wget

""" Get directories for model and seengenes """
module_dir = os.path.dirname(os.path.realpath(__file__))
module_zipped_db_dir = os.path.join(module_dir, "ggCallerdb.tar.bz2")
module_db_dir = os.path.join(module_dir, "ggCallerdb")
module_balrog_model_dir = os.path.join(module_db_dir, "balrog_models")

def download_db(download_db=None):
    if download_db is None:
        zipped_db_path = module_zipped_db_dir
        db_path = module_db_dir
        output_dir = module_dir
    else:
        zipped_db_path = os.path.join(download_db, "ggCallerdb.tar.bz2")
        db_path = os.path.join(download_db, "ggCallerdb")
        output_dir = download_db
    
    if not os.path.exists(zipped_db_path):
        print("Downloading databases...")
        url = "https://ftp.ebi.ac.uk/pub/databases/pp_dbs/ggCallerdb.tar.bz2"
        filename = wget.download(url, out=output_dir)
        print("")
    if not os.path.exists(db_path):
        tar = tarfile.open(zipped_db_path, mode="r:bz2")
        tar.extractall(output_dir)
        tar.close()

    return db_path

def load_balrog_models(db_path):
    balrog_model_dir = os.path.join(db_path, "balrog_models")
    
    # check if directory exists. If not, unzip file
    if not os.path.exists(balrog_model_dir):
        tar = tarfile.open(balrog_model_dir + ".tar.gz", mode="r:gz")
        tar.extractall(db_path)
        tar.close()

    geneTCN = os.path.join(balrog_model_dir, "geneTCN_jit.pt")
    tisTCN = os.path.join(balrog_model_dir, "tisTCN_jit.pt")

    return (geneTCN, tisTCN)
