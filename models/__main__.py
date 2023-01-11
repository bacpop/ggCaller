import os
import tarfile
import wget

""" Get directories for model and seengenes """
module_dir = os.path.dirname(os.path.realpath(__file__))
zipped_db_dir = db_dir = os.path.join(module_dir, "ggCallerdb.tar.bz2")
db_dir = os.path.join(module_dir, "ggCallerdb")
balrog_model_dir = os.path.join(db_dir, "balrog_models")

def download_db():
    if not os.path.exists(zipped_db_dir):
        print("Downloading databases...")
        url = "https://ggcallerdb.blob.core.windows.net/dbs/ggCallerdb.tar.bz2"
        filename = wget.download(url, out=module_dir)
        print("")
    if not os.path.exists(db_dir):
        tar = tarfile.open(db_dir + ".tar.bz2", mode="r:bz2")
        tar.extractall(module_dir)
        tar.close()

    return db_dir

def load_balrog_models():
    # check if directory exists. If not, unzip file
    if not os.path.exists(balrog_model_dir):
        tar = tarfile.open(balrog_model_dir + ".tar.gz", mode="r:gz")
        tar.extractall(db_dir)
        tar.close()

    geneTCN = os.path.join(balrog_model_dir, "geneTCN_jit.pt")
    tisTCN = os.path.join(balrog_model_dir, "tisTCN_jit.pt")

    return (geneTCN, tisTCN)
