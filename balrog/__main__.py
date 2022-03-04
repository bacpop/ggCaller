import os
import tarfile

""" Get directories for model and seengenes """
module_dir = os.path.dirname(os.path.realpath(__file__))
model_dir = os.path.join(module_dir, "balrog_models")

def load_gene_models():
    # check if directory exists. If not, unzip file
    if not os.path.exists(model_dir):
        tar = tarfile.open(model_dir + ".tar.gz", mode="r:gz")
        tar.extractall(module_dir)
        tar.close()

    geneTCN = os.path.join(model_dir, "geneTCN_jit.pt")
    tisTCN = os.path.join(model_dir, "tisTCN_jit.pt")

    return (geneTCN, tisTCN)
