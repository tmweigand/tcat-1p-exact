import os

def remove_file(file):
    """
    Remove files if exists
    """
    if os.path.isfile(file):
        os.remove(file)


def clean_directory():
    """
    Clean the directory
    """
    cwd = os.getcwd()
    remove_file(cwd + "/contraction_energies.txt")
    remove_file(cwd + "/packing_init.xyzd")
    remove_file(cwd + "/packing.nfo")
    remove_file(cwd + "/packing.xyzd")
    remove_file(cwd + "/packing_prev.xyzd")
