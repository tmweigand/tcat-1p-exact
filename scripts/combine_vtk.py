import pyvista as pv
import glob
import os
import argparse

def process_vtk(directory_path):
    """
    Merge the vtk files in the given directory.

    Args:
        directory_path (str): The path to the directory to process.
    """

    # Check if the provided path is a valid directory
    if not os.path.isdir(directory_path):
        print(f"The path {directory_path} is not a valid directory.")
        return
    
    # List of VTK files to combine
    vtk_files = glob.glob(directory_path+'/tcat/local/*vtk')

    # Create an empty mesh
    combined_mesh = pv.UnstructuredGrid()

    # Loop through the files and append each mesh to the combined mesh
    for vtk_file in vtk_files:
        # Read the VTK file
        mesh = pv.read(vtk_file)
        # Append to the combined mesh
        combined_mesh = combined_mesh.merge(mesh)

    # Save the combined mesh to a new VTK file
    combined_mesh.save(directory_path+'/tcat/combined_local.vtk')

def main():
    # Parse the command-line arguments
    parser = argparse.ArgumentParser(description="Process a directory.")
    parser.add_argument("directory", help="The path to the directory to process")
    
    # Get the directory path from the command-line input
    args = parser.parse_args()
    
    # Call the function with the provided directory path
    process_vtk(args.directory)

if __name__ == "__main__":
    main()