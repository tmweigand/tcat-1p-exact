import numpy as np


def volumeIntegrate(file_in, regions, variable):
    """
    Generate the openfoam function to integrate volume
    """
    out_file = open(file_in, "w", encoding="utf-8")
    for r in regions:
        name = str(r)
        out_file.write(f"{name}\n")
        out_file.write("{ \n")
        out_file.write("\t type \t volFieldValue; \n")
        out_file.write("\t libs \t (fieldFunctionObjects); \n")
        out_file.write("\t writeControl \t timeStep; \n")
        out_file.write("\t writeInterval \t 1; \n")
        out_file.write("\t writeFields \t false; \n")
        out_file.write("\t regionType \t cellZone; \n")
        out_file.write(f"\t name \t {r};\n")
        out_file.write("\t operation \t volIntegrate; \n")
        out_file.write("\t fields \n")
        out_file.write("\t( \n")
        for v in variable:
            out_file.write(f"\t \t {v} \n")
        out_file.write("\t ); \n")
        out_file.write("} \n \n")


def surfaceIntegrate(file_in, regions, variable, interface = False):
    """
    Generate the openfoam function to integrate volume
    """
    out_file = open(file_in, "w", encoding="utf-8")
    for r in regions:
        name = str(r)
        if interface:
            out_file.write(f"{name+'_ws'}\n")
        else:
            out_file.write(f"{name}\n")
        out_file.write("{ \n")
        out_file.write("\t type \t surfaceFieldValue; \n")
        out_file.write("\t libs \t (fieldFunctionObjects); \n")
        out_file.write("\t writeControl \t timeStep; \n")
        out_file.write("\t writeInterval \t 1; \n")
        out_file.write("\t writeFields \t false; \n")
        out_file.write("\t regionType \t faceZone; \n")
        if interface:
            out_file.write(f"\t name \t {name+'_ws'};\n")
        else:
            out_file.write(f"\t name \t {name};\n")
        out_file.write("\t operation \t areaIntegrate; \n")
        out_file.write("\t fields \n")
        out_file.write("\t( \n")
        for v in variable:
            out_file.write(f"\t \t {v} \n")
        out_file.write("\t ); \n")
        out_file.write("} \n \n")


def boundarySurfaceIntegrate(file_in, regions, patches, variable):
    """
    Generate the openfoam function to integrate volume
    """
    out_file = open(file_in, "w", encoding="utf-8")
    for region, patch in zip(regions,patches):
        name_region = str(region)
        name_patch = str(patch)

        out_file.write(f"{name_region}\n")
        out_file.write("{ \n")
        out_file.write("\t type \t surfaceFieldValue; \n")
        out_file.write("\t libs \t (fieldFunctionObjects); \n")
        out_file.write("\t writeControl \t timeStep; \n")
        out_file.write("\t writeInterval \t 1; \n")
        out_file.write("\t writeFields \t false; \n")
        out_file.write("\t regionType \t patch; \n")
        out_file.write(f"\t name \t {name_patch};\n")
        out_file.write("\t operation \t areaIntegrate; \n")
        out_file.write("\t fields \n")
        out_file.write("\t( \n")
        for v in variable:
            out_file.write(f"\t \t {v} \n")
        out_file.write("\t ); \n")
        out_file.write("} \n \n")