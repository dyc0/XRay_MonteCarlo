from material_generator import read_attenuation
from spectra_converter import read_spectrum
from efficiency_generator import load_csv_data

from matplotlib import pyplot as plt
import numpy as np
import os
import string
import sys

def get_map(keys: list, values: list, name:str, type1:str, type2:str, comment:str) -> string:
    dict = ""

    dict += "\t// " + comment + "\n"
    dict += "\tstatic const std::map<" + type1 + ", " + type2 + "> " + name + " = {\n"
    elements = ''
    for e, p in zip(keys, values):
        if type1 not in ("std::string", "string"):
            elements += "\t\t{" + str(e) + ", " + str(p) + "},\n"
        else:
            elements += "\t\t{\"" + str(e) + "\", " + str(p) + "},\n"
    elements = elements[:-2]
    elements += "\n"
    dict += elements
    dict += "\t};\n"

    return dict

def get_list(elements: list, name:str, type:str, comment:str) -> string:
    container = ""

    container += "\t// " + comment + "\n"
    container += "\tstatic const std::vector<" + type + "> " + name + " = {\n"
    elements_str = ''
    for e in elements:
        elements_str += "\t\t" + str(e) + ",\n"
    elements_str = elements_str[:-2]
    elements_str += "\n"
    container += elements_str
    container += "\t};\n"

    return container

if __name__ == '__main__':
    argument_list = sys.argv[1:]
    options = "k:d:f:s:m:t:o:"
    long_options = ["kvp=", "deg=", "filt=", "SpectrumPath=", "MaterialsPath=", "Tolerance=", "OutputFile="]

    tube_kvp = "60"
    tube_deg = "17"
    tube_flt = "1Al"
    spectrum_path = "spectra"
    materials_path = "materials_data/attenuations"
    efficiencies_path = "efficiency_data/compound_data"
    tolerance = "1e-6"
    output_file = "../LIBRARY/include/constants.hpp"

    energies = np.arange(6, 120.5, 0.5)

    materials = {}
    energies_MeV = [e*0.001 for e in energies]
    for filename in os.listdir(materials_path):
        mat_energies, mat_attenuation = read_attenuation(os.path.join(materials_path, filename))
        mat_attenuation_interp = np.interp(energies_MeV, mat_energies, mat_attenuation)
        materials[filename[:-4]] = mat_attenuation_interp

    det_effs = {}
    for filename in os.listdir(efficiencies_path):
        det_eff = load_csv_data(os.path.join(efficiencies_path, filename))
        det_eff_interp = np.interp(energies, det_eff[0], det_eff[1])
        det_effs[filename[:-4]] = det_eff_interp

    

    with open(output_file, 'w') as of:
        of.write("#ifndef CONSTANTS_HPP\n#define CONSTANTS_HPP\n\nnamespace xrc {\n\n")

        tolerance_string = "\tconst static double tolerance = " + tolerance + ";\n\n"
        of.write(tolerance_string)

        pi_string = "\tconst static double PI = 3.14159265358979323846;\n\n"
        of.write(pi_string)

        source_e = get_list(energies, "energies", "double", "Possible energies of photons")
        of.write(source_e)
        of.write("\n\n")

        spectra_enum_names = []
        spectra_pointers = []
        for kvp in ['60', '120']:
            for filt in ['1Al', '3Al']:
                spectrum_filename = os.path.join(spectrum_path, "SPECTRA_" + kvp + "kVp_" + \
                "17deg_" + filt + ".txt")

                src_energies, photons = read_spectrum(spectrum_filename)
                # Normalization
                maxphotons = np.sum(photons)
                for i in range(len(photons)): photons[i] = 1.*photons[i]/maxphotons

                photons_tailed = []
                
                for eng in energies:
                    cnt = 0
                    if eng not in src_energies:
                        photons_tailed.append(0)
                    else:
                        i = src_energies.index(eng)
                        photons_tailed.append(photons[i])

                photons_tailed = np.cumsum(photons_tailed)

                spectra_enum_names.append('kVp_' + kvp + "_" + filt)
                spectra_pointers.append("&" + "spectrum_" + kvp + "kVp_" + filt)
                source_s = get_list(photons_tailed, "spectrum_" + kvp + "kVp_" + filt, "double", "Photon count per energy value")
                of.write(source_s)
                of.write("\n")

        spectra_enum = "enum spectra { "
        for key in spectra_enum_names:
            spectra_enum += key.upper() + ", "
        spectra_enum = spectra_enum[:-2]
        spectra_enum += " };"

        spectra_map = get_map(np.arange(0, len(spectra_pointers)), spectra_pointers, "spectra_keys", "int", "const std::vector<double>*", "Map of all spectra")
        of.write(spectra_map)
        of.write("\n\n")

        keys_pointers = []
        materials_enum = "enum materials { "
        for key in materials.keys():
            keys_pointers.append("&" + key)
            attenuation = get_list(materials[key], key, "double", key + " material")
            of.write(attenuation)
            of.write("\n")

            materials_enum += key.upper() + ", "
        materials_enum += "VACUUM = -1 };"

        keys = get_map(np.arange(0, len(materials.keys())), keys_pointers, "materials_keys", "int", "const std::vector<double>*", "Map of all materials")
        of.write(keys)

        densities = [1.050, 9.500e-01, 1.050, 1.040, 1.070, 1.040, 1.920, 1.060, 1.050, 1.020, 1]
        densities_str = get_list(densities, "material_densities", "double", "Densities of materials")
        of.write("\n" + densities_str)
        of.write("\n\n")

        effs_pointers = []
        effs_enum = "enum efficiencies { "
        for key in det_effs.keys():
            effs_pointers.append("&" + "coat_" + key)
            efficiency = get_list(det_effs[key], "coat_" + key, "double", key + " efficiency")
            of.write(efficiency)
            of.write("\n")

            effs_enum += key.upper() + ", "
        effs_enum = effs_enum[:-2]
        effs_enum += " };"

        keys = get_map(np.arange(0, len(det_effs.keys())), effs_pointers, "efficiency_keys", "int", "const std::vector<double>*", "Map of all detector coating efficiencies")
        of.write(keys)
        of.write("\n\n")

        of.write("\n\t" + spectra_enum + "\n")

        of.write("\n\t" + materials_enum + "\n")

        of.write("\n\t" + effs_enum + "\n")
        of.write("\n")


        of.write("} // namespace\n\n#endif")
    
    
    