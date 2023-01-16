#include "XRaySimulator.hpp"

int main()
{
    std::cout << "Generating a bunch of 'em..." << std::endl;

    xru::Vector3D normal = xru::Vector3D(1, 0, 0);
    xru::Vector3D y_axis = xru::Vector3D(0, 1, 0);
    xru::Vector3D origin = xru::Vector3D(0, 0, 0);

    xrp::FiniteSource* fs = xrp::FiniteSource::create_source(origin, normal, y_axis, 10, 5, 0);
    xrp::Photon* ph;

    std::string base_path = "../RESULTS/spectra/";
    std::string filename;
    for (unsigned int spectrum=0; spectrum<4; spectrum++)
    {
        fs->set_spectrum(spectrum);
        filename = "spectrum_" + std::to_string(spectrum) + ".txt";
        std::ofstream out(base_path + filename);
        for (int i = 0; i < 1e6; i++)
        {
            ph = fs->generate_particle();
            out << ph->energy_ << " ";
            delete ph;
        }
        std::cout << "DONE." << std::endl;
        out.close();
    }
}