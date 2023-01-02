#include "XRaySimulator.hpp"

using namespace xrp;

int main(int argc, char* argv[])
{

    xru::Point3D origin = xru::Point3D();
    Photon* photon, * ph;

    std::cout << "\n\n\t------------" << "\n\tPOINT SOURCE" << "\n\t------------\n\n";

    PointSource* ps = PointSource::create_source(origin, xrc::KVP_120_1AL);

    photon = ps->generate_particle();
    if (abs(photon->direction_.dot(photon->direction_) - 1) > xrc::tolerance) return 1;
    else if (photon->energy_ < xrc::energies.front() || photon->energy_ > xrc::energies.back()) return 2;
    std::cout << "Photon (E, ORIGIN, DIRECTION): " << *photon << std::endl;
    delete photon;

    // For checking spacial and energy distributions:

    std::cout << "Generating a bunch of 'em..." << std::endl;

    std::ofstream out("../../HELPERS/ph_vis_data/photons_point.txt");
    for (int i = 0; i < 1e6; i++)
    {
        ph = ps->generate_particle();
        out << *ph << std::endl;
        delete ph;
        if (i % 1000 == 0)
            std::cout << "\r" << int(100. * i * 1e-6) << "\%" << std::flush;
    }
    std::cout << "\rDONE." << std::endl;
    out.close();

    delete ps;



    std::cout << "\n\n\t-------------" << "\n\tFINITE SOURCE" << "\n\t-------------\n\n";

    xru::Vector3D normal = xru::Vector3D(1, 0, 0);
    xru::Vector3D y_axis = xru::Vector3D(0, 1, 0);

    FiniteSource* fs = FiniteSource::create_source(origin, normal, y_axis, 10, 5, xrc::KVP_120_1AL);

    photon = fs->generate_particle();
    if (abs(photon->direction_.dot(photon->direction_) - 1) > xrc::tolerance) return 1;
    else if (photon->energy_ < xrc::energies.front() || photon->energy_ > xrc::energies.back()) return 2;
    std::cout << "Photon (E, ORIGIN, DIRECTION): " << *photon << std::endl;
    delete photon;

    // For checking spacial and energy distributions:

    std::cout << "Generating a bunch of 'em..." << std::endl;

    std::ofstream out2("../../HELPERS/ph_vis_data/photons_finite.txt");
    for (int i = 0; i < 1e6; i++)
    {
        ph = fs->generate_particle();
        out2 << *ph << std::endl;
        delete ph;
        if (i % 1000 == 0)
            std::cout << "\r" << int(100. * i * 1e-6) << "\%" << std::flush;
    }
    std::cout << "\rDONE." << std::endl;
    out2.close();

    std::ofstream out3("../../HELPERS/ph_vis_data/finite_dims.txt");
    out3 << fs->origin_ << " " << fs->normal_ << " " << fs->local_y_ << " " << fs->dx_ << " " << fs->dy_ << std::endl;
    out3.close();

    delete fs;
}