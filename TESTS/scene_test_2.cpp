#include "XRaySimulator.hpp"

int main()
{
    xrp::Source* s = xrp::FiniteSource::create_source(xru::Vector3D(0,-100,0), xru::Vector3D(0, 1, 0), xru::Vector3D(0, 0, 1), 0.2, 0.2, xrc::KVP_60_1AL);
    xrd::Detector* d = new xrd::Detector(xru::Vector3D(0, 7, 0), xru::Vector3D(0, -1, 0), xru::Vector3D(0, 0, 1), 380, 380, 0.1, 0.1);

    xrs::Scene* scene = new xrs::Scene(s, d);
    scene->detector_->add_coating(0.01, xrc::SE);

    xrg::Body* ellipsoid = new xrg::Ellipsoid(30, 6, 30, 0, 0, 0);
    ellipsoid->set_material(xrc::materials::WATER);
    scene->add_body(ellipsoid);

    xrg::Body* sph1 = new xrg::Sphere(2, -12, -1, 5);
    sph1->set_material(xrc::materials::BLOOD);
    scene->add_body(sph1);

    xrg::Body* sph2 = new xrg::Sphere(3, -6, 0, -10);
    sph2->set_material(xrc::materials::MUSCLE);
    scene->add_body(sph2);

    xrg::Body* sph3 = new xrg::Sphere(3,  5,  1, 5);
    sph3->set_material(xrc::materials::BONE);
    scene->add_body(sph3);

    xrg::Body* cyl  = new xrg::Cylinder(2, 8, 0, 0, 0);
    cyl->set_material(xrc::materials::LUNG);
    scene->add_body(cyl);

    scene->shoot_photons(1e6);

    for (auto body: scene->bodies_)
        std::cout << *body << " received the dose of " << body->calculate_dose() << " Gy." << std::endl;

    scene->detector_->save_to_file("../../HELPERS/RESULTS.txt");

    delete scene;
}