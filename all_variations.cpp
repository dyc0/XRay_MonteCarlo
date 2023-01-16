#include "XRaySimulator.hpp"

int main(int argc, char* argv[])
{
    xrp::Source* s = xrp::FiniteSource::create_source(xru::Vector3D(0,-100,0), xru::Vector3D(0, 1, 0), xru::Vector3D(0, 0, 1), 0.1, 0.2, 0);
    xrd::Detector* d = new xrd::Detector(xru::Vector3D(0, 7, 0), xru::Vector3D(0, -1, 0), xru::Vector3D(0, 0, 1), 380, 380, 0.1, 0.1);

    // CREATING SCENE
    xrs::Scene* scene = new xrs::Scene(s, d);

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
    

    std::string base_path = "../RESULTS/";
    std::string current_filename;
    for (unsigned int spectrum = 0; spectrum < 4; spectrum++)
    {
        scene->source_->set_spectrum(spectrum);
        for (int efficiency = -1; efficiency < 2; efficiency++)
        {
            if (efficiency >= 0) scene->detector_->add_coating(2e-4, efficiency);
            else scene->detector_->coating_material_ = nullptr;

            scene->shoot_photons(1e7);

            current_filename = std::to_string(spectrum) + std::to_string(efficiency);
            scene->detector_->save_to_file(base_path + std::string("image_") + current_filename + std::string(".txt"));
            scene->save_absorptions(base_path + std::string("doses/doses_") + current_filename + std::string(".txt"));

            for (auto body: scene->bodies_)
                body->absorbed_energy_ = 0;
        }
    }

    delete scene;
}