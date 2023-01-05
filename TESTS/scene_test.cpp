#include "XRaySimulator.hpp"

using namespace xrs;

void test_photon(Scene* scene, xrp::Photon* ph, bool exp_hit, size_t exp_indices[2], int err)
{
    size_t indices[2];
    bool hit = false;

    std::cout << "-------<USING RAY:" << ph->direction_ << ">-------" << std::endl;
    scene->detector_->check_hit(*ph, indices, hit);
    std::cout << "Got results: { HIT: " << hit << ", COORDINATES: (" << indices[0] << ", " << indices[1] << ") }" << std::endl;
    if (exp_hit != hit) throw err;
    if (hit)
    {
        if (exp_indices[0] != indices[0] || exp_indices[1] != indices[1]) throw err;
        scene->detector_->increment_pixel(indices[0], indices[1]);
    }
}

int main()
{
    xrp::Source* s = xrp::PointSource::create_source(xru::Vector3D(0,0,0), xrc::KVP_60_1AL);
    xrd::Detector* d = new xrd::Detector(xru::Vector3D(1, 0, 0), xru::Vector3D(-1, 0, 0), xru::Vector3D(0, 0, 1), 10, 10, 1, 1);

    Scene* scene = new Scene(s, d);
    
    xrp::Photon* ph = new xrp::Photon(xru::Vector3D(0,0,0), xru::Vector3D(1, 0.5, 0.5), 1, 1, 1);
    size_t exp_indices[2];

    try
    {
        exp_indices[0] = 5;
        exp_indices[1] = 4;
        test_photon(scene, ph, true, exp_indices, 1);
        std::cout << std::endl;

        ph->direction_ = xru::Vector3D(1, 10, 10);
        test_photon(scene, ph, false, exp_indices, 2);
        std::cout << std::endl;

        ph->direction_ = xru::Vector3D(-1, 1, 1);
        test_photon(scene, ph, false, exp_indices, 3);
        std::cout << std::endl;

        exp_indices[0] = 6;
        exp_indices[1] = 4;
        ph->direction_ = xru::Vector3D(1, 1, 1);
        test_photon(scene, ph, true, exp_indices, 4);
        std::cout << std::endl;
    }
    catch (int err)
    {
        std::cout << err << std::endl;
        return err;
    }
    

    for (auto row: scene->detector_->m_)
    {
        for (auto px: row)
        {
            std::cout << px->photons << " ";
        }
        std::cout << std::endl;
    }
}