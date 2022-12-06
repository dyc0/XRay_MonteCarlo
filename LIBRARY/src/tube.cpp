#include "XRaySimulator.hpp"

using namespace xrt;

XRay::XRay(const xru::Vector3D& direction):
    direction_(direction)
{  
    spectrum = xrc::spectrum;
}

static void set_source(const xru::Point3D& origin)
{
    XRay::source_ = origin;
}

std::vector<XRay*> XRay::generate_rays(const Detector& detector)
{
    assert(!detector.pixels.empty());

    std::vector<XRay*> rays;

    for (auto px: detector.pixels)
        rays.push_back(new XRay((px->center_position - XRay::source_).normed()));

    return rays;
}


Pixel::Pixel(const xru::Point3D& center): center_position(center) { };


Detector::Detector(const xru::Point3D& center, const xru::Vector3D& face_direction, const xru::Vector3D& local_y_axis):
    center_(center), face_(face_direction.normed()), y_(local_y_axis.normed()) { };

void Detector::populate_pixels(const int x_number, const int y_number, const double px_width, const double px_height)
{
    assert((x_number % 2 == 0) && (y_number % 2 == 0));

    xru::Vector3D x_ = y_.cross(face_);
    x_.norm();
    
    for (int dx = -x_number/2; dx <= x_number/2; dx++)
        for (int dy = -y_number/2; dy <= y_number/2; dy++)
        {
            pixels.push_back(new Pixel(xru::Vector3D(x_ * px_width/2 * dx) + xru::Vector3D(y_ * px_height/2 * dy)));
            photon_count.push_back(0);
        }
}