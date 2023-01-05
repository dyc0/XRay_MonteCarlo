#include "XRaySimulator.hpp"
#include "detector.hpp"

xrd::Detector::Detector(const xru::Point3D &centre, const xru::Vector3D &normal, const xru::Vector3D local_y,
                        const size_t npx_x, const size_t npx_y, const double dx, const double dy):
    centre_(centre), normal_(normal), local_y_(local_y), npx_x_(npx_x), npx_y_(npx_y), dx_(dx), dy_(dy), xd_(1./dx), yd_(1./dy), coating_material_(nullptr)
{
    normal_.norm();
    local_y_.norm();

    x_span_ = npx_x * dx * 0.5;
    y_span_ = npx_y * dy * 0.5;

    // Populate pixels
    for (size_t px_row = 0; px_row < npx_y; px_row++)
    {
        m_.push_back(std::vector<Pixel*>());
        for(size_t px = 0; px < npx_x; px++)
            m_[px_row].push_back(new Pixel());
    }
}

xrd::Detector::~Detector()
{
    for (auto row: m_)
        for (auto px: row)
            delete px;
}

void xrd::Detector::check_hit(const xrp::Photon &photon, size_t *intersection_index, bool &hit)
{
    // CHECK: This code might be inefficient.

    hit = false;

    // First check plane intersection - the normal should be directed to origin
    if (photon.direction_.dot(centre_) < xrc::tolerance) [[likely]]     // -> likely because this is used in photon generation
        return;

    double t = (centre_ - photon.origin_).dot(normal_) / photon.direction_.dot(normal_);
    // Check if it's behind
    hit = (t > xrc::tolerance);

    // Get local coordinates
    xru::Vector3D local_vec = ( photon.origin_ + photon.direction_ * t - centre_ ) * hit;
    double local_x = local_vec.dot(local_y_.cross(normal_));
    double local_y = local_vec.dot(local_y_);

    // Check if it is hit and not parallel to detector and within the detector's bounds
    hit = (  hit && abs(photon.direction_.dot(normal_)) > xrc::tolerance
          && local_x < x_span_ - xrc::tolerance && local_x > -x_span_ + xrc::tolerance
          && local_y < y_span_ - xrc::tolerance && local_y > -y_span_ + xrc::tolerance );

    intersection_index[0] = static_cast<size_t>((y_span_ + local_y) * yd_) * hit;
    intersection_index[1] = static_cast<size_t>((x_span_ + local_x) * xd_) * hit;

    if (hit)
    {
        assert(intersection_index[0] < npx_y_);
        assert(intersection_index[1] < npx_x_);
    }
}

void xrd::Detector::increment_pixel(const size_t x, const size_t y)
{
    m_[y][x]->photons++;
}

void xrd::Detector::add_coating(const double width, const int material)
{
    coating_layer_ = normal_* (-width);
    coating_material_ = xrc::efficiency_keys.at(material);
}

bool xrd::Detector::coating_hit(const xrp::Photon *photon, const double running_sum, const double treshold)
{
    if (coating_material_ == nullptr) return true;

    double attenuation_coef = photon->lerp_percentage_ * coating_material_->at(photon->lerp_energy_index_) +
                             (1 - photon->lerp_percentage_) * coating_material_->at(photon->lerp_energy_index_ - 1);

    return (attenuation_coef * coating_layer_.dot(coating_layer_) / (coating_layer_.dot(photon->direction_) + running_sum) < treshold);
}

void xrd::Detector::save_to_file(const std::string &path) const
{
    std::ofstream output_file;
    try { output_file.open(path, std::ios::out); }
    catch (...) { std::cerr << "Error opening file " << path << "." << std::endl; }

    std::string row_data;
    for (size_t row = 0; row < npx_y_; row++)
    {   
        row_data = "";
        for (size_t column = 0; column < npx_x_; column++)
            row_data.append(std::to_string(m_[row][column]->photons)).append(" ");
        row_data.pop_back();
        output_file << row_data << std::endl;
    }

    output_file.close();
}
