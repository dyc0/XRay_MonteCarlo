#ifndef DETECTOR_HPP
#define DETECTOR_HPP

namespace xrd
{
    
    struct Pixel
    {
        unsigned int photons: 12;
    };

    typedef std::vector<std::vector<Pixel*>> PixelMatrix;

    class Detector
    {
        public:
        Detector(const xru::Point3D& centre, const xru::Vector3D& normal, const xru::Vector3D local_y,
                 const size_t npx_x, const size_t npx_y, const double dx, const double dy);
        ~Detector();

        void check_hit(const xrp::Photon& photon, size_t* intersection_index, bool& hit);
        void increment_pixel(const size_t x, const size_t y);

        void add_coating(const double width, const int material);
        bool coating_hit(const xrp::Photon* photon, const double running_sum, const double treshold);

        void save_to_file(const std::string& path) const;

        xru::Point3D centre_;
        xru::Vector3D normal_, local_y_;
        size_t npx_x_, npx_y_;
        double dx_, dy_;    // Could be that this is not used
        double xd_, yd_;
        double x_span_, y_span_;

        xru::Vector3D coating_layer_;
        const std::vector<double>* coating_material_;

        PixelMatrix m_;
    };

}


#endif