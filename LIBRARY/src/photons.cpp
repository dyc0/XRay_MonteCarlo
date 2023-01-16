#include "XRaySimulator.hpp"
#include "photons.hpp"

namespace xrp
{

    Photon::Photon(const xru::Point3D &origin, const xru::Vector3D &direction, const float energy, const int lerp_energy_index, const float lerp_percentage):
        origin_(origin), direction_(direction), energy_(energy), lerp_energy_index_(lerp_energy_index), lerp_percentage_(lerp_percentage)
        { }
    
    std::ostream &operator<<(std::ostream &out, const Photon &p)
    {
        out << p.energy_ << " " << p.origin_ << " " << p.direction_;
        return out;
    }

    Source *Source::source_ = nullptr;

    Source::Source(const xru::Point3D& origin, const int spectrum) : origin_(origin)
    {
        rng_ = xru::RandomGenerator::initialize_random_generator();
        spectrum_ = xrc::spectra_keys.at(spectrum);
    }

    inline xru::Vector3D Source::generate_random_vector() const
    {
        // Found this here: https://math.stackexchange.com/questions/44689/how-to-find-a-random-axis-or-unit-vector-in-3d
        // Seems to work faster than generating a normal distribution.
    
        double theta = xru::RandomGenerator::random_range(0, 2*xrc::PI);
        double z     = xru::RandomGenerator::random_scalar();
        double z2    = z*z;

        // No need for norming, as the vector already is normed.
        return xru::Vector3D(sqrt(1-z2)*cos(theta), sqrt(1-z2)*sin(theta), z);
    }

    void xrp::Source::generate_energy(int& index, double& lerp_perc, float& energy) const
    {
        // CHECK: Edge cases when equal to an energy value.
        // The range is to avoid edge cases and branching:
        double prob = xru::RandomGenerator::random_range(xrc::tolerance, 1-xrc::tolerance);
        
        index = xru::index_finder(prob, *spectrum_);
        lerp_perc = (prob - (*spectrum_)[index-1]) / ((*spectrum_)[index] - (*spectrum_)[index - 1]);
        energy = (xrc::energies[index] - xrc::energies[index-1]) * lerp_perc + xrc::energies[index-1];
    }

    xrp::Source::~Source()
    {
        Source::source_ = nullptr;
        delete rng_;
    }

    void Source::set_spectrum(const unsigned int spectrum)
    {
        spectrum_ = xrc::spectra_keys.at(spectrum);
    }

    PointSource::PointSource(const xru::Point3D& origin, const int spectrum): Source(origin, spectrum) {}

    PointSource* PointSource::create_source(const xru::Point3D& origin, const int spectrum)
    {
        if (source_ == nullptr)
            source_ = new PointSource(origin, spectrum);
        else
            std::cout << "Only one source can be created per scene. " << std::endl << std::endl;

        return static_cast<PointSource*> (source_);
    }

    Photon *PointSource::generate_particle() const
    {
        int index;
        double lerp_percentage;
        float energy;

        xru::Vector3D direction = generate_random_vector();
        generate_energy(index, lerp_percentage, energy);
        return new Photon(origin_, direction, energy, index, lerp_percentage);
    }

    FiniteSource::FiniteSource(const xru::Point3D& origin, const xru::Vector3D& normal, const xru::Vector3D& local_y, const float dx, const float dy, const int spectrum):
        Source(origin, spectrum), normal_(normal), local_y_(local_y), dx_(dx), dy_(dy) { };

    FiniteSource* FiniteSource::create_source(const xru::Point3D& origin, const xru::Vector3D& normal, const xru::Vector3D& local_y, const float dx, const float dy, const int spectrum)
    {
        if (source_ == nullptr)
            source_ = new FiniteSource(origin, normal, local_y, dx, dy, spectrum);
        else
            std::cout << "Only one source can be created per scene. " << std::endl << std::endl;

        return static_cast<FiniteSource*>(source_);
    }

    Photon *FiniteSource::generate_particle() const
    {
        xru::Vector3D direction = generate_random_vector();
        xru::Vector3D ypos = local_y_ * xru::RandomGenerator::random_scalar() * dy_;
        xru::Vector3D xpos = local_y_.cross(normal_) * xru::RandomGenerator::random_scalar() * dx_;

        int index;
        double lerp_percentage;
        float energy;
        generate_energy(index, lerp_percentage, energy);

        return new Photon(origin_ + xpos + ypos, direction, energy, index, lerp_percentage);
    }
}